#!/usr/bin/env python3
"""Read-only candidate protein/coordinate diagnostic. No tree input or locus deletion."""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import math
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

EXPECTED_PROTEOME = "4daa19724e54e2e657c4578ac97c5c13997c369ca78d3f2b296b3b3941d95f90"
EXPECTED_GTF = "d9724576f7521700c5ea0770a22792a055cbe93ebe9a918edb046cacb1de6203"
DT = "Dioscorea_tenuipes_historical"
FIELDS = "query target qlen slen identity aligned mismatches gapopens qstart qend sstart send evalue bitscore qseq sseq".split()
INTS = "qlen slen aligned mismatches gapopens qstart qend sstart send".split()
TIERS = (0, 50, 80, 95, 99, 100)
CODONS = {a+b+c: aa for (a,b,c), aa in zip(((a,b,c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
           "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
COMPLEMENT = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def fasta(path):
    records = {}
    name = None
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                assert name not in records, ("duplicate FASTA ID", name)
                records[name] = ""
            elif line.strip():
                assert name is not None
                records[name] += line.strip().upper()
    return records


def write_fasta(path, records):
    with path.open("w") as stream:
        for ident, sequence in records.items():
            stream.write(">"+ident+"\n")
            stream.write("\n".join(sequence[i:i+80] for i in range(0, len(sequence), 80))+"\n")


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fallback=("status", "reason")):
    keys = list(dict.fromkeys(k for r in rows for k in r)) or list(fallback)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, keys, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def run_command(command, receipt, stdout_path=None):
    before = time.monotonic()
    if stdout_path is None:
        result = subprocess.run(command, capture_output=True)
        # Version/build metadata can contain non-UTF8 compiler-path bytes.
        # Escape those bytes visibly; never apply this decoding to scientific files.
        stdout = result.stdout.decode("utf-8", errors="backslashreplace")
    else:
        with stdout_path.open("wb") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.PIPE)
        stdout = None
    stderr = result.stderr.decode("utf-8", errors="backslashreplace")
    receipt["commands"].append({"argv": list(map(str, command)), "exit_code": result.returncode,
                               "elapsed_seconds": time.monotonic()-before,
                               "stderr": stderr, "stderr_raw_sha256": hashlib.sha256(result.stderr).hexdigest(),
                               "stdout_raw_sha256": hashlib.sha256(result.stdout).hexdigest() if stdout_path is None else None,
                               "metadata_decode": "utf-8 with visible backslash escaping of invalid bytes",
                               "scientific_stdout_preserved_as_bytes": stdout_path is not None,
                               "stdout_path": str(stdout_path) if stdout_path else None})
    if result.returncode:
        raise RuntimeError("Command failed: "+str(command)+"\n"+stderr)
    return stdout


def parse_blast(path, queries, subjects):
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, fieldnames=FIELDS, delimiter="\t")
        for row in reader:
            for field in INTS:
                row[field] = int(row[field])
            for field in ("identity", "evalue", "bitscore"):
                row[field] = float(row[field])
            assert row["query"] in queries and row["target"] in subjects
            assert row["qlen"] == len(queries[row["query"]]) and row["slen"] == len(subjects[row["target"]])
            q, s = row["qseq"].upper(), row["sseq"].upper()
            assert len(q) == len(s) == row["aligned"]
            assert q.replace("-", "") == queries[row["query"]][row["qstart"]-1:row["qend"]]
            assert s.replace("-", "") == subjects[row["target"]][row["sstart"]-1:row["send"]]
            assert not any(a == b == "-" for a, b in zip(q, s))
            identical = sum(a == b for a, b in zip(q, s))
            mismatches = sum(a != b and a != "-" and b != "-" for a, b in zip(q, s))
            gaps = sum(sequence[i] == "-" and (i == 0 or sequence[i-1] != "-")
                       for sequence in (q, s) for i in range(len(sequence)))
            assert abs(100*identical/len(q)-row["identity"]) <= .00051
            assert mismatches == row["mismatches"] and gaps == row["gapopens"]
            row["identity_recomputed"] = 100*identical/len(q)
            rows.append(row)
    return rows


def paired_positions(row):
    qpos, spos = row["qstart"], row["sstart"]
    pairs = []
    for q, s in zip(row["qseq"], row["sseq"]):
        if q != "-" and s != "-":
            pairs.append((qpos, spos, q == s))
        qpos += q != "-"
        spos += s != "-"
    assert qpos == row["qend"]+1 and spos == row["send"]+1
    return pairs


def chain_coverage(rows):
    ordered = sorted(rows, key=lambda r: (r["qend"], r["send"], r["qstart"], r["sstart"]))
    states = []
    for i, row in enumerate(ordered):
        weight = len(paired_positions(row))
        eligible = [states[j] for j, previous in enumerate(ordered[:i])
                    if previous["qend"] < row["qstart"] and previous["send"] < row["sstart"]]
        best = max(eligible, key=lambda x: x[0], default=(0, []))
        states.append((best[0]+weight, best[1]+[i]))
    weight, indices = max(states, key=lambda x: x[0], default=(0, []))
    return weight, [[ordered[i][k] for k in ("qstart", "qend", "sstart", "send")] for i in indices]


def coverage_tables(rows, queries):
    target_rows, query_rows, repeats = [], [], []
    by_query = collections.defaultdict(list)
    for row in rows:
        by_query[row["query"]].append(row)
    for query, sequence in queries.items():
        for tier in TIERS:
            selected = [r for r in by_query[query] if r["aligned"] >= 100 and r["identity_recomputed"] >= tier]
            by_target = collections.defaultdict(list)
            for row in selected:
                by_target[row["target"]].append(row)
            pooled = set()
            local = []
            for target, hits in sorted(by_target.items()):
                qcovered, scovered = set(), set()
                spos_queries = collections.defaultdict(set)
                for hit in hits:
                    for qp, sp, _ in paired_positions(hit):
                        qcovered.add(qp); scovered.add(sp); spos_queries[sp].add(qp)
                reused = {sp: qs for sp, qs in spos_queries.items() if len(qs) > 1}
                chain, chain_intervals = chain_coverage(hits)
                row = {"query": query, "target": target, "minimum_hsp_identity_percent": tier,
                       "minimum_hsp_alignment_length": 100, "qualified_hsps": len(hits),
                       "query_length": len(sequence), "subject_length": hits[0]["slen"],
                       "paired_query_positions_union": len(qcovered), "query_union_fraction": len(qcovered)/len(sequence),
                       "paired_subject_positions_union": len(scovered), "subject_union_fraction": len(scovered)/hits[0]["slen"],
                       "subject_positions_reused_by_distinct_query_positions": len(reused),
                       "max_distinct_query_positions_per_subject_position": max(map(len, spos_queries.values()), default=0),
                       "conservative_collinear_chain_query_positions": chain,
                       "conservative_collinear_chain_query_fraction": chain/len(sequence),
                       "chain_intervals_json": json.dumps(chain_intervals), "maximum_hsp_bitscore": max(r["bitscore"] for r in hits)}
                target_rows.append(row); local.append(row); pooled.update(qcovered)
                if tier == 80 and reused:
                    repeats.append({"query": query, "target": target, "hsp_identity_tier": tier,
                                    "reused_subject_positions": len(reused),
                                    "max_distinct_query_positions_per_subject_position": row["max_distinct_query_positions_per_subject_position"],
                                    "examples_first_10_json": json.dumps({str(s): sorted(qs) for s, qs in sorted(reused.items())[:10]})})
            best_union = max(local, key=lambda r: (r["paired_query_positions_union"], r["maximum_hsp_bitscore"]), default=None)
            best_chain = max(local, key=lambda r: (r["conservative_collinear_chain_query_positions"], r["maximum_hsp_bitscore"]), default=None)
            union_fraction = best_union["query_union_fraction"] if best_union else 0
            chain_fraction = best_chain["conservative_collinear_chain_query_fraction"] if best_chain else 0
            pooled_fraction = len(pooled)/len(sequence)
            label = ("AT_LEAST_ONE_SINGLE_TARGET_CHAIN_GE90_PERCENT" if chain_fraction >= .9 else
                     "SINGLE_TARGET_UNION_GE90_PERCENT_WITHOUT_GE90_CHAIN" if union_fraction >= .9 else
                     "ONLY_MULTIPLE_TARGETS_POOLED_GE90_PERCENT" if pooled_fraction >= .9 else
                     "PARTIAL_LOCAL_COVERAGE" if local else "NO_HSP_MEETS_DESCRIPTIVE_TIER")
            query_rows.append({"query": query, "minimum_hsp_identity_percent": tier, "qualified_targets": len(local),
                               "qualified_hsps": len(selected), "query_length": len(sequence),
                               "best_union_target": best_union["target"] if best_union else None,
                               "best_single_target_union_fraction": union_fraction,
                               "best_union_target_reused_subject_positions": best_union["subject_positions_reused_by_distinct_query_positions"] if best_union else None,
                               "best_chain_target": best_chain["target"] if best_chain else None,
                               "best_single_target_chain_fraction": chain_fraction,
                               "all_targets_pooled_query_positions": len(pooled), "all_targets_pooled_query_fraction": pooled_fraction,
                               "descriptive_coverage_label": label})
    return target_rows, query_rows, repeats


def annotation_records(path, wanted):
    models = collections.defaultdict(list)
    selected_lines = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            identifiers = set(re.findall(r'(?:transcript_id|protein_id) "([^"]+)"', f[8]))
            for attribute in f[8].split(";"):
                if "=" in attribute:
                    k, v = attribute.strip().split("=", 1)
                    if k in ("ID", "Parent", "protein_id", "transcript_id"):
                        identifiers.update(v.split(","))
            matches = identifiers & wanted
            if not matches:
                continue
            selected_lines.append(line.rstrip("\n"))
            if f[2] == "CDS":
                for ident in matches:
                    models[ident].append({"contig": f[0], "start": int(f[3]), "end": int(f[4]),
                                          "strand": f[6], "phase": f[7]})
    return dict(models), selected_lines


def translate_model(features, exon_sequences, expected):
    result = {"cds_segments": len(features), "status": "NOT_VERIFIED"}
    try:
        contigs = {f["contig"] for f in features}; strands = {f["strand"] for f in features}
        assert len(contigs) == len(strands) == 1 and next(iter(strands)) in ("+", "-")
        strand = next(iter(strands)); contig = next(iter(contigs))
        ordered = sorted(features, key=lambda f: f["start"], reverse=strand == "-")
        genomic = sorted(features, key=lambda f: f["start"])
        assert all(a["end"] < b["start"] for a, b in zip(genomic, genomic[1:])), "overlapping_CDS_records"
        phases = [int(f["phase"]) for f in ordered]
        assert all(p in (0, 1, 2) for p in phases)
        sequence = ""; phase_checks = []
        for i, feature in enumerate(ordered):
            key = f'{contig}:{feature["start"]}-{feature["end"]}'
            dna = exon_sequences[key]
            assert len(dna) == feature["end"]-feature["start"]+1
            if strand == "-":
                dna = dna.translate(COMPLEMENT)[::-1]
            if i == 0:
                dna = dna[phases[0]:]
            else:
                phase_checks.append(phases[i] == (-len(sequence)) % 3)
            sequence += dna
        protein = "".join(CODONS.get(sequence[i:i+3], "X") for i in range(0, len(sequence)-2, 3))
        exact = protein.rstrip("*") == expected.rstrip("*")
        internal = "*" in protein.rstrip("*")
        result.update({"contig": contig, "strand": strand, "start": min(f["start"] for f in features),
                       "end": max(f["end"] for f in features), "initial_phase": phases[0], "coding_nt": len(sequence),
                       "partial_codon_nt": len(sequence) % 3, "phase_continuity": all(phase_checks),
                       "translated_length_without_terminal_stop": len(protein.rstrip("*")),
                       "expected_length_without_terminal_stop": len(expected.rstrip("*")),
                       "internal_stop": internal, "exact_translation_terminal_stop_normalized": exact,
                       "coding_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                       "translated_sha256": hashlib.sha256(protein.encode()).hexdigest(),
                       "status": "EXACT_COORDINATE_TRANSLATION" if exact and not internal and all(phase_checks) and len(sequence)%3 == 0 else "TRANSLATION_OR_PHASE_DISAGREEMENT"})
    except (AssertionError, ValueError, KeyError) as error:
        result["reason"] = type(error).__name__+":"+str(error)
    return result


def fixtures():
    h = {"qseq": "AB-CD", "sseq": "A-BCD", "qstart": 1, "qend": 4, "sstart": 1, "send": 4}
    assert paired_positions(h) == [(1,1,True), (3,3,True), (4,4,True)]
    def hit(q, s, n):
        return {"qseq": "A"*n, "sseq": "A"*n, "qstart": q, "qend": q+n-1, "sstart": s, "send": s+n-1}
    assert chain_coverage([hit(1,1,10), hit(11,1,10)])[0] == 10
    assert chain_coverage([hit(1,1,10), hit(11,11,10)])[0] == 20
    plus = [{"contig":"c", "start":1,"end":4,"strand":"+","phase":"0"},
            {"contig":"c", "start":10,"end":14,"strand":"+","phase":"2"}]
    assert translate_model(plus, {"c:1-4":"ATGG", "c:10-14":"CCTAA"}, "MA")["status"] == "EXACT_COORDINATE_TRANSLATION"
    minus = [{"contig":"c", "start":1,"end":5,"strand":"-","phase":"2"},
             {"contig":"c", "start":10,"end":13,"strand":"-","phase":"0"}]
    assert translate_model(minus, {"c:1-5":"TTAGG", "c:10-13":"CCAT"}, "MA*")["status"] == "EXACT_COORDINATE_TRANSLATION"
    wrong = [dict(f) for f in plus]; wrong[1]["phase"]="1"
    assert translate_model(wrong, {"c:1-4":"ATGG", "c:10-14":"CCTAA"}, "MA")["status"] == "TRANSLATION_OR_PHASE_DISAGREEMENT"
    command_receipt = {"commands": []}
    version_bytes = b"samtools fixture\ncompiler-path: \xff\xfe\n"
    version = run_command([sys.executable, "-c", "import sys;sys.stdout.buffer.write("+repr(version_bytes)+");sys.stderr.buffer.write(b'warn: \\x80\\n')"], command_receipt)
    assert version == version_bytes.decode("utf-8", errors="backslashreplace")
    assert "\\xff\\xfe" in version and "\\x80" in command_receipt["commands"][-1]["stderr"]
    assert command_receipt["commands"][-1]["stdout_raw_sha256"] == hashlib.sha256(version_bytes).hexdigest()
    scientific_bytes = b">fixture\r\nACGT\r\n"
    with tempfile.TemporaryDirectory(prefix="dt-structure-bytes-") as temporary:
        output = Path(temporary) / "reference.fna"
        run_command([sys.executable, "-c", "import sys;sys.stdout.buffer.write("+repr(scientific_bytes)+")"], command_receipt, output)
        assert output.read_bytes() == scientific_bytes
    return {"paired_positions_through_gaps": "PASS", "reused_vs_collinear_subject": "PASS",
            "plus_minus_translation_phase_and_terminal_stop": "PASS",
            "non_utf8_version_metadata": "PASS", "scientific_stdout_bytes_preserved": "PASS"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-directory", type=Path, default=Path("workspace/wgs/dt_research_validation_20260924/global_qc/screen_v1"))
    parser.add_argument("--candidate-evidence", type=Path)
    parser.add_argument("--wgs-root", type=Path, default=Path("workspace/wgs"))
    parser.add_argument("--blastp", default="workspace/miniconda3/envs/busco/bin/blastp")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--historical-annotation", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    fixture_checks = fixtures()
    if args.self_test:
        print(json.dumps(fixture_checks, indent=2)); return
    if args.output is None:
        parser.error("--output is required")
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output; start = time.monotonic()
    receipt = {"state": "RUNNING", "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "command_argv": sys.argv, "python_version": sys.version, "commands": [], "inputs_sha256": {},
               "fixture_checks": fixture_checks}
    def bind(path):
        receipt["inputs_sha256"][str(path)] = sha(path)
        return receipt["inputs_sha256"][str(path)]
    try:
        bind(Path(__file__).resolve()); bind(Path(__file__).with_name("DESIGN.md"))
        previous = json.loads((args.screen_directory / "RECEIPT.json").read_text())
        assert previous["state"] == "COMPLETE"
        screen_path = args.screen_directory / "ALL_CELLS.tsv"
        assert bind(screen_path) == previous["output_sha256"]["ALL_CELLS.tsv"]
        bind(args.screen_directory / "RECEIPT.json")
        observations = [r for r in read_tsv(screen_path) if r["flag_1_5"] == "True"]
        grouped = collections.defaultdict(list)
        for row in observations:
            assert row["status"] == "CANDIDATE" and row["source_binding"] == "VERIFIED"
            grouped[row["marker"], row["taxon"]].append(row)
        assert len(observations) == 18 and len(grouped) == 13 and len({g for g,t in grouped}) == 12
        queries, source_rows = {}, []
        for i, ((marker, taxon), records) in enumerate(sorted(grouped.items(), key=lambda x: (x[0][1] != DT, x[0])), 1):
            row = records[0]
            assert all(all(r[k] == row[k] for k in ("selected_id", "source_fasta_sha256", "protein_length")) for r in records)
            source = Path(row["source_fasta_path"])
            assert bind(source) == row["source_fasta_sha256"]
            saved = fasta(source)
            assert set(saved) == {row["selected_id"]}
            original = saved[row["selected_id"]]
            sequence = original.rstrip("*")
            assert len(sequence) == int(row["protein_length"]) and "*" not in sequence
            query = f"Q{i:03d}"
            queries[query] = sequence
            source_rows.append({"query": query, "marker": marker, "taxon": taxon, "selected_id": row["selected_id"],
                                "source_type": row["source_type"], "source_fasta_path": str(source),
                                "source_fasta_sha256": row["source_fasta_sha256"], "protein_length": len(sequence),
                                "removed_terminal_stops": len(original)-len(sequence),
                                "query_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                "panels": ";".join(r["panel"] for r in records)})
        write_tsv(output / "QUERY_SOURCES.tsv", source_rows)
        query_path = output / "queries.faa"; write_fasta(query_path, queries)
        if args.candidate_evidence:
            bind(args.candidate_evidence)
            evidence = json.loads(args.candidate_evidence.read_text())
            assert {(r["marker"],r["taxon"]) for r in evidence["candidate_cells"]} == set(grouped)
            for path, info in evidence["files"].items():
                assert hashlib.sha256(info["text"].encode()).hexdigest() == info["sha256"]
            by_pair = {(r["marker"], r["taxon"]):r for r in source_rows}
            bound_domains = []; seen = set()
            for observation in evidence["candidate_cells"]:
                source = by_pair[observation["marker"], observation["taxon"]]
                assert all(str(source[k]) == observation[k] for k in ("selected_id", "source_fasta_sha256", "protein_length"))
                for file in json.loads(observation["hmmer_files_json"]):
                    saved = evidence["files"][file["path"]]
                    assert saved["sha256"] == file["sha256"]
                    selected_rows = [line.split() for line in saved["text"].splitlines() if line.strip()
                                     and not line.startswith("#") and line.split()[0] == source["selected_id"]]
                    assert len(selected_rows) == file["selected_target_rows"]
                    for fields in selected_rows:
                        key = (source["query"], tuple(fields[:22]))
                        if key in seen:
                            continue
                        seen.add(key)
                        bound_domains.append({"query":source["query"],"marker":source["marker"],"taxon":source["taxon"],
                            "selected_id":source["selected_id"],"domain_number":int(fields[9]),"hmm_model_length":int(fields[5]),
                            "independent_evalue":float(fields[12]),"domain_bitscore":float(fields[13]),
                            "hmm_start":int(fields[15]),"hmm_end":int(fields[16]),"protein_start":int(fields[17]),"protein_end":int(fields[18])})
            write_tsv(output / "BOUND_HMM_DOMAINS.tsv", bound_domains)
            (output / "ARCHIVED_HMM_EVIDENCE_BINDING.json").write_text(json.dumps({"sha256": sha(args.candidate_evidence),
                "file_hashes": {p: i["sha256"] for p,i in evidence["files"].items()}, "unique_queries": len(queries)}, indent=2)+"\n")
        proteome = args.wgs_root / "braker.aa"
        gtf = args.wgs_root / "assets_20260914/annotation_20250225/braker.gtf"
        genome = args.wgs_root / "genome.fa"; fai = Path(str(genome)+".fai")
        assert bind(proteome) == EXPECTED_PROTEOME
        assert bind(gtf) == EXPECTED_GTF
        assert fai.exists(), "Pre-existing reference index required; no index creation on inputs"
        bind(genome); bind(fai)
        subjects = fasta(proteome)
        versions = {"blastp": run_command([args.blastp,"-version"], receipt),
                    "samtools": run_command([args.samtools,"--version"], receipt)}
        (output / "VERSIONS.json").write_text(json.dumps(versions, indent=2)+"\n")
        (output / "blastp_help.txt").write_text(run_command([args.blastp,"-help"], receipt))
        for label, subject, n in (("versus_current", proteome, len(subjects)), ("candidate_self", query_path, len(queries))):
            command = [args.blastp, "-query", str(query_path), "-subject", str(subject), "-evalue", "1e-10", "-seg", "yes",
                       "-max_target_seqs", str(n), "-outfmt", "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
                       "-out", str(output / (label+".tsv"))]
            run_command(command, receipt)
        hits = parse_blast(output / "versus_current.tsv", queries, subjects)
        self_hits = parse_blast(output / "candidate_self.tsv", queries, queries)
        (output / "HSP_VALIDATION.json").write_text(json.dumps({"status":"PASS","current_HSPs":len(hits),"self_HSPs":len(self_hits),
            "checks":["exact query and subject FASTA IDs and lengths", "both aligned strings bound to sequence coordinates", "alignment lengths",
                      "identity recomputed from aligned residues", "mismatch count", "gap-open count"],
            "limits":"Sequence/coordinate consistency of local alignments does not certify orthology or full-length structural correspondence"},indent=2)+"\n")
        target_rows, query_rows, repeats = coverage_tables(hits, queries)
        labels = {r["query"]: r for r in source_rows}
        for table_rows in (target_rows, query_rows, repeats):
            for row in table_rows:
                row.update({k: labels[row["query"]][k] for k in ("marker", "taxon", "selected_id")})
        write_tsv(output / "CURRENT_TARGET_COVERAGE.tsv", target_rows)
        write_tsv(output / "QUERY_COVERAGE.tsv", query_rows)
        write_tsv(output / "REUSED_CURRENT_REGIONS.tsv", repeats)
        within_self = [r for r in self_hits if r["query"] == r["target"] and r["aligned"] >= 100
                       and r["identity_recomputed"] >= 80 and (r["qstart"],r["qend"]) < (r["sstart"],r["send"])]
        write_tsv(output / "WITHIN_QUERY_OFFDIAGONAL_HSPS.tsv", [{k:v for k,v in r.items() if k not in ("qseq","sseq")} for r in within_self])
        targets = {h["target"] for h in hits}
        current_models, selected_lines = annotation_records(gtf, targets)
        (output / "matched_current_annotation.gtf").write_text("\n".join(selected_lines)+"\n")
        write_fasta(output / "matched_current_proteins.faa", {t:subjects[t] for t in sorted(targets)})
        annotation_paths = set(args.historical_annotation)
        for pattern in ("*.gtf", "*.gff", "*.gff3", "*.gtf.gz", "*.gff.gz", "*.gff3.gz"):
            annotation_paths.update(args.wgs_root.glob(pattern))
            annotation_paths.update((args.wgs_root / "assets_20260914").rglob(pattern))
        wanted = {r["selected_id"] for r in source_rows if r["taxon"] == DT}
        old_models, inventory, old_lines = [], [], []
        for path in sorted(annotation_paths):
            entry = {"path": str(path), "status": "SEARCHED", "exact_old_identifiers_with_CDS": 0}
            try:
                entry["sha256"] = bind(path)
                models, lines = annotation_records(path, wanted)
                entry["exact_old_identifiers_with_CDS"] = len(models)
                for ident, features in models.items():
                    old_models.append({"path": str(path), "old_id": ident, "features": features})
                old_lines.extend(["# source_annotation="+str(path)] + lines if lines else [])
            except (OSError, ValueError, UnicodeError) as error:
                entry.update(status="NOT_VERIFIED", reason=type(error).__name__+":"+str(error))
            inventory.append(entry)
        write_tsv(output / "HISTORICAL_ANNOTATION_SEARCH.tsv", inventory)
        (output / "historical_exact_id_records.gff").write_text("\n".join(old_lines)+"\n")
        reference_contigs = {line.split("\t")[0] for line in fai.read_text().splitlines()}
        regions = set()
        all_features = [f for features in current_models.values() for f in features] + [f for model in old_models for f in model["features"]]
        for f in all_features:
            if f["contig"] in reference_contigs:
                regions.add(f'{f["contig"]}:{f["start"]}-{f["end"]}')
        region_path = output / "cds_regions.txt"; region_path.write_text("\n".join(sorted(regions))+"\n")
        exons = {}
        if regions:
            run_command([args.samtools,"faidx","-r",str(region_path),str(genome)],receipt,output / "cds_reference_sequences.fna")
            exons = fasta(output / "cds_reference_sequences.fna")
            assert set(exons) == regions
        translations = []
        for target in sorted(targets):
            features = current_models.get(target)
            validation = translate_model(features, exons, subjects[target]) if features else {"status":"NOT_VERIFIED","reason":"no_exact_current_CDS_model"}
            translations.append({"current_protein":target, **validation})
        write_tsv(output / "CURRENT_MODEL_TRANSLATION.tsv", translations)
        historical = []
        for source in source_rows:
            if source["taxon"] != DT:
                historical.append({"query":source["query"],"selected_id":source["selected_id"],"status":"NOT_APPLICABLE_OTHER_SPECIES",
                                   "reason":"Current DT homology comparison does not recover another species original coordinates"})
                continue
            found = [m for m in old_models if m["old_id"] == source["selected_id"]]
            if not found:
                historical.append({"query":source["query"],"selected_id":source["selected_id"],"status":"NOT_RECOVERED_IN_SEARCHED_FILES",
                                   "searched_annotation_files":len(inventory),"historical_file_provenance":"NOT_VERIFIED",
                                   "reason":"No exact historical protein identifier with CDS records; no alias inferred from current gene numbering"})
            for model in found:
                historical.append({"query":source["query"],"selected_id":source["selected_id"],"annotation_path":model["path"],
                                   "historical_file_provenance":"NOT_VERIFIED_UNLESS_INDEPENDENTLY_ESTABLISHED",
                                   **translate_model(model["features"],exons,queries[source["query"]])})
        write_tsv(output / "HISTORICAL_COORDINATE_RECOVERY.tsv", historical)
        receipt["final_input_sha256"] = {str(path): sha(path) for path in (proteome, gtf, genome, fai)}
        assert all(digest == receipt["inputs_sha256"][path] for path, digest in receipt["final_input_sha256"].items()), "Core input changed during diagnostic"
        receipt.update(state="COMPLETE_DIAGNOSTIC",query_count=len(queries),unique_markers=len({r["marker"] for r in source_rows}),
                       historical_DT_queries=sum(r["taxon"]==DT for r in source_rows),current_proteins=len(subjects),
                       current_HSPs=len(hits),self_HSPs=len(self_hits),all_HSPs_sequence_identity_validated=True,
                       current_hit_proteins=len(targets),current_translation_status_counts=dict(collections.Counter(r["status"] for r in translations)),
                       historical_recovery_status_counts=dict(collections.Counter(r["status"] for r in historical)),
                       historical_search_files=len(inventory),historical_search_failures=sum(x["status"]!="SEARCHED" for x in inventory))
    except Exception as error:
        receipt.update(state="FAILED",error_type=type(error).__name__,error=str(error))
        raise
    finally:
        receipt.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_seconds=time.monotonic()-start,
                       limitations=["Local BLAST similarity and pooled coverage do not prove orthology or fusion", "Coverage tiers are descriptive, not calibrated classification",
                                    "Strict nonoverlapping HSP chains can underestimate coherent coverage", "Exact current translation does not establish expression or old model provenance",
                                    "Historical-coordinate recovery is limited to explicitly inventoried files", "No source sequences, annotations, trees or matrices were corrected"])
        receipt["output_sha256"] = {p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file() and p.name!="RECEIPT.json"}
        (output / "RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
        print(json.dumps(receipt,indent=2),flush=True)


if __name__ == "__main__":
    main()
