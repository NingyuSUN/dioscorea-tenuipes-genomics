#!/usr/bin/env python3
"""Full paired FASTQ QC, Python stdlib only; inputs are opened read-only.

Four physical lines per record are required (wrapped FASTQ is unsupported).
LF/CRLF are accepted; a complete final quality line need not end in LF.
Identifier conventions: /1,/2; CASAVA mate fields; ENA/SRA accession.spot
with optional .1/.2 read suffix and optional separate 1/2 mate field.
Counts are paired records, not independently verified unique molecules.
"""
import argparse
from collections import Counter
import gzip
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zlib

CHUNK = 1024 * 1024
VALID_QUALITY = bytes(range(33, 127))
SRA_READ_SUFFIX = re.compile(rb"^([SED]RR[0-9]+\.[0-9]+)\.([12])$")


class AuditError(Exception):
    pass


def snapshot(path):
    s = path.stat()
    return {"path": str(path), "compressed_bytes": s.st_size,
            "mtime_ns": s.st_mtime_ns, "device": s.st_dev, "inode": s.st_ino}


class InputStream:
    def __init__(self, path, backend):
        self.path, self.backend = path, backend
        self.process = self.stderr_file = self.stream = None
        self.completed = False

    def __enter__(self):
        if self.backend == "pigz":
            self.stderr_file = tempfile.TemporaryFile()
            try:
                self.process = subprocess.Popen(
                    [shutil.which("pigz"), "-dc", "--", str(self.path)],
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                    stderr=self.stderr_file, bufsize=CHUNK)
                self.stream = self.process.stdout
            except Exception:
                self.stderr_file.close()
                raise
        else:
            self.stream = io.BufferedReader(gzip.open(self.path, "rb"), CHUNK)
        return self

    def finish(self):
        # Called only after both FASTQ streams have reached EOF. gzip checks CRC
        # and ISIZE while reading; pigz failures must also be checked explicitly.
        if self.process is not None:
            code = self.process.wait()
            self.stderr_file.seek(0)
            detail = self.stderr_file.read(8192).decode("utf-8", "replace").strip()
            if code != 0:
                raise AuditError("pigz decompression failed for %s (exit %s): %s"
                                 % (self.path, code, detail))
        self.completed = True

    def __exit__(self, exc_type, exc, tb):
        if self.stream is not None:
            self.stream.close()
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        if self.stderr_file is not None:
            self.stderr_file.close()


def without_eol(line):
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    return line


def identifier(payload, expected_mate, context):
    fields = payload.split()
    if not fields:
        raise AuditError(context + ": empty read identifier")
    token = fields[0]
    declared = []
    if token.endswith((b"/1", b"/2")):
        declared.append(int(token[-1:]))
        token = token[:-2]
    else:
        match = SRA_READ_SUFFIX.fullmatch(token)
        if match:
            token = match.group(1)
            declared.append(int(match.group(2)))
    if len(fields) > 1:
        field = fields[1]
        if field in (b"1", b"2"):
            declared.append(int(field))
        elif field.startswith((b"1:", b"2:")):
            declared.append(int(field[:1]))
    if not token or any(mate != expected_mate for mate in declared):
        raise AuditError(context + ": mate annotation conflicts with file R%d" % expected_mate)
    return token


def record(stream, mate, index, valid_bases):
    context = "R%d pair %d" % (mate, index)
    header = stream.readline()
    if not header:
        return None
    if not header.startswith(b"@"):
        raise AuditError(context + ": expected @ header; wrapped FASTQ unsupported")
    key = identifier(without_eol(header)[1:], mate, context)
    seq_line, plus_line, qual_line = stream.readline(), stream.readline(), stream.readline()
    if not seq_line or not plus_line or not qual_line:
        raise AuditError(context + ": truncated four-line FASTQ record")
    seq, plus, quality = map(without_eol, (seq_line, plus_line, qual_line))
    if not seq:
        raise AuditError(context + ": empty sequence")
    if not plus.startswith(b"+"):
        raise AuditError(context + ": expected + separator; wrapped FASTQ unsupported")
    if plus[1:].strip():
        repeat = plus[1:]
        if repeat.startswith(b"@"):
            repeat = repeat[1:]
        if identifier(repeat, mate, context + " + separator") != key:
            raise AuditError(context + ": repeated + identifier differs from header")
    if len(seq) != len(quality):
        raise AuditError(context + ": sequence/quality length mismatch (%d/%d)"
                         % (len(seq), len(quality)))
    invalid_bases = seq.translate(None, valid_bases)
    if invalid_bases:
        raise AuditError(context + ": invalid sequence byte(s) " + repr(invalid_bases[:24]))
    invalid_quality = quality.translate(None, VALID_QUALITY)
    if invalid_quality:
        raise AuditError(context + ": invalid Phred+33 byte(s) " + repr(invalid_quality[:24]))
    return key, seq, quality


class Statistics:
    def __init__(self):
        self.reads = self.bases = 0
        self.lengths, self.base_hist, self.quality_hist = Counter(), Counter(), Counter()
        self.base_buffer, self.quality_buffer = bytearray(), bytearray()

    def add(self, seq, quality):
        self.reads += 1
        self.bases += len(seq)
        self.lengths[len(seq)] += 1
        self.base_buffer.extend(seq.upper())
        self.quality_buffer.extend(quality)
        if len(self.base_buffer) >= CHUNK:
            self.flush()

    def flush(self):
        # Large buffered updates avoid two Counter.update calls per read.
        self.base_hist.update(self.base_buffer)
        self.quality_hist.update(self.quality_buffer)
        self.base_buffer.clear()
        self.quality_buffer.clear()

    def summary(self):
        self.flush()
        return summarize(self.reads, self.bases, self.lengths, self.base_hist, self.quality_hist)


def fraction(numerator, denominator):
    return numerator / denominator if denominator else None


def summarize(reads, bases, lengths, base_hist, quality_hist):
    if sum(base_hist.values()) != bases or sum(quality_hist.values()) != bases:
        raise AuditError("Internal base/quality histogram accounting mismatch")
    q20 = sum(n for byte, n in quality_hist.items() if byte >= 53)
    q30 = sum(n for byte, n in quality_hist.items() if byte >= 63)
    gc = base_hist[ord("G")] + base_hist[ord("C")]
    acgt = sum(base_hist[ord(base)] for base in "ACGT")
    return {"mate_records": reads, "bases": bases,
            "read_length_histogram": {str(k): v for k, v in sorted(lengths.items())},
            "base_histogram": {chr(k): v for k, v in sorted(base_hist.items())},
            "phred33_quality_histogram": {str(k - 33): v for k, v in sorted(quality_hist.items())},
            "GC_bases": gc, "N_bases": base_hist[ord("N")], "ACGT_bases": acgt,
            "non_ACGT_bases": bases - acgt,
            "GC_fraction_all_bases": fraction(gc, bases),
            "GC_fraction_ACGT_bases": fraction(gc, acgt),
            "N_fraction_all_bases": fraction(base_hist[ord("N")], bases),
            "Q20_bases": q20, "Q30_bases": q30,
            "Q20_fraction_all_bases": fraction(q20, bases),
            "Q30_fraction_all_bases": fraction(q30, bases),
            "mean_Phred33": fraction(sum((k - 33) * n for k, n in quality_hist.items()), bases)}


def audit(args):
    started = time.monotonic()
    paths = [Path(args.r1).resolve(), Path(args.r2).resolve()]
    stats = [Statistics(), Statistics()]
    pairs, complete, error, inputs, unchanged = 0, False, None, [], None
    backend = "pigz" if args.decompressor == "auto" and shutil.which("pigz") else args.decompressor
    if backend == "auto":
        backend = "python"
    alphabet = b"ACGTN" if args.alphabet == "acgtn" else b"ACGTRYSWKMBDHVN"
    valid_bases = alphabet + alphabet.lower()
    try:
        if backend == "pigz" and not shutil.which("pigz"):
            raise AuditError("pigz requested but not found on PATH")
        inputs = [snapshot(p) for p in paths]
        if os.path.samefile(*paths):
            raise AuditError("R1 and R2 must be different input files")
        with InputStream(paths[0], backend) as left, InputStream(paths[1], backend) as right:
            while True:
                first = record(left.stream, 1, pairs + 1, valid_bases)
                second = record(right.stream, 2, pairs + 1, valid_bases)
                if first is None or second is None:
                    if first is not None or second is not None:
                        raise AuditError("Unequal mate record counts at pair %d" % (pairs + 1))
                    break
                if first[0] != second[0]:
                    raise AuditError("Paired identifier mismatch at pair %d: %r vs %r"
                                     % (pairs + 1, first[0][:160], second[0][:160]))
                stats[0].add(first[1], first[2])
                stats[1].add(second[1], second[2])
                pairs += 1
                if args.progress_every_pairs and pairs % args.progress_every_pairs == 0:
                    elapsed = time.monotonic() - started
                    bases = stats[0].bases + stats[1].bases
                    remaining = None
                    if args.expected_bases is not None and bases and bases <= args.expected_bases:
                        remaining = elapsed * (args.expected_bases - bases) / bases
                    print(json.dumps({"progress": "validated_complete_pairs", "pairs": pairs,
                                      "bases": bases, "elapsed_seconds": round(elapsed, 2),
                                      "estimated_seconds_remaining_from_base_rate": remaining}),
                          file=sys.stderr, flush=True)
            left.finish()
            right.finish()
        complete = True
        unchanged = inputs == [snapshot(p) for p in paths]
        if not unchanged:
            raise AuditError("An input file changed during the audit")
        if pairs == 0:
            raise AuditError("Empty paired FASTQ input")
    except (AuditError, OSError, EOFError, zlib.error) as exc:
        error = "%s: %s" % (type(exc).__name__, exc)
    per_mate = [s.summary() for s in stats]
    combined = summarize(sum(s.reads for s in stats), sum(s.bases for s in stats),
                         stats[0].lengths + stats[1].lengths,
                         stats[0].base_hist + stats[1].base_hist,
                         stats[0].quality_hist + stats[1].quality_hist)
    expected = {"official_expected_bases": args.expected_bases,
                "ENA_reported_read_count": args.ena_read_count,
                "bases_match": None if args.expected_bases is None else combined["bases"] == args.expected_bases,
                "ENA_count_matches_pair_records": None if args.ena_read_count is None else pairs == args.ena_read_count,
                "ENA_count_matches_mate_records": None if args.ena_read_count is None else combined["mate_records"] == args.ena_read_count}
    if error is None and expected["bases_match"] is False:
        error = "Expected total base count does not match the complete files"
    if (error is None and args.ena_read_count is not None and
            not (expected["ENA_count_matches_pair_records"] or expected["ENA_count_matches_mate_records"])):
        error = "ENA read_count matches neither paired records nor mate records"
    passed = error is None
    result = {"schema": "paired_fastq_technical_qc_v1", "label": args.label,
              "status": "PASS_TECHNICAL_QC" if passed else "FAIL_TECHNICAL_QC",
              "scan_reached_both_EOF_and_decompression_passed": complete,
              "decompressor": backend, "python_version": sys.version.split()[0],
              "elapsed_seconds": time.monotonic() - started, "inputs": inputs,
              "input_stat_unchanged": unchanged, "error": error,
              "paired_records": pairs, "mate_records": combined["mate_records"],
              "counts_cover_complete_validated_pairs_only": True,
              "per_mate": {"R1": per_mate[0], "R2": per_mate[1]}, "combined": combined,
              "external_expectations": expected,
              "checks": {"full_paired_identifier_order": complete and error is None,
                         "four_physical_lines_per_record": complete,
                         "sequence_quality_lengths_equal": complete,
                         "base_alphabet_valid": complete,
                         "quality_ASCII_33_through_126_valid": complete,
                         "gzip_end_and_CRC_valid": complete,
                         "global_identifier_uniqueness_checked": False},
              "format_policy": {"wrapped_FASTQ_supported": False, "CRLF_supported": True,
                                "final_quality_line_without_newline_supported": True,
                                "base_alphabet_case_insensitive": alphabet.decode("ascii"),
                                "quality_encoding_assumed": "Phred+33; legal score range 0..93; encoding not independently inferred"},
              "limits": ["Technical integrity and quality only; not contamination, taxonomic identity, biological independence, mapping coverage or marker readiness.",
                         "Read pairs are paired records, not verified unique molecules; duplicate identifiers and sequence duplicates are not globally assessed.",
                         "Q20/Q30 fractions use all validated sequence bases including N; GC is reported with both all-base and ACGT-only denominators.",
                         "On failure counts cover complete validated pairs before failure; they are not complete-file totals."]}
    return result, 0 if passed else 2


def positive_int(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--output", help="JSON destination; default is stdout")
    parser.add_argument("--label", default="paired_FASTQ")
    parser.add_argument("--expected-bases", type=positive_int)
    parser.add_argument("--ena-read-count", type=positive_int)
    parser.add_argument("--decompressor", choices=("auto", "python", "pigz"), default="auto")
    parser.add_argument("--alphabet", choices=("acgtn", "iupac"), default="iupac")
    parser.add_argument("--progress-every-pairs", type=positive_int, default=1000000)
    args = parser.parse_args()
    # Refuse output aliases of either input before opening the output path.
    if args.output:
        output = Path(args.output).resolve()
        for raw in (args.r1, args.r2):
            candidate = Path(raw).resolve()
            if output == candidate or (output.exists() and candidate.exists() and os.path.samefile(output, candidate)):
                parser.error("output must not overwrite either input")
    result, code = audit(args)
    serialized = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if args.output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=output.parent,
                                         prefix=output.name + ".", suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(serialized)
        os.replace(temporary, output)
    else:
        sys.stdout.write(serialized)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
