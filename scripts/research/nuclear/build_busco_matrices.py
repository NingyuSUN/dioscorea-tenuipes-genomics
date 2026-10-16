"""Build audited genome-only and genome+transcriptome nuclear amino-acid matrices.

Orthology is defined by frozen BUSCO lineage profiles. Duplicated and fragmented
hits are missing, never replaced with a best-scoring copy. Decisions precede trees.
"""
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import argparse, hashlib, json, math, subprocess, time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

EXCLUDED_MARKERS = {'2661at3193'}
REQUIRED_FOCAL = {'Dioscorea_tenuipes_historical', 'Dioscorea_gracillima'}
AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
R = Path(__file__).resolve().parents[1]
MAFFT = "mafft"


def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read_fasta(p):
    with p.open() as f:
        records = list(SeqIO.parse(f, 'fasta'))
    assert len({x.id for x in records}) == len(records), p
    return {x.id: str(x.seq).upper() for x in records}


def write_fasta(records, p):
    SeqIO.write([SeqRecord(Seq(s), id=t, description='')
                 for t, s in records.items()], p, 'fasta')


def check_table_binding(fields, ident, genome_mode):
    if not genome_mode:
        assert fields[2] == ident, ('Protein/table mismatch', fields[0], ident, fields[2])
        return
    target,region,strand=ident.rsplit('|',2)
    contig,coords=region.rsplit(':',1)
    start,end=map(int,coords.split('-'))
    assert target.split('_')[0]==fields[0]
    assert contig==fields[2] and strand==fields[5]
    assert sorted([start,end])==sorted([int(fields[3]),int(fields[4])])


def read_busco(run):
    summaries = list(run.glob('short_summary*.json'))
    assert len(summaries) == 1, (run, summaries)
    summary = json.loads(summaries[0].read_text())
    results = run / 'run_embryophyta_odb10'
    table = results / 'full_table.tsv'
    states = {}
    table_rows = {}
    genome_mode = "Gene Start\tGene End" in table.read_text()
    for line in table.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        marker, status = fields[:2]
        if marker in states:
            assert status == states[marker] == 'Duplicated', (run, marker)
        states[marker] = status
        if status == "Complete": table_rows[marker] = fields
    assert len(states) == 1614, (run, len(states))
    complete = {g for g, s in states.items() if s == 'Complete'}
    seqdir = results / 'busco_sequences/single_copy_busco_sequences'
    files = {p.stem: p for p in seqdir.glob('*.faa')}
    assert set(files) == complete, (run, set(files) ^ complete)
    excluded_stops = set()
    excluded_unbound = set()
    stop_audit = None
    if genome_mode:
        from audit_miniprot_stops import audit as audit_stops
        stop_audit = audit_stops(run)
        excluded_stops = set(stop_audit['excluded_markers'])
        excluded_unbound = set(stop_audit['excluded_unbound_markers'])
    sequences, audit = {}, []
    for marker, p in sorted(files.items()):
        records = read_fasta(p)
        assert len(records) == 1, p
        ident, s = next(iter(records.items()))
        check_table_binding(table_rows[marker], ident, genome_mode)
        s = s.rstrip('*')
        reason = None
        if marker in EXCLUDED_MARKERS:
            reason = 'user_excluded_annotation_uncertain_orthogroup'
        elif marker in excluded_unbound:
            reason = 'postprocessed_gene_not_exactly_bound_to_raw_translation'
        elif marker in excluded_stops:
            reason = 'internal_stop_in_raw_miniprot_translation'
        elif '*' in s:
            reason = 'internal_stop'
        elif len(s) < 50:
            reason = 'less_than_50_residues'
        elif sum(c in AA for c in s) / len(s) < .95:
            reason = 'less_than_95_percent_canonical_amino_acids'
        audit.append({'marker': marker, 'id': ident, 'length': len(s),
                      'source': str(p), 'source_sha256': sha(p), 'excluded_reason': reason})
        if reason is None:
            sequences[marker] = ''.join(c if c in AA else 'X' for c in s)
    return sequences, {'summary_file': str(summaries[0]), 'summary_sha256': sha(summaries[0]),
                       'full_table_sha256': sha(table), 'status_counts': dict(Counter(states.values())),
                       'summary_results': summary.get('results'), 'sequence_audit': audit,
                       'usable_single_copy': len(sequences), 'raw_translation_stop_audit': stop_audit}


def trim_gene(alignment, taxa, occupancy=.8):
    """Return trimmed records or a documented exclusion; absent copies stay missing."""
    assert set(alignment) <= set(taxa)
    lengths = {len(s) for s in alignment.values()}
    assert len(lengths) == 1
    n = lengths.pop()
    minimum = math.ceil(occupancy * len(taxa))
    if len(alignment) < minimum:
        return None, {'reason': 'taxon_occupancy', 'present_taxa': len(alignment)}
    cols = [i for i in range(n) if sum(s[i] in AA for s in alignment.values()) >= minimum]
    if len(cols) < 50:
        return None, {'reason': 'fewer_than_50_retained_columns', 'retained_columns': len(cols)}
    records = {t: ''.join(alignment[t][i] if alignment[t][i] in AA else '-' for i in cols)
               if t in alignment else '-' * len(cols) for t in taxa}
    coverage = {t: sum(c in AA for c in records[t]) / len(cols) for t in taxa}
    bad = [t for t in alignment if coverage[t] < .5]
    if bad:
        return None, {'reason': 'present_sequence_under_50_percent_coverage', 'taxa': bad}
    return records, {'reason': None, 'input_columns': n, 'retained_columns': len(cols),
                     'retained_column_indices_0based': cols, 'coverage': coverage,
                     'present_taxa': len(alignment), 'complete_taxon_occupancy': len(alignment) == len(taxa)}


def build_panel(name, source_taxa, dt, threads):
    dest = R / 'matrices' / name
    dest.mkdir(parents=True, exist_ok=True)
    raw, aligned, trimmed = (dest / x for x in ['raw', 'aligned', 'trimmed'])
    for p in [raw, aligned, trimmed]:
        p.mkdir(exist_ok=True)
    seqs, audits, excluded = {}, {}, []
    for row in source_taxa:
        taxon = row['taxon']
        s, a = read_busco(Path(row['busco_run']))
        audits[taxon] = a
        if len(s) < 200:
            excluded.append({'taxon': taxon, 'reason': 'fewer_than_200_usable_single_copy_markers', 'n': len(s)})
        else:
            seqs[taxon] = s
    assert REQUIRED_FOCAL <= set(seqs), ('Required focal taxon failed eligibility', REQUIRED_FOCAL-set(seqs))
    for required in ['Oryza_sativa', 'Arabidopsis_thaliana']:
        assert required in seqs, ('Outgroup failed input gate', required)
    taxa = list(seqs)
    markers = sorted(set().union(*(s.keys() for s in seqs.values())))
    markers = [g for g in markers if all(g in seqs[t] for t in REQUIRED_FOCAL) and g not in EXCLUDED_MARKERS and
               sum(g in seqs[t] for t in taxa) >= math.ceil(.8 * len(taxa))]
    commands = []
    def align(g):
        records = {t: seqs[t][g] for t in taxa if g in seqs[t]}
        input_file, output_file = raw / (g + '.faa'), aligned / (g + '.faa')
        write_fasta(records, input_file)
        cmd = [MAFFT, '--localpair', '--maxiterate', '1000',
               '--thread', '1', '--inputorder', str(input_file)]
        # Files alone are not a cache key: bind alignment to exact input and command.
        stamp = output_file.with_suffix('.receipt.json')
        binding = {'input_sha256': sha(input_file), 'command': cmd}
        if not (stamp.exists() and json.loads(stamp.read_text()) == binding and output_file.exists()):
            with output_file.open('w') as o, output_file.with_suffix('.log').open('w') as e:
                subprocess.run(cmd, stdout=o, stderr=e, check=True)
            stamp.write_text(json.dumps(binding, indent=2))
        a = read_fasta(output_file)
        assert set(a) == set(records)
        for t, s in records.items():
            assert a[t].replace('-', '') == s, (g, t, 'alignment changed residues')
        commands.append(binding)
        return g, a
    print('ALIGN_START', name, len(taxa), len(markers), time.ctime(), flush=True)
    with ThreadPoolExecutor(max_workers=threads) as ex:
        alignments = dict(ex.map(align, markers))
    rounds = []
    while True:
        retained, details = {}, {}
        for g, a in alignments.items():
            a = {t: s for t, s in a.items() if t in taxa}
            # Require DT and at least one outgroup for every retained locus.
            if dt not in a or not ({'Oryza_sativa', 'Arabidopsis_thaliana'} & set(a)):
                continue
            recs, audit = trim_gene(a, taxa)
            details[g] = audit
            if recs is not None:
                retained[g] = recs
        assert len(retained) >= 100, (name, 'fewer than 100 loci', len(retained))
        matrix = {t: ''.join(retained[g][t] for g in sorted(retained)) for t in taxa}
        coverage = {t: sum(c in AA for c in s) / len(s) for t, s in matrix.items()}
        low = [t for t, c in coverage.items() if c < .5]
        rounds.append({'taxa': list(taxa), 'genes': len(retained), 'coverage': coverage, 'removed': low})
        assert not (REQUIRED_FOCAL & set(low)) and 'Oryza_sativa' not in low and 'Arabidopsis_thaliana' not in low
        if not low:
            break
        excluded.extend({'taxon': t, 'reason': 'concatenation_under_50_percent_residue_coverage',
                         'coverage': coverage[t]} for t in low)
        taxa = [t for t in taxa if t not in low]
    assert len(taxa) >= 7
    write_fasta(matrix, dest / 'main.faa')
    write_fasta({t: s for t, s in matrix.items() if t != 'Arabidopsis_thaliana'}, dest / 'no_arabidopsis.faa')
    strict = [g for g in sorted(retained) if details[g]['complete_taxon_occupancy']]
    write_fasta({t: ''.join(retained[g][t] for g in strict) for t in taxa}, dest / 'strict.faa')
    offset, parts = 1, []
    for g, records in sorted(retained.items()):
        write_fasta(records, trimmed / (g + '.faa'))
        n = len(next(iter(records.values())))
        parts.append({'gene': g, 'start_1based': offset, 'end_1based': offset + n - 1,
                      'length': n, 'sha256': sha(trimmed / (g + '.faa'))})
        offset += n
    seq_list = list(matrix.values())
    variable = informative = 0
    for chars in zip(*seq_list):
        counts = Counter(c for c in chars if c in AA)
        variable += len(counts) > 1
        informative += sum(n >= 2 for n in counts.values()) >= 2
    summary = {'required_focal_taxa_per_locus': sorted(REQUIRED_FOCAL), 'excluded_orthogroups_all_taxa': sorted(EXCLUDED_MARKERS), 'panel': name, 'taxa': taxa, 'source_taxa': source_taxa, 'excluded_taxa': excluded,
               'taxon_eligibility_minimum_single_copy': 200, 'gene_and_column_occupancy': .8,
               'minimum_gene_columns': 50, 'minimum_present_sequence_coverage': .5,
               'minimum_final_taxon_coverage': .5, 'selection_rounds': rounds,
               'genes': len(retained), 'columns': offset - 1, 'variable_sites': variable,
               'parsimony_informative_sites': informative, 'strict_genes': len(strict),
               'strict_columns': sum(parts[i]['length'] for i in range(len(parts)) if parts[i]['gene'] in strict),
               'main_sha256': sha(dest / 'main.faa'), 'partitions': parts,
               'strict_marker_ids': strict, 'gene_audit': details, 'input_audit': audits,
               'orthology_limit': 'BUSCO profile-defined complete single-copy hits; no claim of proven orthology or genome completeness from filtered matrix'}
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2))
    (dest / 'commands.json').write_text(json.dumps(commands, indent=2))
    print('MATRIX_DONE', name, len(taxa), len(retained), offset - 1, 'strict', len(strict), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--work-root', type=Path, required=True)
    parser.add_argument('--mafft', default='mafft')
    parser.add_argument('--panel', required=True)
    parser.add_argument('--taxa-json', type=Path, required=True)
    parser.add_argument('--dt', default='Dioscorea_tenuipes_historical')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    R = args.work_root.resolve()
    MAFFT = args.mafft
    build_panel(args.panel, json.loads(args.taxa_json.read_text()), args.dt, args.threads)
