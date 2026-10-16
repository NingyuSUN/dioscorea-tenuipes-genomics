#!/usr/bin/env python3
"""Conservative, UNPHASED partial CP sequence on verified source coordinates.
Preparation only: run explicitly after pilot QC. No tree, reference filling,
assembly, deduplication, or full-plastome claim. A separately declared,
hash-verified SSC reverse-complement creates the normalized path2 output.
"""
import argparse
from array import array
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

EXPECTED_LENGTH = 153650
EXPECTED_SEQUENCE_SHA256 = '2ea9415d91a62bbf3ba0b86b41e0162679c3c5fe6c3ef8936f41ab8c9d6078f9'
NORMALIZED_PATH2_SEQUENCE_SHA256 = '9b650ec46a52bcf6bfe5857e5c583b13106c98ab9dd3c4fb5b03a64f001890cb'
BASES = 'ACGT'
THRESHOLDS = {
    'min_base_quality': 20, 'min_mapping_quality': 20, 'min_ACGT_template_depth': 10,
    'dominant_fraction_numerator': 9, 'dominant_fraction_denominator': 10,
    'min_dominant_support_each_strand': 2,
    'indel_risk_numerator': 1, 'indel_risk_denominator': 10,
    'indel_padding_bp': 5, 'linear_end_mask_bp': 200,
    'pileup_max_depth': 100000, 'threads': 1, 'address_space_limit_bytes': 2 * 1024**3,
    'excluded_flags': 0xF04,  # UNMAP, SECONDARY, QCFAIL, DUP, SUPPLEMENTARY
    'exclude_unknown_MAPQ_255': True, 'exclude_CIGAR_refskip_reads': True,
    'compute_baq': False, 'ignore_orphans': False, 'native_ignore_overlaps': False,
    'overlap_policy': 'RG+QNAME once per position; disagreeing eligible mates -> conflict; agreeing mates -> higher BQ, read1 on ties; one strand only',
    'dominant_fraction_denominator_definition': 'ACGT observations + ambiguous N + deletion + overlap conflicts, after quality filters and template collapse',
    'DP_definition': 'count of A/C/G/T observations after template collapse; excludes N/deletion/conflict',
    'indel_anchor_denominator_definition': 'distinct templates spanning anchor and next aligned query base with both BQ>=20; terminal reads do not count as no-indel',
    'indel_anchor_numerator_definition': 'any qualifying nonzero indel from a template; discordant mate event still counts as masking risk; inserted-base quality is not needed for a risk mask',
    'deletion_site_risk_definition': 'deletion templates / all eligible template observations; flanking query bases both BQ>=20; masks site+-5',
    'indel_mask_extent': 'insertion anchor+-5; deletion anchor through deleted span, padded 5; all triggered masks applied before FASTA output',
    'duplicate_handling': 'flagged duplicates excluded; no duplicate marking/removal performed; unmarked PCR/optical duplicates can remain',
    'normalized_path2_transform': 'path1[:109153] + reverse_complement(path1[109153:128097]) + path1[128097:]; no rotation; preserve N; retain both IR copies',
    'normalization_interval_0based_halfopen': [109153, 128097],
}


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def json_new(path, obj):
    with Path(path).open('x') as handle:
        json.dump(obj, handle, indent=2)
        handle.write('\n')


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def normalized_path2(seq):
    if len(seq) != EXPECTED_LENGTH:
        raise ValueError('Path normalization requires exact verified reference length')
    return seq[:109153] + reverse_complement(seq[109153:128097]) + seq[128097:]


def normalized_source_position(pos):
    return 109153 + 128097 - 1 - pos if 109153 <= pos < 128097 else pos


def verified_reference(path):
    lines = Path(path).read_text().splitlines()
    if sum(line.startswith('>') for line in lines) != 1:
        raise ValueError('Source FASTA must have exactly one record')
    seq = ''.join(line.strip() for line in lines if not line.startswith('>')).upper()
    digest = hashlib.sha256(seq.encode('ascii')).hexdigest()
    if len(seq) != EXPECTED_LENGTH or digest != EXPECTED_SEQUENCE_SHA256:
        raise ValueError(f'Unverified source coordinate system: length={len(seq)}, sequence_sha256={digest}')
    if set(seq) - set(BASES):
        raise ValueError('Unexpected source reference alphabet')
    normalized_hash = hashlib.sha256(normalized_path2(seq).encode('ascii')).hexdigest()
    if normalized_hash != NORMALIZED_PATH2_SEQUENCE_SHA256:
        raise ValueError('Declared SSC transform does not reproduce fixed normalized path2 hash')
    return seq, {'path': str(Path(path).resolve()), 'length': len(seq), 'sequence_sha256': digest,
                 'file_sha256': sha256_file(path), 'sequence_normalization': 'remove FASTA header/line whitespace; uppercase',
                 'raw_path1_preserved': True, 'normalized_path2_reference_sha256': normalized_hash,
                 'normalization_formula': THRESHOLDS['normalized_path2_transform']}


def reject_reason(read):
    if read.flag & THRESHOLDS['excluded_flags']:
        return 'excluded_SAM_flag'
    if read.mapping_quality < 20 or read.mapping_quality == 255:
        return 'MAPQ_below20_or_unknown'
    if any(op == 3 for op, length in (read.cigartuples or [])):
        return 'CIGAR_refskip'
    return None


def summarize_column(pileup_reads):
    """Copy scalars now; no PileupRead/Column retained after iterator advances."""
    groups = defaultdict(lambda: {'observations': [], 'events': []})
    raw = 0
    filtered = Counter()
    for pr in pileup_reads:
        raw += 1
        read = pr.alignment
        reason = reject_reason(read)
        if reason:
            filtered[reason] += 1
            continue
        qualities = read.query_qualities
        sequence = read.query_sequence
        if qualities is None or sequence is None:
            filtered['missing_sequence_or_quality'] += 1
            continue
        key = (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name)
        group = groups[key]
        if pr.is_refskip:
            continue
        if pr.is_del:
            qnext = pr.query_position_or_next
            if qnext is None or qnext <= 0 or qnext >= len(qualities):
                continue
            quality = min(qualities[qnext - 1], qualities[qnext])
            if quality >= 20:
                group['observations'].append(('DEL', quality, read.is_reverse, read.is_read1))
            continue
        q = pr.query_position
        if q is None or qualities[q] < 20:
            filtered['base_quality_below20'] += 1
            continue
        base = sequence[q].upper()
        if base not in BASES:
            base = 'N'
        group['observations'].append((base, qualities[q], read.is_reverse, read.is_read1))
        # An eligible risk denominator must span past the anchor, not merely end here.
        indel = pr.indel
        next_query = q + (indel if indel > 0 else 0) + 1
        if next_query < read.query_alignment_end and qualities[next_query] >= 20:
            group['events'].append(indel)
    counts = {base: [0, 0] for base in BASES}  # forward, reverse
    ambiguous = deletion = conflict = overlap = eval_n = event_n = max_deleted = 0
    for group in groups.values():
        obs = group['observations']
        if len(obs) > 1:
            overlap += len(obs) - 1
        if obs:
            observed_states = {item[0] for item in obs}
            if len(observed_states) != 1:
                conflict += 1
            else:
                state, quality, reverse, read1 = max(obs, key=lambda x: (x[1], x[3]))
                if state in BASES:
                    counts[state][int(reverse)] += 1
                elif state == 'DEL':
                    deletion += 1
                else:
                    ambiguous += 1
        if group['events']:
            eval_n += 1
            nonzero = [x for x in group['events'] if x]
            if nonzero:
                event_n += 1
                max_deleted = max(max_deleted, max([-x for x in nonzero if x < 0] or [0]))
    return {'raw_alignments': raw, 'base_strands': counts, 'ambiguous': ambiguous,
            'deletion': deletion, 'conflict': conflict, 'overlap_observations_removed': overlap,
            'indel_evaluable_templates': eval_n, 'indel_event_templates': event_n,
            'indel_max_deleted_span': max_deleted, 'filtered_alignment_observations': dict(filtered)}


def new_evidence(length):
    fields = ['raw_alignments', 'ambiguous', 'deletion', 'conflict', 'overlap_observations_removed',
              'indel_evaluable_templates', 'indel_event_templates', 'indel_max_deleted_span']
    fields += [base + strand for base in BASES for strand in ('_fwd', '_rev')]
    return {field: array('I', [0]) * length for field in fields}


def collect_evidence(bam, contig, length):
    evidence = new_evidence(length)
    filtered = Counter()
    for col in bam.pileup(contig, 0, length, truncate=True, stepper='nofilter',
                          min_base_quality=0, min_mapping_quality=0, flag_filter=0,
                          max_depth=THRESHOLDS['pileup_max_depth'], ignore_overlaps=False,
                          ignore_orphans=False, compute_baq=False, redo_baq=False,
                          adjust_capq_threshold=0):
        if col.nsegments >= THRESHOLDS['pileup_max_depth']:
            raise RuntimeError('Pileup depth reached cap; no silently truncated consensus allowed')
        pos = col.reference_pos
        row = summarize_column(col.pileups)
        for field, values in evidence.items():
            if field.endswith('_fwd') or field.endswith('_rev'):
                values[pos] = row['base_strands'][field[0]][int(field.endswith('_rev'))]
            else:
                values[pos] = row[field]
        filtered.update(row['filtered_alignment_observations'])
    return evidence, dict(filtered)


def column_depth(evidence, pos):
    return sum(evidence[base + strand][pos] for base in BASES for strand in ('_fwd', '_rev'))


def make_masks(evidence, length):
    mask = bytearray(length)
    triggers = []
    for pos in range(length):
        n = evidence['indel_evaluable_templates'][pos]
        k = evidence['indel_event_templates'][pos]
        if n and k * 10 >= n:
            start, stop = max(0, pos - 5), min(length, pos + evidence['indel_max_deleted_span'][pos] + 6)
            mask[start:stop] = b'\x01' * (stop - start)
            triggers.append({'kind': 'anchor', 'position_1based': pos + 1, 'events': k,
                             'denominator': n, 'mask_start0': start, 'mask_end0_exclusive': stop})
        total = column_depth(evidence, pos) + sum(evidence[x][pos] for x in ('ambiguous', 'deletion', 'conflict'))
        deleted = evidence['deletion'][pos]
        if total and deleted * 10 >= total:
            start, stop = max(0, pos - 5), min(length, pos + 6)
            mask[start:stop] = b'\x01' * (stop - start)
            triggers.append({'kind': 'deletion_site', 'position_1based': pos + 1, 'events': deleted,
                             'denominator': total, 'mask_start0': start, 'mask_end0_exclusive': stop})
    return mask, triggers


def evaluate_position(evidence, pos, length, indel_mask):
    counts = {base: evidence[base + '_fwd'][pos] + evidence[base + '_rev'][pos] for base in BASES}
    dp = sum(counts.values())
    total = dp + sum(evidence[x][pos] for x in ('ambiguous', 'deletion', 'conflict'))
    dominant = max(BASES, key=lambda base: counts[base])
    support = counts[dominant]
    reasons = []
    if dp < 10:
        reasons.append('ACGT_template_DP_below10')
    if not total or support * 10 < 9 * total:
        reasons.append('dominant_fraction_below0.90')
    if evidence[dominant + '_fwd'][pos] < 2 or evidence[dominant + '_rev'][pos] < 2:
        reasons.append('dominant_strand_support_below2')
    if indel_mask[pos]:
        reasons.append('indel_neighborhood')
    if pos < 200 or pos >= length - 200:
        reasons.append('linear_end_mask')
    called = 'N' if reasons else dominant
    if called != 'N':
        assert support > 0 and support == evidence[called + '_fwd'][pos] + evidence[called + '_rev'][pos]
        assert evidence[called + '_fwd'][pos] >= 2 and evidence[called + '_rev'][pos] >= 2
    return called, reasons, dp, total, dominant, support


def alignment_audit(bam, contig):
    counts = Counter()
    for read in bam.fetch(contig):
        counts['CP_alignment_records'] += 1
        for flag, name in [(0x100, 'secondary'), (0x800, 'supplementary'), (0x200, 'QCfailed'), (0x400, 'duplicate_flagged')]:
            if read.flag & flag:
                counts[name] += 1
        if reject_reason(read) is None:
            counts['eligible_alignment_records_before_base_filters'] += 1
        if read.is_paired and not read.is_proper_pair:
            counts['paired_not_proper'] += 1
    return dict(counts)


def competitive_reference_sequence(path, contig):
    pieces = []
    selected = False
    occurrences = 0
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                selected = line[1:].split()[0] == contig
                occurrences += selected
            elif selected:
                pieces.append(line.strip())
    if occurrences != 1:
        raise ValueError('Competitive FASTA must have exactly one named CP record')
    return ''.join(pieces).upper()


def preflight(args):
    import pysam
    seq, reference = verified_reference(args.reference)
    mapped_cp = competitive_reference_sequence(args.competitive_reference, args.contig)
    if mapped_cp != seq:
        raise ValueError('Actual competitive-reference CP sequence differs from bound source')
    mapref_sha = sha256_file(args.competitive_reference)
    source_manifest = args.bam.parent / 'REFERENCE_MANIFEST.json'
    rows = json.loads(source_manifest.read_text())
    matching = [row for row in rows if Path(row['path']).resolve() == args.competitive_reference.resolve()]
    if len(matching) != 1 or matching[0]['sha256'] != mapref_sha:
        raise ValueError('Competitive FASTA does not match the completed pilot reference manifest')
    pysam.quickcheck(str(args.bam))
    with pysam.AlignmentFile(str(args.bam), 'rb', require_index=True, threads=1) as bam:
        bam.check_index()
        if args.contig not in bam.references or bam.get_reference_length(args.contig) != len(seq):
            raise ValueError('BAM CP contig absent or length does not match verified source')
        header = bam.header.to_dict()
        if header.get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('Coordinate-sorted BAM required')
        if len(bam.references) < 2:
            raise ValueError('Competitive reference must include more than the CP contig')
        if str(args.competitive_reference) not in ' '.join(p.get('CL', '') for p in header.get('PG', [])):
            raise ValueError('BAM program header does not bind the requested competitive reference path')
        programs = [{'ID': p.get('ID'), 'PN': p.get('PN'), 'VN': p.get('VN'),
                     'command_sha256': hashlib.sha256(p.get('CL', '').encode()).hexdigest()}
                    for p in header.get('PG', [])]
        dup_programs = [p.get('ID') for p in header.get('PG', [])
                        if any(word in (p.get('PN', '') + ' ' + p.get('CL', '')).lower()
                               for word in ('markdup', 'markduplicates', 'rmdup', 'dedup'))]
    info = {'checked_utc': utcnow(), 'source_reference': reference, 'bam': str(args.bam.resolve()),
            'competitive_reference': str(args.competitive_reference), 'competitive_reference_sha256': mapref_sha,
            'pilot_reference_manifest_sha256': sha256_file(source_manifest),
            'coordinate_proof': {'path': str(args.coordinate_proof), 'sha256': sha256_file(args.coordinate_proof)} if args.coordinate_proof else 'Exact source and transformed sequence hashes independently verified; no external proof file supplied',
            'BAM_reference_contigs': len(header['SQ']), 'CP_contig': args.contig,
            'python': sys.version, 'pysam': pysam.__version__, 'samtools_htslib': pysam.__samtools_version__,
            'header_programs': programs, 'header_duplicate_related_program_IDs': dup_programs,
            'input_duplicate_status': 'Upstream declared not deduplicated; header absence/zero duplicate flags cannot prove independent molecules',
            'output_directory': str(args.output_dir), 'output_exists': args.output_dir.exists()}
    if info['output_exists']:
        raise FileExistsError('Output already exists; choose a new directory, never overwrite')
    return seq, info


def run(args):
    import pysam
    started = time.monotonic()
    limit = THRESHOLDS['address_space_limit_bytes']
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (min(limit, soft) if soft != resource.RLIM_INFINITY else limit, hard))
    seq, pre = preflight(args)
    if args.preflight_only:
        print(json.dumps(pre, indent=2))
        return
    args.output_dir.mkdir(parents=False, exist_ok=False)
    json_new(args.output_dir / 'thresholds.json', {'declared_before_read_processing': True, **THRESHOLDS})
    json_new(args.output_dir / 'preflight.json', pre)
    json_new(args.output_dir / 'STARTED.json', {'started_utc': utcnow(), 'argv': sys.argv})
    before = (args.bam.stat().st_size, args.bam.stat().st_mtime_ns)
    try:
        bam_hash = sha256_file(args.bam)
        with pysam.AlignmentFile(str(args.bam), 'rb', require_index=True, threads=1) as bam:
            read_audit = alignment_audit(bam, args.contig)
        with pysam.AlignmentFile(str(args.bam), 'rb', require_index=True, threads=1) as bam:
            evidence, filters = collect_evidence(bam, args.contig, len(seq))
        if before != (args.bam.stat().st_size, args.bam.stat().st_mtime_ns):
            raise RuntimeError('BAM changed while reading; no final sequence emitted')
        masks, triggers = make_masks(evidence, len(seq))
        output = []
        failure_counts = Counter()
        depth_bins = Counter()
        total_depth = 0
        called_ref = called_alt = mixed = all_observed_support = 0
        columns = ['position_1based', 'reference', 'call', 'ACGT_template_DP', 'fraction_denominator',
                   'dominant_base', 'dominant_support', 'dominant_fraction'] + list(evidence) + ['failure_reasons']
        with (args.output_dir / 'per_base.tsv').open('x') as handle:
            handle.write('\t'.join(columns) + '\n')
            for pos in range(len(seq)):
                call, reasons, dp, denom, dom, support = evaluate_position(evidence, pos, len(seq), masks)
                output.append(call)
                total_depth += dp
                key = '0' if dp == 0 else '1-4' if dp < 5 else '5-9' if dp < 10 else '10-19' if dp < 20 else '20-49' if dp < 50 else '50-99' if dp < 100 else '100-499' if dp < 500 else '500+'
                depth_bins[key] += 1
                failure_counts.update(reasons)
                if sum(evidence[b + '_fwd'][pos] + evidence[b + '_rev'][pos] > 0 for b in BASES) > 1:
                    mixed += 1
                if call != 'N':
                    called_ref += call == seq[pos]
                    called_alt += call != seq[pos]
                    all_observed_support += support > 0
                values = [pos + 1, seq[pos], call, dp, denom, dom, support,
                          f'{support/denom:.8f}' if denom else 'NA'] + [evidence[k][pos] for k in evidence] + [';'.join(reasons) or 'PASS']
                handle.write('\t'.join(map(str, values)) + '\n')
        callable_n = called_ref + called_alt
        assert all_observed_support == callable_n
        sequence = ''.join(output)
        with (args.output_dir / 'Hachi1_CP_reference_coordinate_partial.fasta').open('x') as handle:
            handle.write('>Hachi1_DRR351751_CP_reference_coordinate_PARTIAL_UNPHASED_unknowns_N\n')
            handle.write('\n'.join(sequence[i:i+80] for i in range(0, len(sequence), 80)) + '\n')
        normalized = normalized_path2(sequence)
        normalized_ref = normalized_path2(seq)
        assert normalized.count('N') == sequence.count('N')
        with (args.output_dir / 'Hachi1_CP_normalized_path2_PARTIAL.fasta').open('x') as handle:
            handle.write('>Hachi1_DRR351751_CP_normalized_path2_PARTIAL_UNPHASED_unknowns_N\n')
            handle.write('\n'.join(normalized[i:i+80] for i in range(0, len(normalized), 80)) + '\n')
        with (args.output_dir / 'normalized_path2_callability.tsv').open('x') as handle:
            handle.write('normalized_position_1based\tsource_path1_position_1based\treference\tcall\tcallable\tdominant_support_fwd_normalized\tdominant_support_rev_normalized\n')
            for pos in range(len(seq)):
                source_pos = normalized_source_position(pos)
                raw_call = sequence[source_pos]
                expected = reverse_complement(raw_call) if 109153 <= pos < 128097 else raw_call
                assert normalized[pos] == expected
                fwd = evidence[raw_call + '_fwd'][source_pos] if raw_call in BASES else 0
                rev = evidence[raw_call + '_rev'][source_pos] if raw_call in BASES else 0
                if 109153 <= pos < 128097:
                    fwd, rev = rev, fwd
                handle.write(f'{pos+1}\t{source_pos+1}\t{normalized_ref[pos]}\t{normalized[pos]}\t{int(normalized[pos] != "N")}\t{fwd}\t{rev}\n')
        json_new(args.output_dir / 'indel_risk_masks.json', triggers)
        outputs = {p.name: {'bytes': p.stat().st_size, 'sha256': sha256_file(p)}
                   for p in args.output_dir.iterdir() if p.is_file()}
        metrics = {'status': 'COMPLETE_PARTIAL_REFERENCE_COORDINATE_SEQUENCE', 'completed_utc': utcnow(),
            'elapsed_seconds': round(time.monotonic() - started, 3),
            'max_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            'BAM_sha256': bam_hash, 'BAM_bytes': before[0], 'reference_length_denominator': len(seq),
            'callable_positions': callable_n, 'unknown_N_positions': len(seq) - callable_n,
            'callable_fraction_of_full_reference': callable_n / len(seq),
            'interior_denominator_after_end_mask': len(seq) - 400,
            'callable_fraction_of_interior': callable_n / (len(seq) - 400),
            'called_reference_bases_with_observed_support': called_ref, 'called_nonreference_bases_with_observed_support': called_alt,
            'every_called_base_observed_support_assertion': all_observed_support == callable_n,
            'mean_ACGT_template_depth_full_reference': total_depth / len(seq), 'depth_bins_positions': dict(depth_bins),
            'failure_reason_position_counts_overlapping': dict(failure_counts), 'positions_with_multiple_observed_ACGT_bases': mixed,
            'indel_masked_positions': sum(masks), 'indel_trigger_count': len(triggers),
            'CP_alignment_audit': read_audit, 'filtered_alignment_observations_not_unique_reads': filters,
            'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
            'normalized_path2_sequence_sha256': hashlib.sha256(normalized.encode()).hexdigest(),
            'normalized_callability_preserved': normalized.count('N') == sequence.count('N'),
            'normalization_formula': THRESHOLDS['normalized_path2_transform'],
            'source_and_coordinate_preflight': pre, 'outputs_before_metrics': outputs,
            'limitations': [
                'Reference-coordinate partial unphased sequence; not a complete plastome or de novo assembly.',
                'Unknowns stay N; reference sequence is never used to supply a call.',
                'Circular junction is not validated; first and last 200 bp forced N.',
                'Both IR copies and nuclear CP-like regions create mapping ambiguity; MAPQ filtering causes missingness and cannot prove origin.',
                'No duplicate marking or molecular deduplication; unmarked PCR/optical duplicates may inflate support. Template collapse only handles overlapping mates.',
                'Original path1 and SSC-normalized path2 are separate outputs; exact source/target reference hashes and coordinate formula bind the conversion. Both IR copies retained.',
                'CIGAR indel anchors are not left-normalized; equivalent indels split across repeat anchors may evade 10% anchor threshold. Deletion-site risk masks partly reduce this.',
                'Thresholds establish engineering callability; no biological independence, phasing, heteroplasmy, taxon identity, or full-structure validation.'
            ]}
        json_new(args.output_dir / 'metrics.json', metrics)
        print(json.dumps({k: metrics[k] for k in ('status', 'callable_positions', 'reference_length_denominator', 'elapsed_seconds', 'max_RSS_MiB')}, indent=2))
    except Exception as error:
        json_new(args.output_dir / 'FAILED.json', {'status': 'FAILED_DO_NOT_USE_OUTPUTS', 'error': str(error), 'type': type(error).__name__, 'at_utc': utcnow()})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam', type=Path, default=Path('workspace/wgs/dt_research_p57_20260930/pilot_1m_qualityonly/pilot.bam'))
    parser.add_argument('--reference', type=Path, default=Path('workspace/wgs/organelle_20260914/inputs/D.tenuipes_chloroplast.fasta'))
    parser.add_argument('--contig', default='CP_hime2')
    parser.add_argument('--competitive-reference', type=Path, default=Path('workspace/wgs/dt_research_p57_20260930/pilot_1m/hime2_nuclear_plus_cp.fa'))
    parser.add_argument('--coordinate-proof', type=Path, help='Optional independently prepared path-provenance JSON to bind by SHA256; transformation is also verified against fixed target sequence hash')
    parser.add_argument('--output-dir', type=Path, default=Path('workspace/wgs/dt_research_p57_20260930/cp_callable_1m'))
    parser.add_argument('--preflight-only', action='store_true', help='Check reference, BAM header/index/EOF and output nonexistence; no read analysis or output writes')
    run(parser.parse_args())


if __name__ == '__main__':
    main()
