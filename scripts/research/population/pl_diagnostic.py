#!/usr/bin/env python3
"""Frozen-window PL availability diagnostic; no genotype changes or PCA.

Only indexed regions from the exact P1 BED are queried. PL is interpreted as
likelihood, never posterior probability; disagreement with GT is descriptive.
"""
import argparse
import bisect
import collections
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import struct
import sys
import time

EXPECTED_BED_SHA = '5fdffba9b3afe08abd5418caf6f56927df0e39955fdafc1bfff71f3312d4a288'
EXPECTED_HEADER_SHA = 'f9a99bc5d89c3b69463fd459f015dafa94d3d54284a592e166552651998566f9'
EXPECTED_FAI_SHA = '2c70d093718e9dbe8620a6af948a7c896cf5102b7c8e21717c79ba1ee50d9adb'
HISTORICAL_VCF_SHA = '43472edd28b406493d5745d71a1d7d403d2086d810e59d4f50194858a176adb9'
SAMPLES = 'DRR351754 DRR351755 DRR351756 AIC-F20-S13 AIC-F22-S15 AIC-F23-S16 AIC-M19-S12 AIC-M21-S14 AIC-M6-S11 HCO-F10-S8 HCO-F9-S7 HCO-M6-S5 HCO-M7-S6 IZU-F1-S9 IZU-F3-S10 IZU-F4-S1 IZU-F5-S2 IZU-F6-S3 IZU-F8-S4'.split()
EXCLUDED = {'DRR351754', 'DRR351755', 'DRR351756', 'AIC-F20-S13', 'AIC-M6-S11'}
MAIN14 = [s for s in SAMPLES if s not in EXCLUDED]
QUERY_FORMAT = '%CHROM\t%POS\t%REF\t%ALT\t%MQ\t%QD\t%FS\t%SOR[\t%GT:%AD:%DP:%GQ:%PL]\n'
PL_CLASSES = ('missing', 'wrong_cardinality', 'partial_missing', 'non_integer', 'negative_integer', 'three_nonnegative_min_nonzero', 'three_nonnegative_min_zero')
GT_CLASSES = ('missing_or_partial', 'not_diploid_encoded', 'not_biallelic_0_1', 'complete_0_0', 'complete_0_1', 'complete_1_1')
GQ_STRATA = ('missing', 'invalid', '0-19', '20-29', '>=30')
PRIMARY_REASONS = ('GT_incomplete_or_unsupported', 'DP_missing', 'DP_invalid', 'DP_outside_5_100', 'GQ_missing', 'GQ_invalid', 'GQ_below_20', 'eligible')


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def stat(path):
    s = Path(path).stat()
    return {'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns, 'inode': s.st_ino}


def atomic_json(path, obj):
    temp = Path(str(path) + '.tmp')
    with temp.open('w') as handle:
        json.dump(obj, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def parse_unsigned(raw):
    if raw in ('', '.'):
        return None, 'missing'
    if not re.fullmatch(r'\d+', raw):
        return None, 'invalid'
    return int(raw), 'valid'


def parse_gt(raw):
    alleles = re.split(r'[/|]', raw)
    if any(a in ('', '.') for a in alleles):
        return None, 'missing_or_partial'
    if len(alleles) != 2:
        return None, 'not_diploid_encoded'
    if any(a not in ('0', '1') for a in alleles):
        return None, 'not_biallelic_0_1'
    dosage = sum(map(int, alleles))
    return dosage, ('complete_0_0', 'complete_0_1', 'complete_1_1')[dosage]


def parse_pl(raw):
    """Do not repair malformed values or renormalize a nonzero minimum."""
    if raw in ('', '.'):
        return None, 'missing'
    fields = raw.split(',')
    if len(fields) != 3:
        return None, 'wrong_cardinality'
    if any(x in ('', '.') for x in fields):
        return None, 'partial_missing'
    if any(not re.fullmatch(r'-?\d+', x) for x in fields):
        return None, 'non_integer'
    values = tuple(map(int, fields))
    if min(values) < 0:
        return None, 'negative_integer'
    return values, 'three_nonnegative_min_zero' if min(values) == 0 else 'three_nonnegative_min_nonzero'


def parse_call(raw):
    f = raw.split(':')
    if len(f) != 5:
        raise ValueError('Expected five GT:AD:DP:GQ:PL fields')
    gt, ad, dp, gq, pl = f
    dosage, gt_class = parse_gt(gt)
    depth, dp_status = parse_unsigned(dp)
    quality, gq_status = parse_unsigned(gq)
    values, pl_class = parse_pl(pl)
    ad_values = [parse_unsigned(x) for x in ad.split(',')]
    ad_valid = len(ad_values) == 2 and all(s == 'valid' for _, s in ad_values)
    if dosage is None:
        primary_reason = 'GT_incomplete_or_unsupported'
    elif dp_status != 'valid':
        primary_reason = 'DP_' + dp_status
    elif not 5 <= depth <= 100:
        primary_reason = 'DP_outside_5_100'
    elif gq_status != 'valid':
        primary_reason = 'GQ_' + gq_status
    elif quality < 20:
        primary_reason = 'GQ_below_20'
    else:
        primary_reason = 'eligible'
    gq_stratum = gq_status if gq_status != 'valid' else ('0-19' if quality < 20 else '20-29' if quality < 30 else '>=30')
    return {'dosage': dosage, 'GT_class': gt_class, 'PL': values, 'PL_class': pl_class,
            'float32_preflight': float32_preflight(values) if values is not None else None,
            'GQ_stratum': gq_stratum, 'primary_reason': primary_reason, 'AD_valid_two_counts': ad_valid}


def float32_preflight(pl):
    """Arithmetic-only replay of storing first two normalized likelihood weights.

    These are NOT calibrated posterior probabilities. Actual compiled consumer
    arithmetic is not tested here. Keep the original integer PL unchanged.
    """
    m = min(pl)
    weights = [10.0 ** (-min(x - m, 10000) / 10.0) for x in pl]
    total = sum(weights)
    weights = [x / total for x in weights]
    first = [struct.unpack('f', struct.pack('f', x))[0] for x in weights[:2]]
    reconstructed = 1.0 - first[0] - first[1]
    return {'float32_reconstructed_L2_negative': reconstructed < 0,
            'float32_reconstructed_L2_zero_when_double_positive': reconstructed == 0 and weights[2] > 0,
            'normalized_L2_double_underflow': weights[2] == 0,
            'float32_reconstructed_L2': reconstructed}


class Moments:
    def __init__(self):
        self.n, self.total, self.squares, self.zero = 0, 0, 0, 0
        self.minimum, self.maximum = None, None

    def add(self, value):
        self.n += 1
        self.total += value
        self.squares += value * value
        self.zero += value == 0
        self.minimum = value if self.minimum is None else min(value, self.minimum)
        self.maximum = value if self.maximum is None else max(value, self.maximum)

    def output(self, prefix):
        variance = (self.squares - self.total * self.total / self.n) / (self.n - 1) if self.n > 1 else None
        return {prefix + '_n': self.n, prefix + '_sum': self.total, prefix + '_mean': self.total / self.n if self.n else None,
                prefix + '_sd': math.sqrt(max(0, variance)) if variance is not None else None,
                prefix + '_min': self.minimum, prefix + '_max': self.maximum, prefix + '_zero_n': self.zero}


class Stratum:
    def __init__(self):
        names = ['calls', 'AD_valid_two_counts', 'PL_usable', 'PL_unique_best', 'PL_tied_two_best', 'PL_all_three_tied',
                 'GT_complete_PL_usable', 'GT_complete_PL_unique', 'GT_matches_unique_best', 'GT_differs_unique_best',
                 'GT_complete_PL_tied', 'GT_in_tied_minima', 'GT_outside_tied_minima',
                 'float32_reconstructed_L2_negative', 'float32_reconstructed_L2_zero_when_double_positive', 'normalized_L2_double_underflow']
        self.c = collections.Counter({k: 0 for k in names})
        self.c.update({'PL_' + k: 0 for k in PL_CLASSES})
        self.c.update({'GT_' + k: 0 for k in GT_CLASSES})
        self.m = {k: Moments() for k in ('best_second_gap', 'PL_range', 'GT_to_best_gap')}

    def add(self, call):
        c, pl, gt = self.c, call['PL'], call['dosage']
        c['calls'] += 1
        c['PL_' + call['PL_class']] += 1
        c['GT_' + call['GT_class']] += 1
        c['AD_valid_two_counts'] += call['AD_valid_two_counts']
        if pl is None:
            return
        c['PL_usable'] += 1
        for key, value in call['float32_preflight'].items():
            if key != 'float32_reconstructed_L2':
                c[key] += value
        ordered = sorted(pl)
        minima = [i for i, x in enumerate(pl) if x == ordered[0]]
        c[{1: 'PL_unique_best', 2: 'PL_tied_two_best', 3: 'PL_all_three_tied'}[len(minima)]] += 1
        self.m['best_second_gap'].add(ordered[1] - ordered[0])
        self.m['PL_range'].add(ordered[2] - ordered[0])
        if gt is not None:
            c['GT_complete_PL_usable'] += 1
            self.m['GT_to_best_gap'].add(pl[gt] - ordered[0])
            if len(minima) == 1:
                c['GT_complete_PL_unique'] += 1
                c['GT_matches_unique_best' if gt == minima[0] else 'GT_differs_unique_best'] += 1
            else:
                c['GT_complete_PL_tied'] += 1
                c['GT_in_tied_minima' if gt in minima else 'GT_outside_tied_minima'] += 1

    def output(self):
        c = self.c
        assert sum(c['PL_' + k] for k in PL_CLASSES) == c['calls']
        assert sum(c['GT_' + k] for k in GT_CLASSES) == c['calls']
        assert c['PL_unique_best'] + c['PL_tied_two_best'] + c['PL_all_three_tied'] == c['PL_usable']
        assert c['GT_matches_unique_best'] + c['GT_differs_unique_best'] == c['GT_complete_PL_unique']
        assert c['GT_in_tied_minima'] + c['GT_outside_tied_minima'] == c['GT_complete_PL_tied']
        assert c['GT_complete_PL_unique'] + c['GT_complete_PL_tied'] == c['GT_complete_PL_usable']
        result = dict(c)
        for key, moments in self.m.items():
            result.update(moments.output(key))
        return result


def load_windows(bed, fai):
    reference = {}
    for line in fai.read_text().splitlines():
        f = line.split('\t')
        if len(f) != 5 or f[0] in reference or int(f[1]) <= 0:
            raise ValueError('Invalid FAI')
        reference[f[0]] = int(f[1])
    windows = collections.defaultdict(list)
    for line in bed.read_text().splitlines():
        f = line.split('\t')
        if len(f) != 3:
            raise ValueError('Invalid BED field count')
        chrom, start, end = f[0], int(f[1]), int(f[2])
        if chrom not in reference or not 0 <= start < end <= reference[chrom]:
            raise ValueError('BED outside reference bounds')
        windows[chrom].append((start, end))
    for v in windows.values():
        v.sort()
        if any(v[i][0] < v[i - 1][1] for i in range(1, len(v))):
            raise ValueError('Overlapping windows')
    return dict(windows), reference


class Analyzer:
    def __init__(self, samples, windows):
        self.samples, self.windows = samples, windows
        self.starts = {c: [s for s, _ in v] for c, v in windows.items()}
        self.rows, self.contigs = collections.Counter(), collections.Counter()
        self.seen, self.examples = set(), []
        self.example_counts = collections.Counter()
        self.strata = {(s, scope, gq): Stratum() for s in samples + ['ALL19', 'MAIN14'] for scope in ('all_DP_GQ', 'primary_DP5_100_GQ20') for gq in ('ALL',) + GQ_STRATA}
        self.primary = {s: collections.Counter({k: 0 for k in PRIMARY_REASONS}) for s in samples}

    def consume(self, line):
        f = line.rstrip('\n').split('\t')
        if len(f) != 8 + len(self.samples):
            raise ValueError('Wrong query field count')
        self.rows['queried_rows'] += 1
        chrom, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
        starts = self.starts.get(chrom, [])
        i = bisect.bisect_right(starts, pos - 1) - 1
        if i < 0 or not self.windows[chrom][i][0] <= pos - 1 < self.windows[chrom][i][1]:
            self.rows['outside_window_start'] += 1
            return
        if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
            self.rows['not_biallelic_acgt_SNP'] += 1
            return
        site = chrom, pos, ref, alt
        if site in self.seen:
            raise ValueError('Duplicate selected SNP')
        self.seen.add(site)
        self.rows['analyzed_SNP_rows'] += 1
        self.contigs[chrom] += 1
        for sample, raw in zip(self.samples, f[8:]):
            call = parse_call(raw)
            self.primary[sample][call['primary_reason']] += 1
            groups = [sample, 'ALL19'] + (['MAIN14'] if sample in MAIN14 else [])
            scopes = ['all_DP_GQ'] + (['primary_DP5_100_GQ20'] if call['primary_reason'] == 'eligible' else [])
            for group in groups:
                for scope in scopes:
                    for gq in ('ALL', call['GQ_stratum']):
                        self.strata[group, scope, gq].add(call)
            labels = [call['PL_class']]
            if call['PL'] is not None and call['dosage'] is not None:
                p = call['PL']
                if p.count(min(p)) == 1 and p[call['dosage']] != min(p):
                    labels.append('GT_differs_unique_PL_best')
                elif p.count(min(p)) > 1:
                    labels.append('PL_tied_minimum')
            for label in labels:
                if self.example_counts[label] < 3:
                    self.examples.append({'label': label, 'CHROM': chrom, 'POS': pos, 'REF': ref, 'ALT': alt, 'sample': sample,
                                          'GT_AD_DP_GQ_PL': raw, 'primary_reason': call['primary_reason']})
                    self.example_counts[label] += 1

    def export(self, out):
        rows = []
        for key, value in self.strata.items():
            group, scope, gq = key
            rows.append({'group_or_sample': group, 'scope': scope, 'GQ_stratum': gq, **value.output()})
        n = self.rows['analyzed_SNP_rows']
        for sample in self.samples:
            assert sum(self.primary[sample].values()) == n
            assert self.strata[sample, 'all_DP_GQ', 'ALL'].c['calls'] == n
            assert self.strata[sample, 'primary_DP5_100_GQ20', 'ALL'].c['calls'] == self.primary[sample]['eligible']
        for group in self.samples + ['ALL19', 'MAIN14']:
            for scope in ('all_DP_GQ', 'primary_DP5_100_GQ20'):
                total = self.strata[group, scope, 'ALL'].c
                assert all(sum(self.strata[group, scope, g].c[k] for g in GQ_STRATA) == total[k] for k in total)
        for group, members in [('ALL19', self.samples), ('MAIN14', [s for s in self.samples if s in MAIN14])]:
            for scope in ('all_DP_GQ', 'primary_DP5_100_GQ20'):
                total = self.strata[group, scope, 'ALL'].c
                assert all(sum(self.strata[s, scope, 'ALL'].c[k] for s in members) == total[k] for k in total)
        write_tsv(out / 'PL_strata.tsv', rows)
        write_tsv(out / 'primary_scope_denominators.tsv', [{'sample': s, 'main14': s in MAIN14, 'all_window_SNPs': n, **dict(self.primary[s])} for s in self.samples])
        write_tsv(out / 'examples.tsv', self.examples)
        write_tsv(out / 'contig_counts.tsv', [{'CHROM': c, 'SNPs': n} for c, n in self.contigs.items()])
        result = {'status': 'PL_DIAGNOSTIC_COMPLETE', 'rows': dict(self.rows), 'samples': self.samples, 'main14': MAIN14,
                  'strata_rows': len(rows), 'denominator_checks_pass': True,
                  'overall': [r for r in rows if r['group_or_sample'] in ('ALL19', 'MAIN14') and r['GQ_stratum'] == 'ALL'],
                  'definitions': {'PL_usable': 'Exactly three nonnegative integers for diploid biallelic SNPs, order 0/0, 0/1, 1/1; nonzero minima retained and separately counted.',
                                  'all_DP_GQ': 'Every selected SNP/sample, including incomplete GT and missing/invalid DP/GQ/PL.',
                                  'primary_DP5_100_GQ20': 'Complete diploid biallelic GT, valid DP5-100 inclusive, valid GQ>=20; AD and PL are not entry requirements.',
                                  'unique_comparison_denominator': 'GT_complete_PL_unique; ties are excluded and separately counted, not silently resolved.',
                                  'GT_to_best_gap': 'PL at called GT minus smallest PL; no posterior calculation.',
                                  'float32_preflight': 'From valid triad, normalize 10^(-(PL-minPL)/10) likelihood weights in double precision, cast first two to IEEE float32 with struct, reconstruct third as double 1-f32(L0)-f32(L1). Counters use PL_usable denominator; negative and exactly-zero-despite-positive-double-L2 are separate. Original PL unchanged. This is arithmetic preflight, not compiled PCAngsd validation; normalized weights are not calibrated posteriors.',
                                  'missing_values': 'Numeric zero is a measured zero; JSON null and blank TSV numeric cells mean not estimable (n=0, or SD n<2).',
                                  'primary_exclusions': 'Mutually exclusive first failure in GT, DP, then GQ order. All-data PL summaries are not conditioned on eligibility.'},
                  'limitations': ['Existing ascertained SNP panel, 998 windows; not all-site or whole-genome validation.',
                                  'PL is conditional on upstream read/model/reference assumptions; no independent mapping, ploidy or calibration test.',
                                  'GT need not equal the maximum-likelihood genotype when calling uses priors; disagreement is not automatically error.',
                                  'Likelihood differences and ties are descriptive; no new PL threshold, posterior normalization, significance test or genotype alteration.',
                                  'Per-call counts are dependent across linked loci and samples; they are not independent biological replicates.']}
        atomic_json(out / 'SUMMARY.json', result)
        return result


def write_tsv(path, rows):
    if not rows:
        raise ValueError('Refuse schema-less empty TSV')
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {'status': 'RUNNING', 'started_utc': utc(), 'commands': [], 'script_sha256': sha_file(__file__),
               'historical_VCF_sha256': HISTORICAL_VCF_SHA, 'historical_VCF_hash_date': '2026-09-14', 'large_VCF_body_rehashed': False,
               'paths': {k: str(getattr(args, k).resolve()) for k in ('vcf', 'bed', 'fai')},
               'resource_scope': 'Single indexed BED query plus header read; no full VCF scan, no sample/GT changes',
               'runtime': {'python': sys.version, 'executable': sys.executable, 'executable_sha256': sha_file(sys.executable)}}
    atomic_json(out / 'RUN_RECEIPT.json', receipt)
    process = None
    try:
        before = stat(args.vcf)
        receipt['VCF_stat_before'] = before
        receipt['inputs'] = {}
        for key, path, expected in [('BED', args.bed, EXPECTED_BED_SHA), ('FAI', args.fai, EXPECTED_FAI_SHA)]:
            receipt['inputs'][key] = {'path': str(path.resolve()), 'stat': stat(path), 'sha256': sha_file(path)}
            if receipt['inputs'][key]['sha256'] != expected:
                raise ValueError(key + ' differs from frozen P1 input')
        indexes = [Path(str(args.vcf) + suffix) for suffix in ('.tbi', '.csi') if Path(str(args.vcf) + suffix).is_file()]
        if not indexes:
            raise ValueError('Indexed-window query requires existing VCF index')
        receipt['indexes_before'] = {str(p): {'stat': stat(p), 'sha256': sha_file(p)} for p in indexes}
        windows, reference = load_windows(args.bed, args.fai)
        count = sum(map(len, windows.values()))
        bases = sum(e - s for v in windows.values() for s, e in v)
        receipt['window_count'], receipt['window_reference_bases'] = count, bases
        if (count, bases) != (998, 996926):
            raise ValueError('Unexpected window count/base denominator')
        header_cmd = [args.bcftools, 'view', '--no-version', '-h', str(args.vcf)]
        receipt['commands'].append(header_cmd)
        header = subprocess.check_output(header_cmd)
        (out / 'source.header.vcf').write_bytes(header)
        receipt['header_sha256'] = hashlib.sha256(header).hexdigest()
        if receipt['header_sha256'] != EXPECTED_HEADER_SHA:
            raise ValueError('Source header differs from frozen P1 input')
        text = header.decode()
        samples = next(l for l in text.splitlines() if l.startswith('#CHROM\t')).split('\t')[9:]
        if samples != SAMPLES:
            raise ValueError('Source sample order/identity differs')
        contigs = re.findall(r'^##contig=<ID=([^,]+),length=(\d+)', text, re.M)
        if len(contigs) != len(reference) or {c: int(n) for c, n in contigs} != reference:
            raise ValueError('Source header and FAI dictionaries differ')
        pl_header = next(l for l in text.splitlines() if l.startswith('##FORMAT=<ID=PL,'))
        receipt['PL_header'] = pl_header
        if not pl_header.startswith('##FORMAT=<ID=PL,Number=G,Type=Integer,'):
            raise ValueError('Unexpected PL definition')
        version_cmd = [args.bcftools, '--version']
        receipt['commands'].append(version_cmd)
        receipt['runtime']['bcftools_version'] = subprocess.check_output(version_cmd, text=True).strip()
        tool_path = subprocess.check_output(['which', args.bcftools], text=True).strip()
        receipt['runtime']['bcftools_path'] = tool_path
        receipt['runtime']['bcftools_executable_sha256'] = sha_file(tool_path)
        query = [args.bcftools, 'query', '-R', str(args.bed), '-f', QUERY_FORMAT, str(args.vcf)]
        receipt['commands'].append(query)
        analyzer = Analyzer(samples, windows)
        projection = hashlib.sha256()
        projection.update(('#CHROM\tPOS\tREF\tALT\tMQ\tQD\tFS\tSOR\t' + '\t'.join(s + '[GT:AD:DP:GQ]' for s in samples) + '\n').encode())
        with (out / 'query.stderr').open('w') as errors, gzip.open(out / 'window_calls_PL.tsv.gz', 'wt', compresslevel=6) as raw:
            raw.write('#CHROM\tPOS\tREF\tALT\tMQ\tQD\tFS\tSOR\t' + '\t'.join(s + '[GT:AD:DP:GQ:PL]' for s in samples) + '\n')
            process = subprocess.Popen(query, stdout=subprocess.PIPE, stderr=errors, text=True)
            for line in process.stdout:
                raw.write(line)
                analyzer.consume(line)
                fields = line.rstrip('\n').split('\t')
                projection.update(('\t'.join(fields[:8] + [s.rsplit(':', 1)[0] for s in fields[8:]]) + '\n').encode())
            process.stdout.close()
            receipt['query_exit_code'] = process.wait()
        if receipt['query_exit_code'] != 0:
            raise RuntimeError('bcftools query failed')
        if analyzer.rows['analyzed_SNP_rows'] != 50043:
            raise ValueError('SNP count differs from P1 frozen panel')
        receipt['P1_fields_projected_uncompressed_sha256'] = projection.hexdigest()
        if args.p1_raw:
            old_hash = hashlib.sha256()
            with gzip.open(args.p1_raw, 'rb') as old:
                for chunk in iter(lambda: old.read(1024 * 1024), b''):
                    old_hash.update(chunk)
            receipt['P1_raw_comparison'] = {'path': str(args.p1_raw), 'compressed_sha256': sha_file(args.p1_raw),
                                            'uncompressed_sha256': old_hash.hexdigest(), 'fields_identical': old_hash.hexdigest() == projection.hexdigest()}
            if not receipt['P1_raw_comparison']['fields_identical']:
                raise ValueError('Existing P1 fields are not byte-identical after removing PL')
        summary = analyzer.export(out)
        receipt['VCF_stat_after'] = stat(args.vcf)
        receipt['indexes_after'] = {str(p): {'stat': stat(p), 'sha256': sha_file(p)} for p in indexes}
        receipt['source_stat_unchanged'] = before == receipt['VCF_stat_after'] and receipt['indexes_before'] == receipt['indexes_after']
        if not receipt['source_stat_unchanged']:
            raise RuntimeError('Source VCF/index changed during run')
        receipt['analyzed_SNP_rows'] = summary['rows']['analyzed_SNP_rows']
        receipt['all19_calls'] = summary['overall'][0]['calls']
        receipt['output_files'] = {p.name: {'bytes': p.stat().st_size, 'sha256': sha_file(p), 'server_only_raw': p.name == 'window_calls_PL.tsv.gz'} for p in out.iterdir() if p.is_file() and p.name != 'RUN_RECEIPT.json'}
        receipt['status'] = 'PL_DIAGNOSTIC_COMPLETE'
        return summary
    except Exception as error:
        receipt['status'] = 'FAILED'
        receipt['error'] = type(error).__name__ + ': ' + str(error)
        if process is not None and process.poll() is None:
            process.terminate()
            process.wait()
        raise
    finally:
        receipt['finished_utc'] = utc()
        receipt['elapsed_seconds'] = time.monotonic() - start
        receipt['script_exit_code'] = 0 if receipt['status'] == 'PL_DIAGNOSTIC_COMPLETE' else 1
        atomic_json(out / 'RUN_RECEIPT.json', receipt)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf', type=Path, default=Path('workspace/wgs/population_audit_20260914/inputs/Dtenuipes_population.vcf.gz'))
    p.add_argument('--bed', type=Path, default=Path('workspace/wgs/dt_research_validation_20260924/allele_balance/windows.bed'))
    p.add_argument('--fai', type=Path, default=Path('workspace/wgs/genome.fa.fai'))
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--p1-raw', type=Path)
    p.add_argument('--bcftools', default='bcftools')
    a = p.parse_args()
    r = run(a)
    print(json.dumps({'status': r['status'], 'rows': r['rows'], 'out': str(a.out)}, sort_keys=True))
