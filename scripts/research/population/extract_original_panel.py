#!/usr/bin/env python3
"""Bounded original 37,566-marker PL extraction and BED identity verification.

No PCA, no PL normalization, no genotype repair, no original input writes.
Requires reviewed panel_identity.py next to this file. Run only after review.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import panel_identity as p

FORMAT = '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT:%AD:%DP:%GQ:%PL]\n'
REASONS = ('GT_missing_or_partial', 'DP_missing', 'DP_outside_5_100', 'GQ_missing', 'GQ_below_20', 'baseline_eligible')
PL_STATUSES = ('whole_missing','wrong_cardinality','partial_missing','non_integer','negative_integer','minimum_nonzero','normalized_integer_triad')
COUNTS = ('source_cells', 'source_GT_complete', 'source_GT_missing', 'baseline_called', 'baseline_missing',
          'PL_whole_missing', 'PL_partial_missing', 'PL_wrong_cardinality', 'PL_non_integer', 'PL_negative_integer', 'PL_minimum_nonzero', 'PL_normalized_integer_triad',
          'source_complete_GT_PL_unique', 'source_GT_matches_unique_PL', 'source_GT_differs_unique_PL',
          'source_complete_GT_PL_tied', 'source_GT_in_tied_PL_minima', 'source_GT_outside_tied_PL_minima',
          'PL_unique_minimum', 'PL_two_minima', 'PL_all_three_minima',
          'baseline_called_PL_unavailable_or_nonnormalized', 'candidate_common_called',
          'baseline_missing_PL_normalized_not_restored')


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def stat(path):
    s = Path(path).stat()
    return {'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns, 'inode': s.st_ino}


def atomic_json(path, data):
    t = Path(str(path) + '.tmp')
    with t.open('w') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(t, path)


def write_tsv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def classify_pl(raw):
    """Describe every PL vector without repairing or silently dropping a cell."""
    if raw in ('', '.'):
        return None, 'whole_missing'
    f = raw.split(',')
    if len(f) != 3:
        return None, 'wrong_cardinality'
    if any(x in ('', '.') for x in f):
        return None, 'partial_missing'
    if any(not re.fullmatch(r'-?\d+', x) for x in f):
        return None, 'non_integer'
    pl = tuple(map(int, f))
    if min(pl) < 0:
        return None, 'negative_integer'
    return pl, 'minimum_nonzero' if min(pl) else 'normalized_integer_triad'


class PanelCollector:
    def __init__(self, markers, bed, samples):
        self.markers, self.bed, self.samples = markers, bed, samples
        self.width = (len(samples) + 3) // 4
        if bed[:3] != b'\x6c\x1b\x01' or len(bed) != 3 + self.width * len(markers):
            raise ValueError('Wrong BED magic or size')
        self.coordinates = {(m['CHROM'], m['POS']): m for m in markers}
        if len(self.coordinates) != len(markers):
            raise ValueError('Duplicate panel coordinate')
        self.raw_rows = [None] * len(markers)
        self.reconstructed = bytearray(b'\x6c\x1b\x01' + bytes(self.width * len(markers)))
        self.baseline_mask = bytearray(len(markers) * len(samples))
        self.common_mask = bytearray(len(self.baseline_mask))
        self.calls = {s: collections.Counter({k: 0 for k in COUNTS + REASONS}) for s in samples}
        self.query_rows, self.extra_overlap_rows, self.selected_rows = 0, 0, 0

    def consume(self, line):
        f = line.rstrip('\n').split('\t')
        if len(f) != 5 + len(self.samples):
            raise ValueError('Unexpected query column count')
        self.query_rows += 1
        coordinate = f[0], int(f[1])
        if coordinate not in self.coordinates:
            self.extra_overlap_rows += 1
            return
        m = self.coordinates[coordinate]
        i = m['index']
        if self.raw_rows[i] is not None:
            raise ValueError('Duplicate VCF record at selected coordinate')
        ref, alt = f[2], f[3]
        if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
            raise ValueError('Source marker is not a distinct biallelic ACGT SNP')
        if {ref, alt} != {m['A1'], m['A2']}:
            raise ValueError('Source REF/ALT and original BIM allele set differ; no strand rescue')
        m['source_REF'], m['source_ALT'] = ref, alt
        if f[4] != 'PASS':
            raise ValueError('Original PCA marker no longer PASS')
        packed = self.bed[3 + self.width*i:3 + self.width*(i+1)]
        actual_codes = p.bed_codes(packed, len(self.samples))
        rebuilt = bytearray(self.width)
        for j, (sample, raw) in enumerate(zip(self.samples, f[5:])):
            gt_dosage, masked_dosage, reason = p.gt_and_quality(raw)
            expected_code = p.dosage_to_bed(masked_dosage, m['A1'], m['A2'], ref, alt)
            if expected_code != actual_codes[j]:
                raise ValueError('Original BED differs from source GT plus DP5-100/GQ20 at ' + str(coordinate) + ' sample ' + sample)
            rebuilt[j // 4] |= expected_code << (2*(j % 4))
            pl, pl_status = classify_pl(raw.rsplit(':', 1)[1])
            c = self.calls[sample]
            c['source_cells'] += 1
            c['source_GT_complete' if gt_dosage is not None else 'source_GT_missing'] += 1
            c[reason] += 1
            c['PL_' + pl_status] += 1
            baseline_called = expected_code != 1
            common_called = baseline_called and pl_status == 'normalized_integer_triad'
            c['baseline_called' if baseline_called else 'baseline_missing'] += 1
            c['candidate_common_called'] += common_called
            c['baseline_called_PL_unavailable_or_nonnormalized'] += baseline_called and not common_called
            c['baseline_missing_PL_normalized_not_restored'] += not baseline_called and pl_status == 'normalized_integer_triad'
            self.baseline_mask[i*len(self.samples)+j] = baseline_called
            self.common_mask[i*len(self.samples)+j] = common_called
            if pl is not None:
                minima = [k for k, value in enumerate(pl) if value == min(pl)]
                c[{1:'PL_unique_minimum', 2:'PL_two_minima', 3:'PL_all_three_minima'}[len(minima)]] += 1
                if gt_dosage is not None:
                    if len(minima) == 1:
                        c['source_complete_GT_PL_unique'] += 1
                        c['source_GT_matches_unique_PL' if gt_dosage in minima else 'source_GT_differs_unique_PL'] += 1
                    else:
                        c['source_complete_GT_PL_tied'] += 1
                        c['source_GT_in_tied_PL_minima' if gt_dosage in minima else 'source_GT_outside_tied_PL_minima'] += 1
        self.reconstructed[3+self.width*i:3+self.width*(i+1)] = rebuilt
        self.raw_rows[i] = line
        self.selected_rows += 1

    def finish(self):
        if self.selected_rows != len(self.markers) or any(x is None for x in self.raw_rows):
            raise ValueError('Not every original PCA marker was retrieved exactly once')
        if bytes(self.reconstructed) != self.bed:
            raise ValueError('Reconstructed quality-masked BED differs byte for byte')
        for c in self.calls.values():
            if sum(c[k] for k in REASONS) != len(self.markers):
                raise ValueError('Quality reason denominator failed')
            if c['baseline_called'] + c['baseline_missing'] != c['source_cells']:
                raise ValueError('BED called/missing denominator failed')
            if c['candidate_common_called'] + c['baseline_called_PL_unavailable_or_nonnormalized'] != c['baseline_called']:
                raise ValueError('Common-mask denominator failed')
            if sum(c['PL_'+x] for x in PL_STATUSES) != c['source_cells']:
                raise ValueError('PL denominator failed')
        if any(c > b for c, b in zip(self.common_mask, self.baseline_mask)):
            raise ValueError('Candidate mask restored an original missing GT')

    def export(self, out, header):
        self.finish()
        with gzip.open(out / 'panel_calls_original_order.tsv.gz', 'wt') as f:
            f.write(header)
            f.writelines(self.raw_rows)
        (out / 'baseline_called_mask.bin').write_bytes(self.baseline_mask)
        (out / 'candidate_common_mask.bin').write_bytes(self.common_mask)
        write_tsv(out / 'panel_markers.tsv', self.markers)
        write_tsv(out / 'sample_order.tsv', [{'index0': i, 'sample': s, 'location': s.split('-')[0]} for i, s in enumerate(self.samples)])
        sample_rows = [{'sample': s, 'location': s.split('-')[0], **dict(self.calls[s])} for s in self.samples]
        write_tsv(out / 'sample_denominators.tsv', sample_rows)
        group_rows = []
        for group in ['AIC','HCO','IZU','MAIN14']:
            members = [s for s in self.samples if group == 'MAIN14' or s.split('-')[0] == group]
            if members:
                group_rows.append({'group': group, 'samples': len(members), **{k:sum(self.calls[s][k] for s in members) for k in COUNTS+REASONS}})
        write_tsv(out / 'group_denominators.tsv', group_rows)
        review_flags = []
        for key in ('PL_wrong_cardinality','PL_non_integer','PL_negative_integer','PL_minimum_nonzero','source_GT_differs_unique_PL','source_GT_outside_tied_PL_minima'):
            if sum(c[key] for c in self.calls.values()):
                review_flags.append(key)
        summary = {'status':'EXTRACTED_PANEL_IDENTITY_VERIFIED', 'selected_markers':len(self.markers), 'samples':self.samples,
                   'cells':len(self.markers)*len(self.samples), 'query_rows':self.query_rows, 'extra_overlap_rows_not_selected':self.extra_overlap_rows,
                   'BED_cell_comparisons':len(self.markers)*len(self.samples), 'BED_reconstructed_sha256':hashlib.sha256(self.reconstructed).hexdigest(),
                   'source_REF_ALT_bound_at_every_marker':True,'groups':group_rows,'review_flags':review_flags,
                   'actual_source_allele_orientation':dict(collections.Counter('A1_is_ALT' if m['A1']==m['source_ALT'] else 'A1_is_REF' for m in self.markers)),
                   'historical_ID_direction_relative_to_actual_source':dict(collections.Counter('matches_REF_ALT' if m['historical_ID_alleles']==m['source_REF']+','+m['source_ALT'] else 'reversed_from_REF_ALT' for m in self.markers)),
                   'PL_conversion_gate':'BLOCKED_INVALID_OR_NONNORMALIZED_PL' if any(x.startswith('PL_') for x in review_flags) else 'NOT_EVALUATED_NUMERICAL_METHOD',
                   'not_a_PCA_readiness_approval':True,
                   'definitions':{'row_order':'Original BIM order; source raw query order saved separately.',
                                  'baseline_GT':'Original BED code: 00 A1/A1, 01 missing, 10 heterozygous, 11 A2/A2. Every cell checked against source GT plus original DP5-100/GQ20.',
                                  'PL':'Unmodified source values in actual source VCF REF/ALT order 0/0,0/1,1/1; never swapped to match a BIM ID. Whole/partial missing and malformed cases retained and counted; no equal-likelihood conversion here. Any illegal/nonnormalized vector blocks downstream conversion pending review.',
                                  'historical_ID_alleles':'Historical marker-name allele labels, not a REF/ALT orientation source. source_REF/source_ALT bind actual VCF fields; original BIM allele set and per-cell BED reconstruction remain strict.',
                                  'candidate_common_mask':'Baseline called AND complete three nonnegative integer PL entries with minimum0. Descriptive candidate; no PCA, no site deletion, no restored GTs.',
                                  'mask_layout':'One byte (0/1) per cell, variant-major original BIM order, inner original FAM order. Hash and dimensions in receipt.',
                                  'GT_PL_concordance':'Source complete GT versus raw triad minima, including nonzero-min triads separately counted. Disagreement alone is not a GT error.',
                                  'quality_reasons':'Sequential first failure GT then DP then GQ; AD/PL do not determine original BED mask.'},
                   'limits':['Variant-only ascertained, LD-pruned original PCA panel; not callable whole-genome denominator.',
                             'Internal identity/replay checks, not genotype truth, calibration or reference-bias validation.',
                             'No PL normalization, float32 conversion, numerical repair, frequency model, imputation or PCA.',
                             'No proof of full per-sample upstream PL-processing history; historical provenance gaps remain.']}
        atomic_json(out / 'SUMMARY.json', summary)
        return summary


def run(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    receipt = {'status':'RUNNING','started_utc':utc(),'commands':[], 'script_sha256':p.sha(__file__),
               'helper_sha256':p.sha(Path(__file__).with_name('panel_identity.py')), 'python':sys.version,'python_executable':sys.executable,
               'historical_source_VCF_sha256':p.HISTORICAL_VCF_SHA,'historical_source_hash_date':'2026-09-14','VCF_body_rehashed_current':False,
               'scope':'One indexed query at exactly original PCA positions; no full VCF scan and no PCA', 'output_directory':str(out)}
    atomic_json(out / 'EXTRACTION_RECEIPT.json', receipt)
    process = None
    try:
        source_files = {label:{e:Path(str(prefix)+'.'+e) for e in ('bed','bim','fam')} for label,prefix in [('original',args.original_prefix),('refiltered',args.refiltered_prefix)]}
        receipt['panel_inputs'] = {}
        for label, files in source_files.items():
            receipt['panel_inputs'][label] = {e:{'path':str(f),'sha256':p.sha(f),'stat':stat(f)} for e,f in files.items()}
            if any(receipt['panel_inputs'][label][e]['sha256'] != p.EXPECTED[label][e] for e in files):
                raise ValueError('Panel files changed since frozen P3 and reviewed preflight')
        original, refiltered = source_files['original'], source_files['refiltered']
        samples = p.parse_fam(original['fam'].read_text())
        if original['fam'].read_bytes() != refiltered['fam'].read_bytes():
            raise ValueError('Original and refiltered FAM differ')
        markers = p.bind_bims(original['bim'].read_text(), refiltered['bim'].read_text())
        if len(markers) != 37566:
            raise ValueError('Expected original 37566 markers')
        before = stat(args.vcf)
        receipt['VCF_path'], receipt['VCF_stat_before'] = str(args.vcf), before
        indexes = [Path(str(args.vcf)+s) for s in ('.tbi','.csi') if Path(str(args.vcf)+s).is_file()]
        if not indexes:
            raise ValueError('No source index')
        receipt['indexes_before'] = {str(x):{'stat':stat(x),'sha256':p.sha(x)} for x in indexes}
        fai_bytes = args.fai.read_bytes()
        if hashlib.sha256(fai_bytes).hexdigest() != p.FAI_SHA:
            raise ValueError('Reference FAI changed')
        reference = {x.split('\t')[0]:int(x.split('\t')[1]) for x in fai_bytes.decode().splitlines()}
        ref_order = {name:i for i,name in enumerate(reference)}
        command = [args.bcftools,'view','--no-version','-h',str(args.vcf)]
        receipt['commands'].append(command)
        header_bytes = subprocess.check_output(command)
        receipt['source_header_sha256'] = hashlib.sha256(header_bytes).hexdigest()
        if receipt['source_header_sha256'] != p.HEADER_SHA:
            raise ValueError('Source header changed')
        (out / 'source.header.vcf').write_bytes(header_bytes)
        header_text = header_bytes.decode()
        source_samples = next(x for x in header_text.splitlines() if x.startswith('#CHROM\t')).split('\t')[9:]
        if source_samples != p.SOURCE_SAMPLES or [x for x in source_samples if x in samples] != samples:
            raise ValueError('Source sample identity/order differs')
        contigs = re.findall(r'^##contig=<ID=([^,]+),length=(\d+)',header_text,re.M)
        if len(contigs) != len(reference) or {c:int(n) for c,n in contigs} != reference:
            raise ValueError('Header/FAI dictionary differs')
        for m in markers:
            if m['CHROM'] not in reference or not 1 <= m['POS'] <= reference[m['CHROM']]:
                raise ValueError('Original marker outside FAI bounds')
        regions = out / 'original_positions.bed'
        regions.write_text(''.join(m['CHROM']+'\t'+str(m['POS']-1)+'\t'+str(m['POS'])+'\n' for m in sorted(markers,key=lambda x:(ref_order[x['CHROM']],x['POS']))))
        samplefile = out / 'original_samples.txt'
        samplefile.write_text('\n'.join(samples)+'\n')
        receipt['generated_query_inputs'] = {x.name:{'sha256':p.sha(x),'bytes':x.stat().st_size} for x in (regions,samplefile)}
        receipt['FAI'] = {'path':str(args.fai),'sha256':p.FAI_SHA,'stat':stat(args.fai)}
        command = [args.bcftools,'--version']; receipt['commands'].append(command)
        receipt['bcftools_version'] = subprocess.check_output(command,text=True).strip()
        receipt['bcftools_executable'] = {'path':args.bcftools,'sha256':p.sha(args.bcftools)}
        query = [args.bcftools,'query','-R',str(regions),'-S',str(samplefile),'-f',FORMAT,str(args.vcf)]
        receipt['commands'].append(query)
        panel = PanelCollector(markers, original['bed'].read_bytes(), samples)
        raw_header = '#CHROM\tPOS\tREF\tALT\tFILTER\t'+'\t'.join(s+'[GT:AD:DP:GQ:PL]' for s in samples)+'\n'
        with (out / 'query.stderr').open('w') as err, gzip.open(out / 'query_order_calls.tsv.gz','wt') as raw:
            raw.write(raw_header)
            process = subprocess.Popen(query,stdout=subprocess.PIPE,stderr=err,text=True)
            for line in process.stdout:
                raw.write(line)
                panel.consume(line)
                if panel.selected_rows and panel.selected_rows % 5000 == 0:
                    receipt['progress'] = {'selected':panel.selected_rows,'queried':panel.query_rows,'elapsed_seconds':time.monotonic()-started}
                    atomic_json(out / 'EXTRACTION_RECEIPT.json',receipt)
            process.stdout.close()
            receipt['query_exit_code'] = process.wait()
        if receipt['query_exit_code']:
            raise RuntimeError('bcftools query failed')
        summary = panel.export(out,raw_header)
        receipt['VCF_stat_after'] = stat(args.vcf)
        receipt['indexes_after'] = {str(x):{'stat':stat(x),'sha256':p.sha(x)} for x in indexes}
        receipt['panel_inputs_sha256_after'] = {label:{e:p.sha(f) for e,f in files.items()} for label,files in source_files.items()}
        if before != receipt['VCF_stat_after'] or receipt['indexes_before'] != receipt['indexes_after']:
            raise ValueError('Source VCF/index changed during query')
        if any(receipt['panel_inputs_sha256_after'][label] != p.EXPECTED[label] for label in source_files):
            raise ValueError('Original/refiltered panel changed during query')
        receipt['original_inputs_unchanged'] = True
        receipt['status'] = summary['status']
        receipt['selected_markers'],receipt['BED_cell_comparisons'] = summary['selected_markers'],summary['BED_cell_comparisons']
        receipt['review_flags'] = summary['review_flags']
        receipt['output_files'] = {x.name:{'sha256':p.sha(x),'bytes':x.stat().st_size,'exclude_from_Git':x.suffix in ('.gz','.bin')} for x in out.iterdir() if x.is_file() and x.name != 'EXTRACTION_RECEIPT.json'}
        return summary
    except Exception as error:
        receipt['status'],receipt['error'] = 'FAILED',type(error).__name__+': '+str(error)
        if process is not None and process.poll() is None:
            process.terminate();process.wait()
        raise
    finally:
        receipt['finished_utc'],receipt['elapsed_seconds'] = utc(),time.monotonic()-started
        receipt['script_exit_code'] = 0 if receipt['status']=='EXTRACTED_PANEL_IDENTITY_VERIFIED' else 1
        atomic_json(out / 'EXTRACTION_RECEIPT.json',receipt)


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--out',type=Path,required=True)
    cli.add_argument('--original-prefix',type=Path,default=Path(p.ORIGINAL_PREFIX))
    cli.add_argument('--refiltered-prefix',type=Path,default=Path(p.REFILTERED_PREFIX))
    cli.add_argument('--vcf',type=Path,default=Path(p.VCF))
    cli.add_argument('--fai',type=Path,default=Path(p.FAI))
    cli.add_argument('--bcftools',default='/usr/bin/bcftools')
    args = cli.parse_args()
    result = run(args)
    print(json.dumps({'status':result['status'],'markers':result['selected_markers'],'cells':result['cells'],'review_flags':result['review_flags']}))
