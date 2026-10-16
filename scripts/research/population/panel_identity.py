#!/usr/bin/env python3
"""Small exact-panel identity helpers. No PL conversion, PCA, or source edits."""
import hashlib
from pathlib import Path
import re

SAMPLES = 'AIC-F22-S15 AIC-F23-S16 AIC-M19-S12 AIC-M21-S14 HCO-F10-S8 HCO-F9-S7 HCO-M6-S5 HCO-M7-S6 IZU-F1-S9 IZU-F3-S10 IZU-F4-S1 IZU-F5-S2 IZU-F6-S3 IZU-F8-S4'.split()
SOURCE_SAMPLES = 'DRR351754 DRR351755 DRR351756 AIC-F20-S13 AIC-F22-S15 AIC-F23-S16 AIC-M19-S12 AIC-M21-S14 AIC-M6-S11 HCO-F10-S8 HCO-F9-S7 HCO-M6-S5 HCO-M7-S6 IZU-F1-S9 IZU-F3-S10 IZU-F4-S1 IZU-F5-S2 IZU-F6-S3 IZU-F8-S4'.split()
ORIGINAL_PREFIX = 'workspace/wgs/population_analysis_20260914/structure/main14/pruned'
REFILTERED_PREFIX = 'workspace/wgs/dt_research_p2_20260924/population/analysis_full_v1/pca/refiltered/baseline/pruned'
VCF = 'workspace/wgs/population_audit_20260914/inputs/Dtenuipes_population.vcf.gz'
FAI = 'workspace/wgs/genome.fa.fai'
HEADER_SHA = 'f9a99bc5d89c3b69463fd459f015dafa94d3d54284a592e166552651998566f9'
FAI_SHA = '2c70d093718e9dbe8620a6af948a7c896cf5102b7c8e21717c79ba1ee50d9adb'
HISTORICAL_VCF_SHA = '43472edd28b406493d5745d71a1d7d403d2086d810e59d4f50194858a176adb9'
EXPECTED = {'original': {'bed': '71da3d9c5081d4100ba1718126ee6a4396c1797fc725cc36eacc28349ab410a0', 'bim': '09b3cefd6cec30746ab7d4f856749d5724ff46e13d515d99d45cd6765dcd9b64', 'fam': 'c2a43d3f6821b3195277fbae34250c5aee93e5dd5ab9f09b7b12cf90be79e5ac'},
            'refiltered': {'bed': '71da3d9c5081d4100ba1718126ee6a4396c1797fc725cc36eacc28349ab410a0', 'bim': '229f212c7134ca383366b43838c372552b2a511cf52b530fd0d41c4c771a642f', 'fam': 'c2a43d3f6821b3195277fbae34250c5aee93e5dd5ab9f09b7b12cf90be79e5ac'}}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def parse_bim(text):
    rows, seen = [], set()
    for line in text.splitlines():
        f = line.split()
        if len(f) != 6 or not f[3].isdigit() or int(f[3]) < 1:
            raise ValueError('Malformed BIM')
        key = f[0], int(f[3])
        if key in seen:
            raise ValueError('Duplicate BIM coordinate, ambiguous VCF match')
        if f[4] not in 'ACGT' or f[5] not in 'ACGT' or len(f[4]) != 1 or len(f[5]) != 1 or f[4] == f[5]:
            raise ValueError('BIM is not a distinct biallelic ACGT SNP')
        seen.add(key)
        rows.append(f)
    if not rows:
        raise ValueError('Empty BIM')
    return rows


def bind_bims(original_text, refiltered_text):
    old, new = parse_bim(original_text), parse_bim(refiltered_text)
    if len(old) != len(new):
        raise ValueError('BIM counts differ')
    markers = []
    for i, (a, b) in enumerate(zip(old, new)):
        # IDs may differ; all other fields including allele order must agree.
        if a[:1] + a[2:] != b[:1] + b[2:]:
            raise ValueError('Original/refiltered BIM row or allele order differs at ' + str(i))
        parts = b[1].rsplit(':', 3)
        if len(parts) != 4 or parts[:2] != [b[0], b[3]] or set(parts[2:]) != set(b[4:]):
            raise ValueError('Refiltered marker ID coordinate/allele-label set differs from BIM')
        chrom, pos, label1, label2 = parts
        if label1 == label2 or len(label1) != 1 or len(label2) != 1:
            raise ValueError('Invalid historical ID allele labels')
        # PLINK 1.9 $1/$2 names labels in ASCII-sort order, not REF/ALT.
        # Real REF/ALT is bound only when the source VCF record is consumed.
        markers.append({'index': i, 'CHROM': chrom, 'POS': int(pos),
                        'historical_ID_alleles': label1 + ',' + label2,
                        'source_REF': None, 'source_ALT': None,
                        'original_ID': a[1], 'refiltered_ID': b[1], 'A1': a[4], 'A2': a[5]})
    return markers


def parse_fam(text):
    rows = [line.split() for line in text.splitlines()]
    if any(len(x) != 6 for x in rows):
        raise ValueError('Malformed FAM')
    names = [x[1] for x in rows]
    if names != SAMPLES:
        raise ValueError('Original 14 sample identity/order differs')
    return names


def gt_and_quality(raw):
    """Original GT -> DP5-100/GQ20 mask, with first-failure reason.

    Partial/missing GT is missing, never imputed reference. Unsupported complete
    encodings and malformed quality values fail rather than silently recode.
    """
    f = raw.split(':')
    if len(f) != 5:
        raise ValueError('Expected GT:AD:DP:GQ:PL')
    gt, ad, dp, gq, pl = f
    alleles = re.split(r'[/|]', gt)
    if any(a in ('', '.') for a in alleles):
        return None, None, 'GT_missing_or_partial'
    if len(alleles) != 2 or any(a not in ('0', '1') for a in alleles):
        raise ValueError('Unsupported complete GT encoding')
    dosage = sum(map(int, alleles))
    if dp == '.':
        return dosage, None, 'DP_missing'
    if not re.fullmatch(r'\d+', dp):
        raise ValueError('Invalid DP')
    if not 5 <= int(dp) <= 100:
        return dosage, None, 'DP_outside_5_100'
    if gq == '.':
        return dosage, None, 'GQ_missing'
    if not re.fullmatch(r'\d+', gq):
        raise ValueError('Invalid GQ')
    if int(gq) < 20:
        return dosage, None, 'GQ_below_20'
    return dosage, dosage, 'baseline_eligible'


def pl_details(raw):
    if raw == '.':
        return None, 'whole_missing'
    fields = raw.split(',')
    if len(fields) != 3:
        raise ValueError('PL wrong cardinality for diploid biallelic SNP')
    if any(x == '.' for x in fields):
        return None, 'partial_missing'
    if any(not re.fullmatch(r'\d+', x) for x in fields):
        raise ValueError('PL not three nonnegative integers')
    values = tuple(map(int, fields))
    if min(values) != 0:
        return values, 'minimum_nonzero'
    return values, 'normalized_integer_triad'


def dosage_to_bed(dosage, a1, a2, ref, alt):
    if {a1, a2} != {ref, alt} or a1 == a2:
        raise ValueError('BIM/VCF allele identity differs; no strand rescue')
    if dosage is None:
        return 1
    if dosage == 1:
        return 2
    if dosage not in (0, 2):
        raise ValueError('Invalid ALT dosage')
    allele = ref if dosage == 0 else alt
    return 0 if allele == a1 else 3


def bed_codes(packed, n=14):
    if len(packed) != (n + 3) // 4:
        raise ValueError('BED row byte count differs')
    codes = [(packed[i // 4] >> (2 * (i % 4))) & 3 for i in range(n)]
    if n % 4 and packed[-1] >> (2 * (n % 4)):
        raise ValueError('Nonzero unused BED padding bits')
    return codes
