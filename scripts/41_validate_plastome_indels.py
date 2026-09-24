#!/usr/bin/env python3
"""Read-level, alignment-independent indel assay using diagnostic junction k-mers.

Both hypotheses retain the candidate's nearby SNPs/other indels; only the tested
gap is reversed. Counts are read-name groups (mates deduplicated), not independent
PCR molecules. Exact BQ20 31-mers with >=8 bases on either side of a breakpoint.
Long indels need both breakpoints; evidence does not establish long-range phase.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def canonical(s):
    return min(s, rc(s))


def fasta(path):
    return ''.join(x.strip() for x in Path(path).read_text().splitlines() if not x.startswith('>')).upper()


def circular(seq, start, end):
    return ''.join(seq[i % len(seq)] for i in range(start, end))


def markers(seq, boundary, k=31, anchor=8):
    return {canonical(seq[i:i+k]) for i in range(max(0,boundary-k+anchor), min(boundary-anchor, len(seq)-k)+1)
            if len(seq[i:i+k]) == k and set(seq[i:i+k]) <= set('ACGT')}


def kmers(seq, k):
    return {canonical(seq[i:i+k]) for i in range(len(seq)-k+1)}


def shift_range(genome, s, e, allele):
    """All equivalent placements of replacement within circular local sequence.

    Test by equality of the entire mutated 1 kb context, so reported ambiguous
    positions reflect actual sequence equivalence, not just homopolymer length.
    Range is capped at 200 bases each direction and cap status is explicit.
    """
    left = s-500
    context = circular(genome, left, e+500)
    target = context[:500]+allele+context[500+e-s:]
    hits = []
    d = e-s
    for delta in range(-200,201):
        q = 500+delta
        # An insertion/deletion may rotate its motif under an equivalent shift.
        if d:
            trial = context[:q]+context[q+d:]
        else:
            inserted = target[q:q+len(allele)]
            trial = context[:q]+inserted+context[q:]
        if trial == target:
            hits.append(delta)
    assert 0 in hits
    return min(hits), max(hits)


def build_assays(genome, events, k):
    counts = collections.Counter(canonical(circular(genome,i,i+k)) for i in range(len(genome)))
    lookup = collections.defaultdict(set)
    for index, event in enumerate(events):
        s,e = int(event['own_start0']), int(event['own_end0'])
        allele = event['published_allele']
        assert genome[s:e] == event['own_allele']
        flank = 250
        left,right = circular(genome,s-flank,s), circular(genome,e,e+flank)
        own, other = left+genome[s:e]+right, left+allele+right
        ownall, otherall = kmers(own,k), kmers(other,k)
        assays = {}
        for side,boundary in [('left',flank),('right',flank+e-s)]:
            assays['own_'+side] = markers(own,boundary,k)-otherall
        for side,boundary in [('left',flank),('right',flank+len(allele))]:
            assays['public_'+side] = {x for x in markers(other,boundary,k)-ownall if counts[x] == 0}
        for label,keys in assays.items():
            event[label+'_diagnostic_kmers'] = len(keys)
            for key in keys:
                lookup[key].add((index,label))
                lookup[rc(key)].add((index,label))
        ownkeys = assays['own_left'] | assays['own_right']
        event['own_marker_max_genome_copies'] = max([counts[x] for x in ownkeys] or [0])
        event['own_marker_min_genome_copies'] = min([counts[x] for x in ownkeys] or [0])
        lo,hi = shift_range(genome,s,e,allele)
        event['equivalent_gap_shift_min_bp'] = lo
        event['equivalent_gap_shift_max_bp'] = hi
        event['gap_shift_search_capped'] = lo == -200 or hi == 200
    return lookup


def scan_fastq(paths, lookup, events, k):
    supports = [collections.defaultdict(set) for _ in events]
    reads = 0
    evidence_reads = 0
    manifests = []
    for path in paths:
        digest = hashlib.sha256()
        with open(path,'rb') as raw:
            for block in iter(lambda:raw.read(1024*1024), b''):
                digest.update(block)
        count = 0
        with gzip.open(path, 'rt') as fh:
            while True:
                name = fh.readline().strip()
                if not name:
                    break
                seq,plus,quality = fh.readline().strip().upper(),fh.readline().strip(),fh.readline().strip()
                assert name.startswith('@') and plus.startswith('+') and len(seq)==len(quality)
                read_id = name[1:].split()[0]
                if read_id.endswith('/1') or read_id.endswith('/2'):
                    read_id = read_id[:-2]
                hits = set()
                for i in range(len(seq)-k+1):
                    values = lookup.get(seq[i:i+k])
                    if values and min(quality[i:i+k]) >= '5':  # Phred+33, Q20
                        hits.update(values)
                for event,label in hits:
                    supports[event][label].add(read_id)
                evidence_reads += bool(hits)
                reads += 1
                count += 1
        manifests.append(dict(file=str(path), reads=count, sha256=digest.hexdigest()))
        print('Scanned %s: %d reads' % (path.name,count), flush=True)
    return supports, dict(total_reads=reads, reads_with_diagnostic_marker=evidence_reads, input_files=manifests)


def classify(event, support):
    own = support['own_left'] | support['own_right']
    other = support['public_left'] | support['public_right']
    mixed = own & other
    for label in ('own_left','own_right','public_left','public_right'):
        event[label+'_read_names'] = len(support[label]-mixed)
    event['own_total_read_names'] = len(own-mixed)
    event['public_total_read_names'] = len(other-mixed)
    event['conflicting_read_names'] = len(mixed)
    event['own_both_junctions_read_names'] = len((support['own_left'] & support['own_right'])-mixed)
    nown,nother = len(own-mixed),len(other-mixed)
    fraction = nown/(nown+nother) if nown+nother else None
    event['own_fraction_among_informative_names'] = fraction if fraction is not None else ''
    enough = all(event[x+'_diagnostic_kmers'] > 0 for x in ('own_left','own_right','public_left','public_right'))
    if not enough:
        status = 'unresolved_no_discriminating_junction'
    elif nother >= 3 and fraction is not None and fraction < 0.90:
        status = 'conflicting_or_public_allele_evidence'
    elif min(event['own_left_read_names'], event['own_right_read_names']) >= 5 and fraction is not None and fraction >= 0.95 and len(mixed) <= 0.05*max(1,len(own|other)):
        status = 'candidate_supported'
    else:
        status = 'unresolved_low_or_mixed_support'
    event['validation_status'] = status
    event['repeat_copy_resolution'] = ('no_marker' if event['own_marker_min_genome_copies'] == 0 else ('copy_unresolved' if event['own_marker_min_genome_copies'] > 1 else 'has_copy_specific_marker'))
    return status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    genome = fasta(args.out/'Dtenuipes.path2.fasta')
    events = list(csv.DictReader(open(args.out/'indel_events.annotated.tsv'), delimiter='\t'))
    lookup = build_assays(genome, events, 31)
    supports, meta = scan_fastq(sorted((args.root/'inputs').glob('extended_*fq.gz')), lookup, events, 31)
    assert meta['total_reads'] == 578379
    counts = collections.Counter(classify(e,s) for e,s in zip(events,supports))
    with open(args.out/'indel_validation.tsv','w') as fh:
        w=csv.DictWriter(fh,fieldnames=list(events[0]),delimiter='\t')
        w.writeheader()
        w.writerows(events)
    meta.update(dict(status_counts=dict(counts), event_count=len(events),
                     kmer_size=31, minimum_quality=20, minimum_anchor=8,
                     candidate_support_criteria='Both own junctions >=5 distinct read names, own fraction >=0.95; both hypotheses have markers; mixed <=5% of evidence names.',
                     method='Exact junction markers against both local hypotheses, scanned across all four existing recruited FASTQs without BAM alignment filtering. Each counterfactual changes only the tested indel.',
                     limitations=['Recruited assembly reads, not independent full WGS.', 'Read names merge overlapping mates but do not remove PCR duplicates.',
                                  'Marker counts/fractions are not unbiased allele-frequency estimates.', 'Two IR copies may share evidence.',
                                  'Long allele breakpoint support does not establish complete allele phase.', 'Gap placement can be nonunique in repeats.',
                                  'Reference-favoring evidence concerns a counterfactual local allele, not independent resequencing of PV612491.1.']))
    (args.out/'indel_validation_summary.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta,indent=2))


if __name__ == '__main__':
    main()
