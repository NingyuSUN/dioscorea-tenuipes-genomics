#!/usr/bin/env python3
"""
HISTORICAL, UNSUITABLE FOR ANNOTATION: reconstruct hints from incomplete records.

Use 29_prepare_rna_evidence.py with complete transcript structures instead.
The old reconstruction below is retained for inspection, not endorsed as a fix.

The original file writes every record with feature type "hint", which BRAKER
does not recognise, so every hint was discarded with a neutral bonus/malus.
The records are exons (minimum length 4 bp rules out introns), so the introns
are the gaps between consecutive exons of the same transcript.

Introns are the valuable hint type: they pin down splice sites, which is what
BRAKER primarily uses RNA-seq for. This script derives them and validates the
call against canonical GT-AG splice sites before writing anything.
"""
import sys, re, random
from collections import defaultdict

# Usage: 24_repair_stringtie_hints.py <hintsfile.gff> <genome.fa> [out_dir]
if len(sys.argv) < 3:
    sys.exit("Usage: 24_repair_stringtie_hints.py <hintsfile.gff> <genome.fa> [out_dir]")
HINTS, GENOME = sys.argv[1], sys.argv[2]
OUT_DIR = sys.argv[3] if len(sys.argv) > 3 else "."
OUT_INTRON = f"{OUT_DIR}/hints_fixed_introns.gff"
OUT_EXONPART = f"{OUT_DIR}/hints_fixed_exonparts.gff"

# ---------- 1. parse ----------
groups = defaultdict(list)   # (contig, strand, grp) -> [(start, end, mult)]
n_rec = 0
with open(HINTS) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        contig, start, end, strand, attr = f[0], int(f[3]), int(f[4]), f[6], f[8]
        m = re.search(r"grp=Parent=([^;]+)", attr)
        if not m:
            continue
        mult = 1
        mm = re.search(r"mult=(\d+)", attr)
        if mm:
            mult = int(mm.group(1))
        groups[(contig, strand, m.group(1))].append((start, end, mult))
        n_rec += 1

print(f"parsed {n_rec} exon records in {len(groups)} transcript groups", flush=True)

# ---------- 2. derive introns ----------
introns = []   # (contig, start, end, strand, mult)
for (contig, strand, grp), exons in groups.items():
    exons.sort()
    for i in range(len(exons) - 1):
        istart = exons[i][1] + 1
        iend   = exons[i + 1][0] - 1
        if iend < istart:          # overlapping/abutting exons -> no intron
            continue
        mult = min(exons[i][2], exons[i + 1][2])
        introns.append((contig, istart, iend, strand, mult))

print(f"derived {len(introns)} candidate introns", flush=True)
if not introns:
    sys.exit("no introns derived - aborting")

lens = sorted(e - s + 1 for _, s, e, _, _ in introns)
n = len(lens)
print(f"intron length: min={lens[0]} p25={lens[n//4]} median={lens[n//2]} "
      f"p75={lens[3*n//4]} p95={lens[int(n*0.95)]} max={lens[-1]}", flush=True)

# ---------- 3. validate with canonical splice sites ----------
need = defaultdict(list)
sample = random.Random(0).sample(introns, min(5000, len(introns)))
for contig, s, e, strand, _ in sample:
    need[contig].append((s, e, strand))

seqs = {}
want = set(need)
cur = None
buf = []
with open(GENOME) as fh:
    for line in fh:
        if line.startswith(">"):
            if cur in want:
                seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]
            buf = []
        elif cur in want:
            buf.append(line.strip())
    if cur in want:
        seqs[cur] = "".join(buf)

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
canonical = other = missing = 0
motifs = defaultdict(int)
for contig, lst in need.items():
    seq = seqs.get(contig)
    if seq is None:
        missing += len(lst)
        continue
    for s, e, strand in lst:
        donor = seq[s - 1: s + 1].upper()
        acceptor = seq[e - 2: e].upper()
        if strand == "-":
            donor, acceptor = (acceptor.translate(COMP)[::-1],
                               donor.translate(COMP)[::-1])
        motif = donor + "-" + acceptor
        motifs[motif] += 1
        if motif in ("GT-AG", "GC-AG", "AT-AC"):
            canonical += 1
        else:
            other += 1

tot = canonical + other
print(f"\nsplice-site check on {tot} sampled introns ({missing} skipped):")
print(f"  canonical (GT-AG / GC-AG / AT-AC): {canonical} ({100*canonical/tot:.1f}%)")
for motif, c in sorted(motifs.items(), key=lambda x: -x[1])[:6]:
    print(f"    {motif}: {c} ({100*c/tot:.1f}%)")

if canonical / tot < 0.80:
    sys.exit("\nABORT: canonical splice-site rate too low - the records are "
             "probably not exons after all. Nothing written.")

# ---------- 4. write ----------
with open(OUT_INTRON, "w") as out:
    for contig, s, e, strand, mult in sorted(introns):
        out.write(f"{contig}\tStringTie\tintron\t{s}\t{e}\t.\t{strand}\t.\t"
                  f"mult={mult};src=E\n")

with open(OUT_EXONPART, "w") as out:
    for (contig, strand, grp), exons in sorted(groups.items()):
        for s, e, mult in exons:
            out.write(f"{contig}\tStringTie\texonpart\t{s}\t{e}\t.\t{strand}\t.\t"
                      f"mult={mult};src=E\n")

print(f"\nwrote {OUT_INTRON}")
print(f"wrote {OUT_EXONPART}")
