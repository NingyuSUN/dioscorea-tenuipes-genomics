#!/usr/bin/env bash
set -euo pipefail
ROOT=${1:?Usage: 36_validate_organelle_evidence.sh ROOT GENOME_FASTA}
GENOME=$(realpath "${2:?Usage: 36_validate_organelle_evidence.sh ROOT GENOME_FASTA}")
cd "$ROOT"
MM=${MINIMAP2:-minimap2}
SAM=${SAMTOOLS:-samtools}
PY=${PYTHON:-python3}
mkdir validation
"$PY" 35_recover_organelle_sequences.py "$ROOT" > validation/recovery.stdout.json
"$MM" --version > validation/minimap2.version.txt
"$SAM" --version > validation/samtools.version.txt
cat inputs/D.tokoro_chro.fasta inputs/D.tokoro_mt.fasta > validation/comparative_refs.fasta
"$MM" -t 8 -x asm20 --secondary=yes -N 100 validation/comparative_refs.fasta "$GENOME" > validation/genome_vs_organelle.paf 2> validation/genome_vs_organelle.log
"$MM" -t 4 -x asm20 --secondary=yes -N 100 validation/comparative_refs.fasta recovered/Dtenuipes_Izu_reference_candidate.plastome.path1.fasta > validation/plastome_vs_refs.paf 2> validation/plastome_vs_refs.log
for TAG in native rotated; do
  if [[ "$TAG" == native ]]; then REF=recovered/Dtenuipes_Izu_reference_candidate.plastome.path1.fasta; else REF=recovered/plastome.rotated_validation.fasta; fi
  "$SAM" faidx "$REF"
  "$MM" -t 8 -ax sr "$REF" inputs/extended_1_paired.fq.gz inputs/extended_2_paired.fq.gz 2> "validation/$TAG.paired.log" | "$SAM" sort -@ 2 -m 512M -o "validation/$TAG.paired.bam"
  "$MM" -t 8 -ax sr "$REF" inputs/extended_1_unpaired.fq.gz 2> "validation/$TAG.unpaired1.log" | "$SAM" sort -@ 2 -m 512M -o "validation/$TAG.unpaired1.bam"
  "$MM" -t 8 -ax sr "$REF" inputs/extended_2_unpaired.fq.gz 2> "validation/$TAG.unpaired2.log" | "$SAM" sort -@ 2 -m 512M -o "validation/$TAG.unpaired2.bam"
  "$SAM" merge -@ 2 "validation/$TAG.all.bam" "validation/$TAG.paired.bam" "validation/$TAG.unpaired1.bam" "validation/$TAG.unpaired2.bam"
  "$SAM" index "validation/$TAG.all.bam"
  "$SAM" quickcheck -v "validation/$TAG.all.bam"
  "$SAM" flagstat "validation/$TAG.all.bam" > "validation/$TAG.flagstat.txt"
  "$SAM" depth -aa -q 20 -Q 0 "validation/$TAG.all.bam" > "validation/$TAG.depth.tsv"
  "$SAM" depth -aa -q 20 -Q 20 "validation/$TAG.all.bam" > "validation/$TAG.unique_depth.tsv"
done
"$PY" - <<'PY'
import pysam
b=pysam.AlignmentFile('validation/rotated.all.bam','rb'); length=b.lengths[0]; junction=length-length//2
spans=[]
for r in b.fetch(b.references[0],junction-1,junction+1):
 if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality<20: continue
 blocks=r.get_blocks()
 if any(a<=junction-30 and z>=junction+30 for a,z in blocks):spans.append((r.query_name,r.mapping_quality))
with open('validation/circular_junction_support.tsv','w') as f:
 f.write('read_name\tmapq\n');f.writelines(f'{a}\t{q}\n' for a,q in spans)
print('junction_spanning_alignments',len(spans),'unique_fragments',len(set(a for a,q in spans)))
PY
touch validation/ALIGNMENT_COMPLETE
