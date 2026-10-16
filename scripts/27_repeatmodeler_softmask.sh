#!/bin/bash
# TEMPLATE: review RepeatMasker -pa batch count versus actual CPU allocation before execution; no completed repeat analysis is claimed here.
# De novo repeat library + soft-masking of the Kita1-scaffolded assembly.
# The soft-masked output is the required input for re-running BRAKER: the
# previous annotation was run with --softmasking on a genome that contained
# zero lowercase bases, i.e. with no repeat masking at all.
# Usage: 27_repeatmodeler_softmask.sh <scaffold.fasta> [work_dir] [threads]
# Requires RepeatModeler and RepeatMasker on PATH.
set -e
GENOME=$(realpath "${1:?Usage: 27_repeatmodeler_softmask.sh <scaffold.fasta> [work_dir] [threads]}")
WORK=${2:-repeat}
THREADS=${3:-16}

mkdir -p $WORK
WORK=$(realpath "$WORK")
cd $WORK

echo "### start: $(date)"
echo "### genome: $GENOME"

echo "=== [1/3] BuildDatabase ==="
BuildDatabase -name Dtenuipes "$GENOME"

echo "=== [2/3] RepeatModeler (de novo library, with LTR structural search) ==="
RepeatModeler -database Dtenuipes -threads $THREADS -LTRStruct

echo "=== [3/3] RepeatMasker (soft-masking with -xsmall) ==="
LIB=$(ls -t $WORK/Dtenuipes-families.fa 2>/dev/null | head -1)
if [ -z "$LIB" ]; then
  echo "ERROR: RepeatModeler library not found" >&2
  exit 1
fi
echo "using library: $LIB ($(grep -c '^>' $LIB) families)"

mkdir -p $WORK/masked
RepeatMasker -lib "$LIB" -xsmall -gff -pa $THREADS -dir $WORK/masked "$GENOME"

echo "### done: $(date)"

echo "=== masking summary ==="
cat $WORK/masked/*.tbl 2>/dev/null | head -40

echo "=== verify soft-masking actually applied ==="
MASKED=$WORK/masked/$(basename $GENOME).masked
if [ -f "$MASKED" ]; then
  grep -v '^>' "$MASKED" | awk '
  { n=length($0); for(i=1;i<=n;i++){ c=substr($0,i,1)
      if (c ~ /[acgt]/) lo++; else if (c ~ /[ACGT]/) up++; else ot++ } }
  END { t=lo+up+ot; printf "lowercase (masked): %d (%.2f%%)\nuppercase: %d (%.2f%%)\nN/other: %d\n", lo,100*lo/t,up,100*up/t,ot }'
else
  echo "masked fasta not found at $MASKED"
  ls -la $WORK/masked/
fi
