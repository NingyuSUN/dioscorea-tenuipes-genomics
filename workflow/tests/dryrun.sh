#!/usr/bin/env bash
# Check that the Snakemake DAG resolves, with empty placeholder inputs.
# Nothing is executed (-n). Usage: bash workflow/tests/dryrun.sh [snakemake]
set -euo pipefail
SMK=${1:-snakemake}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

mkdir -p "$T/in/audit" "$T/in/fast5" "$T/in/RNA_Seq"
touch "$T/in/pop.vcf.gz" "$T/in/manifest.tsv" "$T/in/r1.fq.gz" "$T/in/r2.fq.gz" \
      "$T/in/adapters.fa" "$T/in/p1.faa" "$T/in/p2.faa"
cat > "$T/dryrun.yaml" <<EOF
results_dir: $T/results
frozen_dirs: [$T/frozen]
inputs:
  population_vcf: $T/in/pop.vcf.gz
  sample_manifest: $T/in/manifest.tsv
  audit_package: $T/in/audit
assembly:
  outdir: $T/assembly
  fast5_dir: $T/in/fast5
  illumina_r1: $T/in/r1.fq.gz
  illumina_r2: $T/in/r2.fq.gz
  adapters: $T/in/adapters.fa
  rnaseq_dir: $T/in/RNA_Seq
  proteins: [$T/in/p1.faa, $T/in/p2.faa]
EOF

cd "$REPO"
# --configfile takes several values, so it goes last.
run() { "$SMK" "$@" --cores 4 --configfile config/config.yaml config/assembly.yaml "$T/dryrun.yaml"; }
run all -n --quiet rules
run assembly_all -n --quiet rules
# Lint is advisory: accepted findings are listed in workflow/README.md. Fail
# only on rules whose shell reads globals, which hides them from provenance.
run --lint > "$T/lint.txt" 2>&1 || true
if grep -q "Shell command directly uses variable" "$T/lint.txt"; then
  grep -B3 "Shell command directly uses variable" "$T/lint.txt" >&2; exit 1
fi

# The frozen-directory guard must refuse a directory listed in frozen_dirs.
if "$SMK" -n --cores 1 --config="results_dir=$T/frozen/x" --configfile config/config.yaml "$T/dryrun.yaml" >"$T/guard.log" 2>&1; then
  echo "FAIL: frozen-directory guard did not stop the run" >&2; exit 1
fi
grep -q "completed evidence is read-only" "$T/guard.log"
echo "Dry run, lint and frozen-directory guard: OK"
