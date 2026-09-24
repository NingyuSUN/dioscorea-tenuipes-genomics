#!/bin/bash
set -euo pipefail

ASSEMBLY="${1:-final_assembly.fa}"
LINEAGE="${2:-embryophyta_odb10}"
OUTPUT_DIR="${3:-busco_out}"
THREADS="${4:-24}"

mkdir -p "${OUTPUT_DIR}"

busco \
    -i "${ASSEMBLY}" \
    -l "${LINEAGE}" \
    -m genome \
    -c "${THREADS}" \
    -o busco_run \
    --out_path "${OUTPUT_DIR}"
