#!/bin/bash
set -euo pipefail

# ============================================================
# 19_novoplasty_mitochondria.sh
# Run NOVOPlasty for mitochondrial genome assembly.
#
# Usage:
#   bash 19_novoplasty_mitochondria.sh <config.txt> [novoplasty_script]
#
# Example:
#   bash 19_novoplasty_mitochondria.sh config/novoplasty_mito_config.txt NOVOPlasty4.3.5.pl
#
# Notes:
# - This is a clean template because the original mito script
#   was not available.
# - Provide a mitochondrial seed sequence from a related species
#   or a conserved mitochondrial gene region.
# ============================================================

CONFIG_FILE="${1:-config/novoplasty_mito_config.txt}"
NOVOPLASTY_BIN="${2:-NOVOPlasty4.3.5.pl}"

[[ -f "${CONFIG_FILE}" ]] || { echo "Error: config file not found: ${CONFIG_FILE}"; exit 1; }
[[ -f "${NOVOPLASTY_BIN}" ]] || command -v "${NOVOPLASTY_BIN}" >/dev/null 2>&1 || { echo "Error: NOVOPlasty not found."; exit 1; }

perl "${NOVOPLASTY_BIN}" -c "${CONFIG_FILE}"

echo "NOVOPlasty mitochondrial assembly completed."
