#!/bin/bash
set -euo pipefail

# ============================================================
# 18_novoplasty_chloroplast.sh
# Run NOVOPlasty for chloroplast genome assembly.
#
# Usage:
#   bash 18_novoplasty_chloroplast.sh <config.txt> [novoplasty_script]
#
# Example:
#   bash 18_novoplasty_chloroplast.sh config/novoplasty_chloroplast_config.txt NOVOPlasty4.3.5.pl
#
# Notes:
# - Edit the config file paths before running
# - The seed sequence can be the reference-guided cp consensus
# ============================================================

CONFIG_FILE="${1:-config/novoplasty_chloroplast_config.txt}"
NOVOPLASTY_BIN="${2:-NOVOPlasty4.3.5.pl}"

[[ -f "${CONFIG_FILE}" ]] || { echo "Error: config file not found: ${CONFIG_FILE}"; exit 1; }
[[ -f "${NOVOPLASTY_BIN}" ]] || command -v "${NOVOPLASTY_BIN}" >/dev/null 2>&1 || { echo "Error: NOVOPlasty not found."; exit 1; }

perl "${NOVOPLASTY_BIN}" -c "${CONFIG_FILE}"

echo "NOVOPlasty chloroplast assembly completed."
