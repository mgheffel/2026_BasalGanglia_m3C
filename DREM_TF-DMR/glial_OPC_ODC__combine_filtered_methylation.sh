#!/usr/bin/env bash
set -euo pipefail

DIR=~/epigenome/bican/protein_coding
OUT=~/epigenome/bican/protein_coding/OPC_ODC_dmr_beta.matrix.tsv

# OPC_ODC lineage cluster order (edit to match your actual labels/files)
ORDER=(
  2T_GPC
  GPC
  OPC
  tOPC
  tODC
  ODC
  adult_ODC
)

FIRST="${ORDER[0]}"

# Input per-cluster files produced by your aggregation script:
#   $DIR/dmrs_OPDC_hypo.<LABEL>.tsv
# Expected columns:
#   1 dmr_id
#   2 beta
#   3 mC_sum
#   4 C_total_sum
#   5 CpG_count

# Header
{
  printf "dmr"
  for L in "${ORDER[@]}"; do printf "\t%s" "$L"; done
  printf "\n"
} > "$OUT"

# Build paste args:
# first column = DMR IDs from FIRST,
# then beta column (col 2) from each cluster in ORDER.
declare -a PASTE_ARGS
PASTE_ARGS+=("<(cut -f1 \"$DIR/dmrs_OPDC_hypo.$FIRST.clean.tsv\")")
for L in "${ORDER[@]}"; do
  PASTE_ARGS+=("<(cut -f2 \"$DIR/dmrs_OPDC_hypo.$L.clean.tsv\")")
done

# shellcheck disable=SC2206
eval paste ${PASTE_ARGS[@]} >> "$OUT"

echo "Wrote $OUT"
