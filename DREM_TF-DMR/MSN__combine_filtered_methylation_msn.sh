#!/usr/bin/env bash
set -euo pipefail

DIR=~/epigenome/bican/protein_coding
OUT=~/epigenome/bican/protein_coding/MSN_DRD1_EPHA4_dmr_beta.matrix.tsv

# Column display names in order
ORDER=(
  2Te
  2T
  3T
  1m
  4-7m
  adult
)

# Map display names to file suffixes (the part after "DRD1-EPHA4_dms2_hypo_merged.")
declare -A SUFFIX
SUFFIX[2Te]="2Te_EPHA4"
SUFFIX[2T]="2T_EPHA4_paths.txt"
SUFFIX[3T]="3T_paths.txt"
SUFFIX[1m]="1m_paths.txt"
SUFFIX[4-7m]="4-7mpaths.txt"
SUFFIX[adult]="adult_paths.txt"

FIRST="${ORDER[0]}"

# Input per-cluster files produced by your aggregation script:
#   $DIR/DRD1-EPHA4_dms2_hypo_merged.<SUFFIX>.tsv
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
PASTE_ARGS+=("<(cut -f1 \"$DIR/DRD1-EPHA4_dms2_hypo_merged.${SUFFIX[$FIRST]}.tsv\")")
for L in "${ORDER[@]}"; do
  PASTE_ARGS+=("<(cut -f2 \"$DIR/DRD1-EPHA4_dms2_hypo_merged.${SUFFIX[$L]}.tsv\")")
done

# shellcheck disable=SC2206
eval paste ${PASTE_ARGS[@]} >> "$OUT"

echo "Wrote $OUT"
