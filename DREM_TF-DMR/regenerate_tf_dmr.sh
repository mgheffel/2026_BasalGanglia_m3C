#!/usr/bin/env bash
# Regenerate lineage-specific TF-DMR file using the chr22-corrected
# ChIP-Atlas top-3-per-TF experiment list, intersected with the union of
# (a) the paper's curated SRX whitelist and
# (b) 14 hand-curated new-only SRXs from the chr22-corrected list.
#
# Usage: bash regenerate_tf_dmr.sh <lineage>
#   lineage: mge | bach2 | epha4 | opc | ecge
#
# Defaults:
#   SUBSET, PAPER_LABELS, ADDENDUM_LABELS resolve relative to this script.
#   DMR_DIR defaults to $SCRIPT_DIR/dmrs (a subdirectory next to this script).
#
# Override any of these via environment, e.g.:
#   DMR_DIR=/path/to/my/dmrs bash regenerate_tf_dmr.sh mge
#   PAPER_LABELS=/elsewhere/paper.tsv bash regenerate_tf_dmr.sh mge
#
# DMR BED filenames the script expects under DMR_DIR:
#   MGE_P_all_dms2_hypo_merged.bed              (mge)
#   DRD1-BACH2_dms2_hypo_merged.bed             (bach2)
#   DRD1-EPHA4_dms2_hypo_merged.bed             (epha4)
#   OPDC_hypo_dmrs.200-250bp.id.bed             (opc)
#   CGE_LAMP5_allLAMP5_dms2_hypo_merged.bed     (ecge)

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <lineage>" >&2
  echo "       lineage: mge | bach2 | epha4 | opc | ecge" >&2
  exit 2
fi

LINEAGE="$1"

# Resolve script directory so default input paths are self-contained.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# DMR BED location (override with DMR_DIR=... in environment)
DMR_DIR="${DMR_DIR:-$SCRIPT_DIR/dmrs}"

case "$LINEAGE" in
  mge)
    DMRS="$DMR_DIR/MGE_P_all_dms2_hypo_merged.bed"
    OUT_NAME="TFcellACC_to_dmrs_mge_alltf_new.tsv"
    ;;
  bach2)
    DMRS="$DMR_DIR/DRD1-BACH2_dms2_hypo_merged.bed"
    OUT_NAME="TFcellACC_to_dmrs_bach2_alltf_new.tsv"
    ;;
  epha4)
    DMRS="$DMR_DIR/DRD1-EPHA4_dms2_hypo_merged.bed"
    OUT_NAME="TFcellACC_to_dmrs_epha4_alltf_new.tsv"
    ;;
  opc)
    DMRS="$DMR_DIR/OPDC_hypo_dmrs.200-250bp.id.bed"
    OUT_NAME="TFcellACC_to_dmrs_opc_alltf_new.tsv"
    ;;
  ecge)
    DMRS="$DMR_DIR/CGE_LAMP5_allLAMP5_dms2_hypo_merged.bed"
    OUT_NAME="TFcellACC_to_dmrs_ecge_alltf_new.tsv"
    ;;
  *)
    echo "ERROR: unknown lineage '$LINEAGE'" >&2
    echo "Valid: mge | bach2 | epha4 | opc | ecge" >&2
    exit 2
    ;;
esac

# --- Inputs (default to script directory; override via environment) ---
SUBSET="${SUBSET:-$SCRIPT_DIR/chip-atlas_top3_alltf_new.txt}"
PAPER_LABELS="${PAPER_LABELS:-$SCRIPT_DIR/paper_ACC_to_tf_cell.tsv}"
ADDENDUM_LABELS="${ADDENDUM_LABELS:-$SCRIPT_DIR/new_addendum_ACC_to_tf_cell.tsv}"

# --- Validate ---
[[ -f "$SUBSET" ]]          || { echo "ERROR: SUBSET missing: $SUBSET" >&2; exit 1; }
[[ -f "$DMRS" ]]            || { echo "ERROR: DMRS missing: $DMRS (override with DMR_DIR=...)" >&2; exit 1; }
[[ -f "$PAPER_LABELS" ]]    || { echo "ERROR: PAPER_LABELS missing: $PAPER_LABELS" >&2; exit 1; }
[[ -f "$ADDENDUM_LABELS" ]] || { echo "ERROR: ADDENDUM_LABELS missing: $ADDENDUM_LABELS" >&2; exit 1; }
command -v bedtools >/dev/null || { echo "ERROR: bedtools not on PATH (module load bedtools?)" >&2; exit 1; }

echo "[config] SCRIPT_DIR=$SCRIPT_DIR" >&2
echo "[config] DMR_DIR=$DMR_DIR" >&2
echo "[config] DMRS=$DMRS" >&2

# --- Work dir ---
WORKDIR="tf_cell_acc_dmrs_${LINEAGE}_alltf_new"
mkdir -p "$WORKDIR"

DMRS_SORTED="$WORKDIR/dmrs.sorted.bed"
bedtools sort -i "$DMRS" > "$DMRS_SORTED"

# 0) Build combined whitelist (paper + addendum), then intersect with new SRX list
COMBINED_LABELS="$WORKDIR/combined_ACC_to_tf_cell.tsv"
cat "$PAPER_LABELS" "$ADDENDUM_LABELS" | sort -u > "$COMBINED_LABELS"

ACC_LIST="$WORKDIR/acc_list.txt"
awk '
  NF==0 {next}
  $0 ~ /^#/ {next}
  {
    gsub(/^[ \t]+|[ \t]+$/, "", $0)
    n=split($0, a, "/")
    f=a[n]
    sub(/\.bed\.gz$/, "", f)
    print f
  }' "$SUBSET" | sort -u > "$WORKDIR/new_acc.txt"

cut -f1 "$COMBINED_LABELS" | sort -u > "$WORKDIR/whitelist_acc.txt"
comm -12 "$WORKDIR/new_acc.txt" "$WORKDIR/whitelist_acc.txt" > "$ACC_LIST"

n_new=$(wc -l < "$WORKDIR/new_acc.txt")
n_paper=$(wc -l < "$PAPER_LABELS")
n_addendum=$(wc -l < "$ADDENDUM_LABELS")
n_whitelist=$(wc -l < "$WORKDIR/whitelist_acc.txt")
n_use=$(wc -l < "$ACC_LIST")
echo "[0/3] SRX sets: new=$n_new, paper=$n_paper, addendum=$n_addendum, combined_whitelist=$n_whitelist, intersection=$n_use" >&2

# 1) ACC -> DMR overlaps (only for whitelisted SRXs)
echo "[1/3] Intersecting peaks with DMRs (lineage=$LINEAGE)..." >&2
> "$WORKDIR/acc_dmr_pairs.tsv"

while IFS= read -r bedpath; do
  [[ -z "${bedpath:-}" || "${bedpath:0:1}" == "#" ]] && continue
  acc=$(basename "$bedpath" .bed.gz)
  if ! grep -qx "$acc" "$ACC_LIST"; then
    continue
  fi
  [[ -f "$bedpath" ]] || { echo "WARNING: missing peak bed: $bedpath" >&2; continue; }
  echo "    [+] $acc" >&2

  zcat -f "$bedpath" \
    | bedtools intersect -a stdin -b "$DMRS_SORTED" -wa -wb \
    | awk -v ACC="$acc" 'BEGIN{OFS="\t"} {print ACC, $NF}' \
    >> "$WORKDIR/acc_dmr_pairs.tsv"

done < "$SUBSET"

sort -u -k1,1 -k2,2 "$WORKDIR/acc_dmr_pairs.tsv" -o "$WORKDIR/acc_dmr_pairs.tsv"
echo "    -> $(wc -l < "$WORKDIR/acc_dmr_pairs.tsv") unique ACC-DMR pairs" >&2

# 2) Use combined labels (filter to ACC_LIST)
echo "[2/3] Loading combined labels..." >&2
ACC_META="$WORKDIR/ACC_to_tf_cell.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==FNR{keep[$1]=1; next} ($1 in keep)' \
    "$ACC_LIST" "$COMBINED_LABELS" \
  | sort -u > "$ACC_META"

echo "    -> $(wc -l < "$ACC_META") ACCs with TF+Cell labels" >&2

# 3) Join and write outputs (with DREM header)
echo "[3/3] Joining and writing outputs..." >&2

OUT1="$WORKDIR/$OUT_NAME"
OUT3="$WORKDIR/TF__dmr_3col_${LINEAGE}_alltf_new.tsv"

printf 'TF\tGene\tInput\n' > "$OUT1"
printf 'TF\tGene\tInput\n' > "$OUT3"

join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 "$WORKDIR/acc_dmr_pairs.tsv") \
  <(sort -k1,1 "$ACC_META") \
| awk 'BEGIN{FS="\t"; OFS="\t"}
       {
         acc=$1; dmr=$2; tf=$3; cell=$4;
         key=tf "|" cell "|" acc
         print key, dmr, 1 >> "'"$OUT1"'"
         print tf, dmr, 1 >> "'"$OUT3"'"
       }'

# Summary
n_data=$(($(wc -l < "$OUT1") - 1))
n_tfs=$(awk -F'\t' 'NR>1{split($1,a,"|"); print a[1]}' "$OUT1" | sort -u | wc -l)
n_accs=$(awk -F'\t' 'NR>1{split($1,a,"|"); print a[3]}' "$OUT1" | sort -u | wc -l)
echo "Done."
echo "  $OUT1"
echo "    data rows:   $n_data"
echo "    unique TFs:  $n_tfs"
echo "    unique ACCs: $n_accs"