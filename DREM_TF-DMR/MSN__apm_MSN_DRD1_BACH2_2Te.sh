#!/usr/bin/env bash
set -euo pipefail

# --- INPUTS ---
FILE=/u/project/jflint/heffel/BICAN3/DMR/MSN_DRD1-EPHA4_age/bican_2025_2T_eMSN_paths.txt.tsv_sort_strandmerged_CpG.tsv.gz
DMRS=~/epigenome/bican/protein_coding/DRD1-EPHA4_dms2_hypo_merged.bed   # chr start end id
OUTDIR=~/epigenome/bican/protein_coding
mkdir -p "$OUTDIR"

LABEL="2Te_EPHA4"
CPGS="$OUTDIR/cpgs.${LABEL}.bed"
OUT="$OUTDIR/DRD1-EPHA4_dms2_hypo_merged.${LABEL}.tsv"

# Sorting controls (tune if needed)
TMPDIR_SORT="${TMPDIR:-$OUTDIR/tmp_sort}"
SORT_MEM="${SORT_MEM:-2G}"
SORT_THREADS="${SORT_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
mkdir -p "$TMPDIR_SORT"

# Ensure DMRs are sorted in a bedtools-compatible way
DMRS_SORTED="$OUTDIR/DRD1-EPHA4_dms2_hypo_merged.sorted.bed"
if [[ ! -s "$DMRS_SORTED" ]]; then
  echo "[*] Sorting DMRs -> $DMRS_SORTED" >&2
  bedtools sort -i "$DMRS" > "$DMRS_SORTED"
fi

echo "[*] Processing $LABEL" >&2

# 1) Create CpG BED if it doesn't already exist (or is empty)
if [[ ! -s "$CPGS" ]]; then
  echo "    [+] Building CpGs: $CPGS" >&2
  zcat -f "$FILE" \
    | awk 'BEGIN{OFS="\t"} $4 ~ /^CG/ {print $1,$2-1,$2,$5,$6}' \
    | LC_ALL=C sort -T "$TMPDIR_SORT" -S "$SORT_MEM" --parallel="$SORT_THREADS" -k1,1 -k2,2n \
    > "$CPGS"
else
  echo "    [=] CpGs exists, skipping: $CPGS" >&2
fi

# 2) Aggregate CpGs over DMRs (DMRs col4 is the ID to carry through)
echo "    [+] Aggregating over DMRs -> $OUT" >&2
bedtools map -sorted \
  -a "$DMRS_SORTED" \
  -b "$CPGS" \
  -c 4,5,1 \
  -o sum,sum,count \
  -null 0 \
| awk -v OFS="\t" '{
    id=$4;
    m=$(NF-2); cov=$(NF-1); n=$NF;
    beta = (cov==0 || n==0) ? "NA" : (m/cov);
    print id, beta, m, cov, n
  }' > "$OUT"

echo "[✓] Done: $CPGS and $OUT" >&2
