#!/usr/bin/env bash
set -euo pipefail

# --- INPUTS ---
MANIFEST=~/epigenome/bican/protein_coding/MSN_DRD1_BACH2_order.txt
DMRS=~/epigenome/bican/protein_coding/DRD1-BACH2_dms2_hypo_merged.bed   # chr start end id
OUTDIR=~/epigenome/bican/protein_coding
mkdir -p "$OUTDIR"

# Sorting controls (tune if needed)
TMPDIR_SORT="${TMPDIR:-$OUTDIR/tmp_sort}"
SORT_MEM="${SORT_MEM:-2G}"
SORT_THREADS="${SORT_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
mkdir -p "$TMPDIR_SORT"

# Ensure DMRs are sorted in a bedtools-compatible way
DMRS_SORTED="$OUTDIR/DRD1-BACH2_dms2_hypo_merged.sorted.bed"
bedtools sort -i "$DMRS" > "$DMRS_SORTED"

# --- PROCESS EACH file ---
while read -r FILE; do
  [[ -z "${FILE:-}" || "${FILE:0:1}" == "#" ]] && continue

  if [[ ! -f "$FILE" ]]; then
    echo "WARNING: missing file: $FILE" >&2
    continue
  fi

  # Derive LABEL from filename (handles bican_2025_ and bican_2025_merge_)
  base=$(basename "$FILE")
  tmp=${base#bican_2025_merge_}
  tmp=${tmp#bican_2025_}
  LABEL=${tmp%%.tsv*}

  CPGS="$OUTDIR/cpgs.${LABEL}.bed"
  OUT="$OUTDIR/DRD1-BACH2_dms2_hypo_merged.${LABEL}.tsv"

  echo "[*] Processing $LABEL" >&2

  # 0) If DMR TSV already exists, skip everything for this label
  if [[ -s "$OUT" ]]; then
    echo "    [=] DMR TSV exists, skipping: $OUT" >&2
    continue
  fi

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

done < "$MANIFEST"

echo "[✓] Per-timepoint DMR summaries written to $OUTDIR (DRD1-BACH2_dms2_hypo_merged.*.tsv)" >&2
