#!/bin/bash
#$ -N gene_contact_score
#$ -cwd
#$ -V
#$ -j y
#$ -o logs/gene_contact_score.$JOB_ID.out
#$ -pe shared 24
#$ -l h_rt=24:00:00,mem_free=2G           # per-slot memory; adjust per your site’s policy
#$ -m abe
#$ -M hex002@ucsd.edu

set -euo pipefail
mkdir -p logs

echo "Job ID: $JOB_ID"
echo "Host: $(hostname)"
echo "CWD:  $(pwd)"
echo "NSLOTS: ${NSLOTS:-<unset>}"
echo "Date: $(date)"
echo "----"

# Conda (optional)
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
  conda activate schic || true
fi
echo "Conda env: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "----"

# Inputs
dataset="${dataset:-${1:-}}"
if [[ -z "$dataset" ]]; then
  echo "ERROR: dataset not provided."
  echo "Usage:"
  echo "  qsub gene_contact_score.sge -- <dataset>"
  echo "  # or"
  echo "  qsub -v dataset=<dataset> gene_contact_score.sge"
  exit 2
fi
echo "Dataset: $dataset"

# Paths
path="/u/project/jflint/heffel/Heng/hic_impute"
genes_bed="/u/home/h/hex002/project-cluo/BICAN/geneslop2k.bed"
chrom_sizes="/u/home/h/hex002/project-cluo/BICAN/hg38.chrom_1-22.sizes"
out_hdf="/u/project/jflint/heffel/Heng/hic_gene_score/${dataset}.geneimputescore.hdf"
cell_table_txt="${path}/${dataset}_cell_table_10k.txt"
cell_table_tsv="${path}/${dataset}_cell_table_10k.tsv"

threads="${NSLOTS:-30}"

# Checks
[[ -d "${path}/${dataset}/10k" ]] || { echo "Missing: ${path}/${dataset}/10k"; exit 3; }
[[ -f "${genes_bed}" ]]          || { echo "Missing: ${genes_bed}"; exit 4; }
[[ -f "${chrom_sizes}" ]]        || { echo "Missing: ${chrom_sizes}"; exit 5; }

# 1) list .cool files (stable order)
echo "[1/3] Scanning .cool files…"
find "${path}/${dataset}/10k" -type f -name "*.cool" | sort > "${cell_table_txt}"
nfiles=$(wc -l < "${cell_table_txt}" || echo 0)
(( nfiles > 0 )) || { echo "No .cool files found"; exit 6; }

# 2) build TSV: <cell_id>\t<path>
echo "[2/3] Building cell table TSV…"
paste <(awk -F'/' '{print $NF}' "${cell_table_txt}" | sed 's/\.cool$//') \
      "${cell_table_txt}" \
 | sort -k1,1 > "${cell_table_tsv}"

# 3) hicluster gene-score
echo "[3/3] Running hicluster gene-score…"
set -x
hicluster gene-score \
  --cell_table_path "${cell_table_tsv}" \
  --gene_meta_path  "${genes_bed}" \
  --resolution      10000 \
  --output_hdf_path "${out_hdf}" \
  --chrom_size_path "${chrom_sizes}" \
  --cpu             "${threads}" \
  --mode            impute
set +x

echo "Done → ${out_hdf}"
