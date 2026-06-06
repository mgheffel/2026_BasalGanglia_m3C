#!/bin/bash

#SBATCH --job-name=schic_rmbl   # Job name
#SBATCH --output=rmbkl-%j.out    # Standard output and error log (%j will be replaced with the job ID)
#SBATCH --error=rmbkl-%j.err     # Error log
#SBATCH --time=10:00:00          # Time limit (hh:mm:ss)
#SBATCH --partition=general      # Partition (queue) to submit the job
#SBATCH --ntasks=1               # Number of tasks
#SBATCH --cpus-per-task=60       # Number of CPU cores per task
#SBATCH --mem=60G                # Memory pool for the job (e.g., 64 GB)

echo "Current environment: ${CONDA_DEFAULT_ENV:-<not set>}"

dataset="$1"
echo "Processing dataset ${dataset}"
path="/tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C"


output_file="${path}/hic_pairs_lcov_dedup/${dataset}_contact_table.tsv"
counts="${path}/hic_pairs_lcov_dedup/${dataset}_contact_counts.tsv"

# Remove existing outputs if they exist
rm -f "$output_file" "$counts"

echo -e "cellID\tCount" > "$counts"

export output_file
export counts

find "${path}/hic_pairs_lcov_dedup/${dataset}/" -maxdepth 1 -mindepth 1 -type f | \
  parallel -j8 '
    file="{}"
    id=$(basename "${file%%.srt.dedup.lcov.pairs.gz}")
    id="${id//_map3C/}"
    count=$(zgrep -c -v "^#" "$file")
    if [ "$count" -gt 10000 ]; then
        echo -e "${id}\t${file}" >> "$output_file"
        echo -e "${id}\t${count}" >> "$counts"
    fi
  '


# # Create the output directory for the filtered contacts
mkdir -p "${path}/hic_rmbkl/${dataset}"

# Run the hicluster filter-contact command
hicluster filter-contact \
    --output_dir "${path}/hic_rmbkl/${dataset}" \
    --blacklist_1d_path "${path}/hg38-blacklist.v2.bed.gz" \
    --blacklist_2d_path "${path}/hg38_m3C_2d_blacklist.bedpe.gz" \
    --cpu 60 \
    --chr1 1 \
    --pos1 2 \
    --chr2 3 \
    --pos2 4 \
    --contact_table "${output_file}" \
    --chrom_size_path "${path}/hg38.chrom_1-22.sizes"
