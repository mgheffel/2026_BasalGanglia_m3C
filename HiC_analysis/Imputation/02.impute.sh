#!/bin/bash
echo "Current environment: ${CONDA_DEFAULT_ENV:-<not set>}"

dataset="$1"
path="/tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C"
# path="/tuba/cndd2/hex002/Tian_snm3C"
# path="/tuba/cndd2/hex002/Luo_IGVF"

# Find the path for all files
find "${path}/hic_rmbkl/${dataset}" -type f > "${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.txt"

paste <(awk -F'/' '{print $NF}'  ${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.txt | cut -d. -f1)  ${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.txt | sort -k1,1 >  ${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.tsv


# Imputation for 100k
mkdir -p ${path}/hic_impute/${dataset}/100k

hicluster prepare-impute --cell_ ${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.tsv \
                         --batch_size 1536 \
                         --pad 1 \
                         --cpu_per_job 70 \
                         --chr1 1 \
                         --pos1 2 \
                         --chr2 3 \
                         --pos2 4 \
                         --output_dir ${path}/hic_impute/${dataset}/100k \
                         --chrom_size_path "${path}/hg38.chrom_1-22.sizes" \
                         --output_dist 500000000 \
                         --window_size 500000000 \
                         --step_size 500000000 \
                         --resolution 100000


while read -r cmd; do
    sbatch submit_snakemake.sh "$cmd"
done < ${path}/hic_impute/${dataset}/100k/snakemake_cmd.txt


# Imputation for 25k
# mkdir -p ${path}/hic_impute/${dataset}/25k

# hicluster prepare-impute --cell_ ${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.tsv \
#                          --batch_size 1536 \
#                          --pad 2 \
#                          --cpu_per_job 60 \
#                          --chr1 1 \
#                          --pos1 2 \
#                          --chr2 3 \
#                          --pos2 4 \
#                          --output_dir ${path}/hic_impute/${dataset}/25k \
#                          --chrom_size_path "${path}/hg38.chrom_1-22.sizes" \
#                          --output_dist 5050000 \
#                          --window_size 500000000 \
#                          --step_size 500000000 \
#                          --resolution 25000


# while read -r cmd; do
#     sbatch submit_snakemake.sh "$cmd"
# done < ${path}/hic_impute/${dataset}/25k/snakemake_cmd.txt  


# Imputation for 10k
# mkdir -p ${path}/hic_impute/${dataset}/10k
# mkdir -p /tuba/cndd2/hex002/Luo_Development/hic_impute/${dataset}/10k

# hicluster prepare-impute --cell_ ${path}/hic_rmbkl/${dataset}_contact_table_rmbkl.tsv \
#                         --batch_size 700 \
#                         --pad 2 \
#                         --cpu_per_job 70 \
#                         --chr1 1 \
#                         --pos1 2 \
#                         --chr2 3 \
#                         --pos2 4 \
#                         --output_dir /tuba/cndd2/hex002/Luo_Developement/hic_impute/${dataset}/10k \
#                        --chrom_size_path "${path}/hg38.chrom_1-22.sizes" \
#                         --output_dist 5050000 \
#                         --window_size 30000000 \
#                         --step_size 10000000 \
#                         --resolution 10000


# cmdfile=/tuba/cndd2/hex002/Luo_Developement/hic_impute/${dataset}/10k/snakemake_cmd.txt
# while read -r cmd; do
#     sbatch submit_snakemake.sh "$cmd"
# done < "$cmdfile"