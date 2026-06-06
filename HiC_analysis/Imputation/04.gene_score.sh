#!/bin/bash
#SBATCH --job-name=gene_contact_score
#SBATCH --output=gene_contact_score_%j.out  # %j will be replaced with the job ID
#SBATCH --error=gene_contact_score_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10         # Adjust CPU per task
#SBATCH --mem=1G                   # Adjust memory based on your requirements
#SBATCH --time=1:00:00             # Adjust time (2 hours example)
#SBATCH --partition=general       # Adjust partition as necessary


dataset="$1"
path=/tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hic_impute
ls ${path}/${dataset}/10k/chunk*/*.cool |   awk '{printf("%s\n", $0)}' > "${path}/${dataset}_cell_table_10k.txt"


paste <(awk -F'/' '{print $NF}' ${path}/${dataset}_cell_table_10k.txt | cut -d. -f1) ${path}/${dataset}_cell_table_100k.txt | sort -k1,1 > ${path}/${dataset}_cell_table_10k.tsv
hicluster gene-score --cell_table_path ${path}/${dataset}_cell_table_10k.tsv \
                     --gene_meta_path /tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/gencode.v37.annotation.intragenic.bed \
                     --resolution 10000 \
                     --output_hdf_path /tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hic_gene_score/${dataset}.geneimputescore.hdf \
                     --chrom_size_path /tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hg38.chrom_1-22.sizes \
                     --cpu 70 \
                     --mode impute
