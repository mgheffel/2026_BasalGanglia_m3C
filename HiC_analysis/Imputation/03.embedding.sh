#!/bin/bash

#SBATCH --job-name=embedding
#SBATCH --output=embedding_%j.out  # %j will be replaced with the job ID
#SBATCH --error=embedding_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25          # Adjust CPU per task
#SBATCH --mem=100G                   # Adjust memory based on your requirements
#SBATCH --time=72:00:00             # Adjust time 
#SBATCH --partition=general       # Adjust partition as necessary



echo "Current environment: $CONDA_DEFAULT_ENV"
path=/tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hic_impute
dataset=$1
mkdir -p /tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hic_embedding/100k_schicluster/${dataset}

#### Get list of imputed cells
> "${path}/${dataset}_cell_table_100k.txt"
for dir in "${path}/${dataset}/100k/chunk*" ; do
    find $dir -type f -name "*.cool" >> "${path}/${dataset}_cell_table_100k.txt"
done
paste <(awk -F'/' '{print $NF}' ${path}/${dataset}_cell_table_100k.txt | cut -d. -f1) ${path}/${dataset}_cell_table_100k.txt | sort -k1,1 > ${path}/${dataset}_cell_table_100k.tsv


hicluster embedding --cell_table_path  ${path}/${dataset}_cell_table_100k.tsv --output_dir /tuba/datasets/Public_Datasets/Luo_BICAN_U01_human_brain_dev/snm3C_3C/hic_embedding/100k_schicluster/${dataset} --dim 50 --dist 1000000 --resolution 100000 --scale_factor 100000 --norm_sig --save_model --save_raw --cpu 25