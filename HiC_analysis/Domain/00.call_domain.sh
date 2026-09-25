#!/bin/bash

#SBATCH --job-name=domain_insulation
#SBATCH --output=domain_%j.out 
#SBATCH --error=domain_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30     
#SBATCH --mem=240G                   
#SBATCH --time=40:00:00            
#SBATCH --partition=general     

echo "Current environment: ${CONDA_DEFAULT_ENV:-<not set>}"

dataset="$1"
path="/tuba/cndd2/hex002/Luo_Development"

mkdir -p ${path}/hic_domain/${dataset}/


find "${path}/hic_impute/${dataset}/25k" -type f -path "*/chunk*/*.cool" > "${path}/hic_impute/${dataset}_cell_table_25k.txt"
paste <(awk -F'/' '{print $NF}' "${path}/hic_impute/${dataset}_cell_table_25k.txt" | cut -d. -f1) "${path}/hic_impute/${dataset}_cell_table_25k.txt" | sort -k1,1 > "${path}/hic_impute/${dataset}_cell_table_25k.tsv"


## 25kb
hicluster domain --cell_table_path ${path}/hic_impute/${dataset}_cell_table_25k.tsv \
                --output_prefix ${path}/hic_domain/${dataset} \
                --resolution 25000 \
                --window_size 10 \
                --save_count  \
                --cpu 30