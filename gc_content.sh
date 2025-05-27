#!/bin/bash
#SBATCH --job-name=gc_content
#SBATCH --partition=prod
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --nodes=1
 
input_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/gc_content_results/"
output_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/processed_gc_results/"
 
for file in "${input_folder}"*_merged_gc_content.txt; do
    awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1, $2, $3, $7}' "$file" > "$output_folder/$(basename "$file" .txt)_processed.tsv"
done
