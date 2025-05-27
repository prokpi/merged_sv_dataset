#!/bin/bash
#SBATCH --job-name=extract_coverage
#SBATCH --partition=prod
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=coverage-%j.out
#SBATCH --error=coverage-%j.err
 
input_base="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/mosdepth_results"
output_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/mosdepth_extracted"
mkdir -p "$output_folder"
 
find "$input_base" -type f -name "*.regions.bed.gz" | while read -r regions_file; do
    sample_name=$(basename "$(dirname "$regions_file")" _results)
    zcat "$regions_file" | awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $3, $5}' > "$output_folder/${sample_name}_coverage.tsv"
done
