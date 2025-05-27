#!/bin/bash
#SBATCH --job-name=merge_gc_cov
#SBATCH --partition=prod
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=merge-%j.out
#SBATCH --error=merge-%j.err

processed_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/output"
gc_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/processed_gc_results"
cov_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/mosdepth_extracted"
output_folder="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/merged"
mkdir -p "$output_folder"

for processed_file in "$processed_folder"/*_processed.tsv; do
    sample=$(basename "$processed_file" _processed.tsv)
    gc_file="$gc_folder/${sample}_merged_gc_content_processed.tsv"
    cov_file="$cov_folder/${sample}_coverage.tsv"
    output_file="$output_folder/${sample}_merged.tsv"

    if [[ -f "$gc_file" && -f "$cov_file" ]]; then
        awk -F'\t' -v OFS='\t' \
        -v gcfile="$gc_file" -v covfile="$cov_file" \
        '
        BEGIN {
            while ((getline < gcfile) > 0) {
                key = $1 "\t" $2 "\t" $3
                gc[key] = $4
            }
            close(gcfile)

            while ((getline < covfile) > 0) {
                key = $1 "\t" $2 "\t" $3
                cov[key] = $4
            }
            close(covfile)
        }
        NR == 1 {
            print $0, "GC_CONTENT", "COVERAGE_MOSDEPTH"
            next
        }
        {
            # Adjust start coordinate (column 3) from 1-based to 0-based for matching
            key = $2 "\t" ($3 - 1) "\t" $4
            gc_val = (key in gc) ? gc[key] : "NA"
            cov_val = (key in cov) ? cov[key] : "NA"
            print $0, gc_val, cov_val
        }
        ' "$processed_file" > "$output_file"
    else
        echo "Missing GC or coverage files for sample $sample"
    fi
done
