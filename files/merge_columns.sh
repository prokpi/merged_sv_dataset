#!/bin/bash
#SBATCH --job-name=combine_bed_gc_cov
#SBATCH --partition=prod
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=combine_bed-%j.out
#SBATCH --error=combine_bed-%j.err
 
GC_RESULTS_PATH="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/processed_gc_results"
SV_RESULTS_PATH="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/output" 
MOSDEPTH_PATH="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/mosdepth_extracted"
OUTPUT_PATH="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/concatenated"
 
mkdir -p "$OUTPUT_PATH"
 
for gc_file in "$GC_RESULTS_PATH"/*_merged_gc_content_processed.tsv; do 
    base_name=$(basename "$gc_file" _merged_gc_content_processed.tsv) 
    sv_file="$SV_RESULTS_PATH/${base_name}_processed.tsv"
    mosdepth_file="$MOSDEPTH_PATH/${base_name}_coverage.tsv"
  
    output_file="$OUTPUT_PATH/${base_name}_merged.bed"
    header="SAMPLE_NAME\tCHROM_CALLER\tPOS_CALLER\tEND_CALLER\tSVTYPE_CALLER\tSVLEN_CALLER\tMANTA\tDELLY\tSMOOVE\tOVERLAPS_REPEATS\tOVERLAPS_SEG_DUP\tGC_CONTENT\tCOVERAGE_MOSDEPTH"
    echo -e "$header" > "$output_file"
 
    declare -A gc_data
    while IFS=$'\t' read -r chrom start end gc_value; do
        gc_key="$chrom:$start:$end"
        gc_data["$gc_key"]=$gc_value
    done < "$gc_file"
 
    declare -A coverage_data
    while IFS=$'\t' read -r chrom start end coverage; do
        coverage_key="$chrom:$start:$end"
        coverage_data["$coverage_key"]=$coverage
    done < "$mosdepth_file"
 
    while IFS=$'\t' read -r sample chrom start end sv_type sv_len manta delly smoove overlap_rep seg_dup; do 
        if [[ "$chrom" == "CHROM_CALLER" ]]; then
            continue
        fi
 
        bed_start=$((start - 1))
        gc_key="$chrom:$bed_start:$end"
        coverage_key="$chrom:$bed_start:$end"
 
        gc_content="${gc_data[$gc_key]:-NA}"
        coverage="${coverage_data[$coverage_key]:-NA}"

        echo -e "$sample\t$chrom\t$bed_start\t$end\t$sv_type\t$sv_len\t$manta\t$delly\t$smoove\t$overlap_rep\t$seg_dup\t$gc_content\t$coverage" >> "$output_file"
 
    done < "$sv_file"
 
    unset gc_data
    unset coverage_data
    echo "Processed $base_name"
 
done
