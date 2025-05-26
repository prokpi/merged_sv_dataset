#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="final_bed"
#SBATCH --error="job-%j.err"
#SBATCH --output="job-%j.out"
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

TSV_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/output"
GC_BASE_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/gc_content"
COV_BASE_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/mosdepth_results"
OUT_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/output_bed"
mkdir -p "$OUT_DIR"
for tsv in "$TSV_DIR"/*_processed.tsv; do
    sample=$(basename "$tsv" "_processed.tsv")
    gc_dir="$GC_BASE_DIR/${sample}_merged_results"
    gc_bed="$gc_dir/${sample}_merged.bed"
    cov_dir="$COV_BASE_DIR/${sample}_results"
    cov_bed_gz="$cov_dir/${sample}.regions.bed.gz"
    if [[ ! -f "$gc_bed" ]]; then
        echo "WARNING: GC file not found for sample $sample: $gc_bed"
        continue
    fi
    if [[ ! -f "$cov_bed_gz" ]]; then
        echo "WARNING: Coverage file not found for sample $sample: $cov_bed_gz"
        continue
    fi
    cov_bed="$cov_dir/${sample}.regions.bed"
    if [[ ! -f "$cov_bed" ]]; then
        echo "Uncompressing coverage for sample $sample"
        zcat "$cov_bed_gz" > "$cov_bed"
    fi
    output_bed="$OUT_DIR/${sample}.bed"
    python3 bed_file.py "$tsv" "$gc_bed" "$cov_bed" "$output_bed"
done
