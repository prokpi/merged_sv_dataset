#!/bin/bash
#SBATCH --job-name=mosdepth_rerun
#SBATCH --partition=prod
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

source /shared/conda/miniconda3/etc/profile.d/conda.sh
conda activate mosdepth_env
 
BASE_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/gc_content"
OUTPUT_DIR_BASE="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/mosdepth_results"
REFERENCE_FASTA="/shared/archive/ngsbo/migrated-from-ngsra/db/trioCEU_1KGP_resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
 
for SAMPLE_DIR in "$BASE_DIR"/*_merged_results; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR" "_merged_results")
    REGION_BED="$SAMPLE_DIR/${SAMPLE_NAME}_merged.bed"
    if [ ! -f "$REGION_BED" ]; then
        echo "BED file for sample ${SAMPLE_NAME} not found at ${REGION_BED}. Skipping."
        continue
    fi

    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${SAMPLE_NAME}_results"
    mkdir -p "$OUTPUT_DIR"
    ALIGNMENT_FILE="/shared/archivenew/PI-tommaso.pippucci/RF-WGS/cram/${SAMPLE_NAME}/${SAMPLE_NAME}.cram"
    if [ ! -f "$ALIGNMENT_FILE" ]; then
        echo "CRAM file for sample ${SAMPLE_NAME} not found at ${ALIGNMENT_FILE}. Skipping."
        continue
    fi

    mosdepth --threads 20 \
             --by "$REGION_BED" \
             --use-median \
             --mapq 30 \
             --fasta "$REFERENCE_FASTA" \
             "$OUTPUT_DIR/${SAMPLE_NAME}" \
             "$ALIGNMENT_FILE"
done
