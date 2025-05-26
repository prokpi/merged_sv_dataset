#!/bin/bash
 
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --error="job-%j.err"
#SBATCH --output="job-%j.out"
#SBATCH --partition=prod
#SBATCH --job-name="process_vcfs"
 
source /shared/conda/miniconda3/etc/profile.d/conda.sh
conda activate pandas
 
echo "Using Python from: $(which python)"
echo "Python version: $(python --version)"
echo "Pandas version: $(python -c 'import pandas as pd; print(pd.__version__)')"
 
INPUT_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/vcf_files"
OUTPUT_FOLDER="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/output"
 
mkdir -p "$OUTPUT_FOLDER"
 
for VCF_PATH in "$INPUT_DIR"/*_merged.vcf; do
    if [[ -f "$VCF_PATH" ]]; then
        SAMPLE_NAME=$(basename "$VCF_PATH" "_merged.vcf")
        OUTPUT_PATH="$OUTPUT_FOLDER/${SAMPLE_NAME}_processed.tsv"
        echo "Processing $VCF_PATH -> $OUTPUT_PATH"
        python3 conc.py --repeat_masker repeat_masker.tsv \
                        --duplications segmental_duplications.tsv \
                        --vcf_file "$VCF_PATH" \
                        --output_file "$OUTPUT_PATH"
    else
        echo "No VCF files found in $INPUT_DIR"
    fi
done
