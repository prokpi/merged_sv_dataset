#!/bin/bash
#SBATCH --job-name=gc_content
#SBATCH --partition=prod
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=nextflow-%j.out 
#SBATCH --error=nextflow-%j.err 
#SBATCH --nodes=1
 
FILE_LIST="/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/ndd2024.list"
 
while IFS= read -r VCF_FILE; do
    SAMPLE_NAME=$(basename "$VCF_FILE" .vcf)
    OUTPUT_DIR="/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/gc_content/${SAMPLE_NAME}_results"
 
    mkdir -p "$OUTPUT_DIR"
 
    #Converting VCF to BED
    bcftools query -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%SVLEN\n' "$VCF_FILE" > "$OUTPUT_DIR/$SAMPLE_NAME.bed"
 
    #Calculating GC content
    bedtools nuc -fi /shared/archive/ngsbo/migrated-from-ngsra/db/trioCEU_1KGP_resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
                 -bed "$OUTPUT_DIR/$SAMPLE_NAME.bed" > "$OUTPUT_DIR/$SAMPLE_NAME_gc_content.txt"
 
done < "$FILE_LIST"
