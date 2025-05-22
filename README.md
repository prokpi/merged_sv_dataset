# merged_sv_dataset
This repository contains scripts to process genomic VCF files, calculate GC content, Mosdepth coverage and integrate predictions from multiple structural variant callers (Manta, Delly, Smoove). The final output is a consolidated dataset for machine learning, enabling assessment of prediction of structural variations.

## Structure

```
merged_sv_dataset/
│
├── README.md                    
├── files/                        
│   └── run_gc_content.sh        # GC content
│   └── run_cov_mosdepth.sh      # mosdepth
│   └── conc.py                  # concatenating 
│   └── slurm.sh                 # run conc.py
├── notebook/                   
│   └── water_quality_analysis.ipynb   # The code that was ran in Google Colab
└── requirements.txt             # List of dependencies (for installing via pip)
```
## Requirements
- Packages that will need to be activated:
```
conda activate mypython
conda activate pandas
......
```
- Submitting `.sh' files to slurm:
```
chmod +x 'filename.sh'
sbatch 'filename.sh'
```

## 1. GC content
- Files to be processed are listed here: `ndd2024.list`
- Each file is VCF file

bcftools query: creates `.bed` file containing structural variants in the tabular format with the following columns:
```
<CHROM> <POS0> <END> <SVTYPE> <SVLEN>
```
bedtools query: calculates **GC** content for each genomic region specified in the `.bed` file, output:
```
#1-based_Start: The start position (1-based).
#1-based_End: The end position (1-based).
GC%: The percentage of guanine (G) and cytosine (C) nucleotides.
```

## 2. Mosdepth (coverage)
Input requirements:
- `.bed` files produced by `run_gc_content.sh` with genomic regions 
- CRAM files containing aligned sequencing reads for each sample.
- Reference FASTA file for the genome assembly (GRCh38)

This script executes Mosdepth with following parameters:
```
--threads 20
--by "$REGION_BED"
--use-median
--mapq 30
--fasta "$REFERENCE_FASTA"
```
Output: Coverage statistics for each genomic region, including median coverage, are generated in the specified output directory.






