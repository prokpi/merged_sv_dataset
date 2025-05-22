# merged_sv_dataset
This repository contains scripts to process genomic VCF files, calculate GC content, Mosdepth coverage and integrate predictions from multiple structural variant callers (Manta, Delly, Smoove). The final output is a consolidated dataset for machine learning, enabling assessment of prediction of structural variations.

## Project Structure

```
Machine_Learning/
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
## Activate packages:
```
conda activate mypython
conda activate pandas
```

## GC content
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






