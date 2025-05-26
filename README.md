# merged_sv_dataset
This repository contains scripts to process genomic VCF files, calculate GC content, Mosdepth coverage and integrate predictions from multiple structural variant callers (Manta, Delly, Smoove). The final output is a consolidated dataset for machine learning, enabling assessment of prediction of structural variations.

## Structure

```
merged_sv_dataset/
│
├── README.md                    
├── files/
│   └── process_vcfs.py          
│   └── process_vcfs.sh          # initial table             
│   └── run_gc_content.sh        # GC content
│   └── run_cov_mosdepth.sh      # mosdepth
│   └── bed_file.py              
│   └── bed_file.sh              # create bed files for each sample
│  
└──    
```

## 1. Create initial table 
```
conda activate mypython
conda activate pandas
chmod +x 'process_vcfs.sh'
sbatch 'process_vcfs.sh'
```

## 2. GC content
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

## 3. Mosdepth (coverage)
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

## 4. Final BED file
Produce bed file for each sample:
```
sbatch bed_file.sh
```

## 5. Merge all bed files into one table:
In the directory of `.bed` files (/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/output_bed)
```
cat *.bed >> merged.bed
```




