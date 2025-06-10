# merged_sv_dataset
This repository contains scripts to process genomic VCF files, calculate GC content, Mosdepth coverage and integrate predictions from multiple structural variant callers (Manta, Delly, Smoove). The final output is a consolidated dataset for machine learning, enabling assessment of prediction of structural variations.

# Structure
- Part 1: Concatenating dataset
- Part 2: Creating a Dockerfile
- Part 3: Nextflow configuration

```
merged_sv_dataset/
│
├── README.md
│   └── requirements.txt         # process_vcfs.py packages 
│   └── Dockerfile
│   └──                  
├── files/
│   └── process_vcfs.py          
│   └── process_vcfs.sh          # initial table             
│   └── run_gc_content.sh        # GC content
│   └── run_cov_mosdepth.sh      # mosdepth
│   └── gc_content.sh            # extracting gc values
│   └── mosdepth_coverage.sh     # extracting mosdepth values
│   └── merge_columns.sh         # merge GC content and Mosdepth results with initial table
│   └── Dockerfile
├── nextflow/
│   └──    
└──    
```

# Part 1: Consolidating dataset

### 1. Create initial table 
```
conda activate mypython
conda activate pandas
chmod +x 'process_vcfs.sh'
sbatch 'process_vcfs.sh'
```

### 2. GC content `run_gc_content.sh`
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

### 3. Mosdepth (coverage) `run_cov_mosdepth.sh`
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

### 4. Extracting values for GC content and Mosdepth into a separate TSV files
Produce files for each sample:
```
sbatch gc_content.sh
sbatch mosdepth_coverage.sh
```

### 5. Merge GC and Mosdepth columns with initial table:
To add extra columns for GC content and Modepth execute this:
```
sbatch merge_columns.sh
```

### 6. Merge all .bed files into one bed file (with all samples):
In the directory of `.bed` file (/shared/work/PI-tommaso.pippucci/RF-WGS/SVs/results/concatenatingfiles/merged)
```
cat *.bed >> merged.bed
```

# Part 2: Creating a Dockerfile
Dockerfile is created (containing `process_vcfs.py` packages) to be pulled using Singularity. 
```
#Create a dockerfile:
docker login
docker tag python-packages kornelipro/python-packages:latest

#push it
docker push kornelipro/python-packages:latest

#pull it using Docker:
docker pull kornelipro/python-packages:latest

#pull it using Singularity:
singularity pull docker://kornelipro/python-packages:latest

#run Dockerfile using Singularity:
singularity run python-packages_latest.sif
```

# Part 3: Nextflow 
Requirements for publishing pipeline:
```
https://nf-co.re/docs/contributing/tutorials/creating_with_nf_core#pipeline-git-repo
```

