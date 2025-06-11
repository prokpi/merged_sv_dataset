nextflow.enable.dsl = 2
/*
 * Pipeline SVs detection for SR WGS sequencing data
 * Autor: Emanuela Iovino emanuela.iovino@aosp.bo.it
 * This pipeline relies on Nexflow and it works using Nexflow version >= 21.10.6.5660
 */


if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    🧬  Structural Variant Detection Pipeline for SR-WGS data
    ----------------------------------------------------------
    This pipeline detects, filters, and merges SVs from short-read
    WGS data using multiple tools and custom processing steps.

    Usage :
    ______________________________________________________________________

    Required:
    --sample_alignments_tsv      PATH              Path to a TSV file with two columns: the sample ID and the full (absolute) path to the CRAM file
    --reference_fasta            PATH              Reference genome to which the reads are aligned.
    --singularity_cache          PATH              Path to singularity cache directory
    --delly_exclude_regions_bed	 PATH		       Path to bed regiosn used from Delly
    --smoove_exclude_regions_bed PATH		       Path to bed regions used from Smoove
    --expansion_hunter_variant_catalog_json        Path to STR variant catalog file: choose between ExpansionHunter_variant_catalog.json (174k polymorphic loci) or variant_catalog_genes.json (30 known pathogenic loci)
     
    --account_name               NAME              BC+user name
    --profile           NAME              slurm

    """.stripIndent()
    exit 0
}


println '''
🧬  Structural Variant Detection Pipeline 🧬
-------------------------------------------
Detects structural variants from SR-WGS data using multiple tools.
'''


data = channel
        .fromPath(params.sample_alignments_tsv, type: "file", checkIfExists: true)
        .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
        .map { row -> tuple(row.sample_id, row.alignment_file) }

reference_fasta = file(params.reference_fasta, type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)
delly_exclude_regions_bed = file(params.delly_exclude_regions_bed, type: "file", checkIfExists: true)
smoove_exclude_regions_bed = file(params.smoove_exclude_regions_bed, type: "file", checkIfExists: true)
expansion_hunter_variant_catalog_json = file(params.expansion_hunter_variant_catalog_json, type: "file", checkIfExists: true)
recode_delly_python_script = file(params.recode_delly_python_script, type: "file", checkIfExists: true)
mosdepth_segmental_duplications_bed = file(params.mosdepth_segmental_duplications_bed, type:"file", checkIfExists: true)
samtools_path = params.samtools_path

manta_inv = file(params.manta_inv, type: "file", checkIfExists: true)

// Run mosdepth for getting coverage statistics.
process COVERAGE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file) 
    output:
    path("${sample_id}.mosdepth.summary.txt")
    path("${sample_id}.segdups.regions.bed.gz") 
    script:
    """
    mosdepth --threads $task.cpus --no-per-base --fast-mode --fasta $reference_fasta $sample_id $alignment_file
    mosdepth --threads $task.cpus --by $mosdepth_segmental_duplications_bed --use-median --mapq 30 --fasta $reference_fasta ${sample_id}.segdups $alignment_file
    """
    stub:
    """
    touch ${sample_id}.mosdepth.summary.txt
    """
}

process cnvpytor {                                                              
    publishDir "results/$sample_id/", mode:"copy"                               
    input:                                                                      
    tuple val(sample_id), val(alignment_file)                                   
    output:                                                                     
    tuple val(sample_id), path("cnvpytor/$sample_id*.tsv")                      
    script:                                                                     
    """                                                                         
    cnvpytor -root ${sample_id}.pytor -chrom echo chr{1..22} -rd $alignment_file -T $reference_fasta
    cnvpytor -root ${sample_id}.pytor -his 1000                                 
                                                                                
    cnvpytor -root ${sample_id}.pytor -partition 1000                           
    cnvpytor -root ${sample_id}.pytor  -call 1000 > ${sample_id}.tsv            
    cut -f 2  ${sample_id}.tsv | cnvpytor -root ${sample_id}.pytor  -genotype 1000
    #perl Hg38/cnvnator2VCF.pl -prefix $sample_id -reference GRCh38 ${sample_id}.tsv $reference_fasta > ${sample_id}_CNV.vcf
    """                                                                         
    stub:                                                                       
    """                                                                         
    #touch ${sample_id}_CNV.vcf  
    touch ${sample_id}.tsv                                                
    """                                                                         
}     

// Run the Delly variant caller.
process DELLY {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-delly.bcf")
    script:
    """
    delly call -g $reference_fasta -x $delly_exclude_regions_bed  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
    delly call -g $reference_fasta -v tmp.bcf -x $delly_exclude_regions_bed -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5
    rm tmp.bcf
    """
    stub:
    """
    touch $sample_id-delly.bcf
    """
}
//// Recode the Delly .bcf file and convert it into a .vcf file.
process RECODE {
    input:
    tuple val(sample_id), path(output_delly)
    output:
    tuple val(sample_id), path("$sample_id-delly-recode.vcf")
    script:
    """
    python $recode_delly_python_script $output_delly ${sample_id}-delly-recode.vcf
    """
    stub:
    """
    touch ${sample_id}-delly-recode.vcf
    """
}


// Run the Manta variant caller.
//process MANTA {
//    publishDir "results/$sample_id/", mode: "copy"
//    input:
//    tuple val(sample_id), val(alignment_file)
//    output:
//    tuple val(sample_id), path("$sample_id-manta.vcf.gz")
//    script:
//    """
//    configManta.py --bam $alignment_file --referenceFasta $reference_fasta --runDir run_folder/ 
//    cd run_folder
//    python runWorkflow.py
//    python  $manta_inv $samtools_path $reference_fasta results/variants/diploidSV.vcf.gz > results/variants/diploidSV_inv.vcf 
//    mv results/variants/diploidSV_inv.vcf ../$sample_id-manta.vcf.gz
//    #mv results/variants/diploidSV.vcf.gz ../$sample_id-manta.vcf.gz
//    """
//    stub:
//    """
//    touch $sample_id-diploidSV.vcf.gz
//    """
//}


process MANTA {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-manta.vcf.gz")
    script:
    """
    configManta.py --bam $alignment_file --referenceFasta $reference_fasta --runDir run_folder/ 
    cd run_folder
    python runWorkflow.py
    python  $manta_inv $samtools_path $reference_fasta results/variants/diploidSV.vcf.gz > results/variants/diploidSV_inv.vcf 
    mv results/variants/diploidSV_inv.vcf ../$sample_id-manta.vcf.gz
    """
    stub:
    """
    touch $sample_id-diploidSV.vcf.gz
    """
}

// Run the Smoove variant caller.
process SMOOVE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-smoove.genotyped.vcf.gz")
    script:
    """
    smoove call --name $sample_id --exclude $smoove_exclude_regions_bed --fasta $reference_fasta  -p $task.cpus --duphold --genotype $alignment_file 
    """
    stub:
    """
    touch $sample_id-smoove.genotyped.vcf.gz
    """
}

// Filter VCF with bcftools (filter on SVLEN and PASS flag).
process FILTER {
    publishDir "results/$sample_id/", mode: "copy"  
    input:
    tuple val(sample_id), path(vcf)
    output: 
    tuple val(sample_id), path("${sample_id}_$vcf-filtered.vcf")
    script:
    if (vcf.name =~ /smoove/ )
    """                                                                         
    #cat $sample_id > sample.txt
    echo "$sample_id" > sample.txt 
    bcftools reheader --sample sample.txt -o ${sample_id}_name.vcf  $vcf 
    bcftools view --threads $task.cpus --include "SVLEN>=50 || SVLEN<=-50" ${sample_id}_name.vcf > ${sample_id}_$vcf-filtered.vcf
    """
    else
    """
    
    echo "$sample_id" > sample.txt
    bcftools reheader --sample sample.txt -o ${sample_id}_name.vcf  $vcf 
    bcftools view --thread $task.cpus --include "SVLEN>=50 || SVLEN<=-50" ${sample_id}_name.vcf | bcftools view --include "FILTER='PASS'" > ${sample_id}_$vcf-filtered.vcf
    """
    stub:
    """
    touch ${sample_id}_$vcf-filtered.vcf 
    """
}

// Create INV vcf file from Manta
process INV {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), path(vcf)
    output:
    tuple val(sample_id), path("${sample_id}_$vcf-inv.vcf")
    script:
    """
    bcftools view --threads $task.cpus --include 'INFO/SVTYPE="INV"' $vcf > ${sample_id}_$vcf-inv.vcf
    """
    stub:
    """
    touch ${sample_id}_$vcf-inv.vcf
    """
}

// Run SURVIVOR to merge the within-sample variant calls.
process MERGE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    val(sample_id)
    tuple path(vcf1), path(vcf2), path(vcf3)
    output: 
    tuple val(sample_id), path("${sample_id}_merged.vcf")
    script:
    """
    ls $vcf2    >  vcf_list.txt
    ls $vcf1    >> vcf_list.txt
    ls $vcf3    >> vcf_list.txt
    SURVIVOR merge vcf_list.txt 1000 1 1 1 1 50 ${sample_id}_merged.vcf
    """
    stub:
    """
    echo $vcf2 >  vcf_list.txt
    echo $vcf1 >> vcf_list.txt
    echo $vcf3 >> vcf_list.txt
    sort vcf_list.txt
    touch ${sample_id}_merged.vcf
    """
}


 process EXPANSION_HUNTER {
   publishDir "results/$sample_id/", mode: "copy"
   input:                                                                   
   tuple val(sample_id), val(alignment_file)
   output: 
   path("${sample_id}.expansion_hunter.vcf")
   script:
    """
    ExpansionHunter --reads $alignment_file --reference $reference_fasta --variant-catalog $expansion_hunter_variant_catalog_json  --output-prefix ${sample_id}.expansion_hunter --analysis-mode streaming --threads $task.cpus
    """
    stub:
    """
    touch ${sample_id}.eh.vcf
    """
 }

process SV_COVERAGE_EXTRACT {
                                                                        
    publishDir "results/$sample_id/", mode:"copy"                               
    input:
    tuple val(sample_id), path(merged_vcf)
    output:
    tuple val(sample_id), path("${sample_id}.sv.bed")                                 
    script:
    """
    bcftools query  -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%SVLEN\n' $merged_vcf > ${sample_id}.sv.bed
    
    """
    
    stub:
    """
    ${sample_id}.sv.bed
    """

 
 }

process SV_COVERAGE {
    publishDir "results/$sample_id/", mode: "copy"                              
    input:                                                                      
    tuple val(sample_id), path(vcf_file), path(repeat_masker), path(duplications)
    output:
    path("${sample_id}.sv.regions.bed.gz")                                 
    script:
    """
    mosdepth --threads $task.cpus  --by $sample_sv_bed --use-median --mapq 30  --fasta $reference_fasta ${sample_id}.sv $alignment_file
    """
    stub:
    
    """
    ${sample_id}.sv.regions.bed.gz
    """
 
 }


process CG_CONTENT {
    publishDir "results/$sample_id/", mode: "copy"                              
    input:                                                                      
    tuple val(sample_id), path(sample_sv_bed)
    output:
    path("${sample_id}.sv.CG.txt")                                 
    script:
    """
    bedtools nuc -fi $reference_fasta -bed $sample_sv_bed > ${sample_id}.sv.CG.txt
    """
    stub:
    
    """
    ${sample_id}.sv.CG.txt
    """
 
 }


process INITIAL_TABLE {
    publishDir "results/$sample_id/", mode: "copy"                              
    input:                                                                      
    tuple val(sample_id), path(vcf_file), path(params.repeat_masker), path(params.duplications)
    output:
    path("${sample_id}_processed.tsv")                             
    script:
    """
    python process_vcfs.py \ --repeat_masker $repeat_masker \ --duplications $duplications \ --vcf_file $vcf_file \ --output_file ${sample_id}_processed.tsv
    """
    stub:
    
    """
    ${sample_id}_processed.tsv
    """
 
 }

workflow {
    // Initial table 
    // Get coverage using Mosdepth 
    COVERAGE(data)
    //cnvpytor(data)

    // Run variant calling (using Delly, Manta, and Smoove)
    delly_output = RECODE(DELLY(data))
    manta_output = MANTA(data)
    smoove_output = SMOOVE(data)
    INV_manta = INV(manta_output) 

    // Filter all the VCFs, and then group tuples (by sample ID), then split into 2 channels using multimap (the first channel is the sample ID, and the second channel is the tuple of filtered VCFs for each sample)
    filtered_vcf_output = FILTER(delly_output.mix(manta_output, smoove_output))
    filtered_vcf_output.groupTuple().multiMap{ it -> sample_id: it[0]; vcfs: it[1] }.set{ filtered_vcfs_to_merge }

    merged_vcf_output = MERGE(filtered_vcfs_to_merge.sample_id, filtered_vcfs_to_merge.vcfs)
    sample_alignment_svvcf = data.join(merged_vcf_output)
    sample_sv_bed = SV_COVERAGE_EXTRACT(merged_vcf_output)
    //sample_sv_bed.view { "BED OUTPUT: ${it}" }
    sample_alignment_svbed = data.join(sample_sv_bed)
    SV_COVERAGE(sample_alignment_svbed)
    CG_CONTENT(sample_alignment_svbed)
    //sample_sv_bed_map = sample_sv_bed.map { it.sample_id, it.bed } 
    //sample_cram_bed = data.join(sample_sv_bed_map)
    //SV_COVERAGE(sample_cram_bed)
    //sample_alignment_svvcf_with_bed = sample_alignment_svvcf.join(sample_sv_bed)
    //SV_COVERAGE(sample_alignment_svvcf_with_bed)
    expansion_hunter_output = EXPANSION_HUNTER(data)
}
