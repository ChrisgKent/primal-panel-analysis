#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\
    P R I M A L  P A N N E L  A N A L Y S I S   
    =========================================
    """
    .stripIndent()

process minimap_align {
    label 'process_medium'
    container 'biocontainers/minimap2:2.26--he4a0461_1'
    conda 'conda-forge::minimap2==2.26--he4a0461_1'

    input:
        tuple val(unique_id), path(barcode_dir), path(ref)
    
    output:
        tuple val(unique_id), path("${unique_id}_align.sam")

    script:
    """
    minimap2 -a -x map-ont -t $task.cpus ${barcode_dir}/*.fastq.gz ${ref} > ${unique_id}_align.sam
    """
}

process samtools_depth {
    label 'process_single'
    container 'biocontainers/samtools:1.17--hd87286a_1'
    conda 'conda-forge::samtools==1.17--hd87286a_1'
    publishDir "${params.outdir}/${unique_id}_depth", mode: 'copy'

    input:
        tuple val(unique_id), path(align_sam)
    
    output:
        tuple val(unique_id), path("${unique_id}_depth")

    script:
    """
    samtools sort -O bam ${align_sam} | samtools depth - > ${unique_id}_depth 
    """
}


workflow {
    // Get the referance files
    referance_fasta = file(params.referance_fasta)

    // Find all barcodes from fastq_pass
    barcode_ch = Channel.fromPath("${params.fastq_pass}/*", type: 'dir')
        .map{tuple(it.getName(), it, referance_fasta)}

    
    results = minimap_align(barcode_ch) | samtools_depth
    view(results)
        

}
