#!/usr/bin/env nextflow

/*
 * Alignment Module
 * Aligns Nanopore reads using minimap2 and sorts BAM
 */

process ALIGN_SAMPLE {
    tag "sample: ${sample}"
    container 'quay.io/biocontainers/mulled-v2-66534bcbd7031a14d93639a505b8508eec9b9cd2:16124d6233b3a62309f7a7d4a6520fc84af1d222-0'
    publishDir "${params.outdir}/bam/${sample}", mode: 'copy'
    
    input:
    val sample
    path ref
    val threads
    
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    
    script:
    def fastq = file("${sample}.fastq.gz")
    
    """
    # Check if FASTQ exists
    if [ ! -f ${fastq} ]; then
        echo "ERROR: FASTQ file ${fastq} not found!"
        exit 1
    fi
    
    # Align with minimap2 and sort
    minimap2 -x map-ont --MD -t ${threads} \\
        -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:nanopore' \\
        -a ${ref} ${fastq} | \\
        samtools sort -@ ${threads} -o ${sample}.bam -
    
    # Index BAM
    samtools index ${sample}.bam
    
    # Generate flagstat
    samtools flagstat ${sample}.bam > ${sample}.flagstat.txt
    
    # Validate BAM
    samtools quickcheck ${sample}.bam && echo "BAM valid" || echo "BAM corrupted"
    """
}