#!/usr/bin/env nextflow

/*
 * Coverage Calculation Module
 * Generates flagstat statistics for BAM files
 */

process CALCULATE_COVERAGE {
    tag "flagstat: ${sample}"
    publishDir "${params.outdir}/flagstat", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai) 
    path(bed)
    
    output:
    path "${sample}_flagstat.txt", emit: flagstat
    path "${sample}_coverage_mean.txt" , emit:coverage
    
    script:
    """
    # Generate flagstat statistics
    samtools flagstat ${bam} > ${sample}_flagstat.txt
    bedtools coverage -a ${bed} -b ${bam} -mean > ${sample}_coverage_mean.txt
    
    # Display basic stats (optional, for logging)
    echo "=== Flagstat for ${sample} ==="
    cat ${sample}_flagstat.txt
    """
}