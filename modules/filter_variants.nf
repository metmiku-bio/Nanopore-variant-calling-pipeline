#!/usr/bin/env nextflow

/*
 * Filter Variants Module
 * Applies quality filters to VCF
 */

process FILTER_VARIANTS {
    tag "filtering"
    publishDir "${params.outdir}/variants/filtered", mode: 'copy'
    
    input:
    path raw_vcf
    path ref
    path bed
    val threads
    
    output:
    path "filtered.vcf.gz", emit: filtered_vcf
    path "filtering_stats.txt",emit: filtered_stats
    // path "filtering_stats.txt", emit: stats
    
    script:
    """
    # Apply filters: DP > 10 in at least one sample and QUAL > 30
    bcftools filter -i 'FMT/DP>10' -S . ${raw_vcf} | \\
        bcftools view --threads ${threads} -i 'QUAL>30' | \\
        bcftools sort -Oz -o filtered.vcf.gz
    
    # Index filtered VCF
    bcftools index filtered.vcf.gz
    
    # Generate filtering stats
    bcftools stats filtered.vcf.gz > filtering_stats.txt
    """
}