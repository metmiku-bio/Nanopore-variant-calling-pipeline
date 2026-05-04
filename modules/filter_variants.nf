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
    path "combined.genotyped_filtered_FMTDP10.vcf.gz"
    
    script:
    """
    # Apply filters: DP > 10 in at least one sample and QUAL > 30
    bcftools filter -i 'FMT/DP>10' -S . ${raw_vcf} | \\
        bcftools view --threads ${threads} -i 'QUAL>30' | \\
        bcftools sort | \\
        bcftools norm -m - -Oz -o combined.genotyped_filtered_FMTDP10.vcf.gz
    
    # Index filtered VCF
    bcftools index combined.genotyped_filtered_FMTDP10.vcf.gz
    
    # Generate filtering stats
    bcftools stats combined.genotyped_filtered_FMTDP10.vcf.gz > filtering_stats.txt
    """
}