#!/usr/bin/env nextflow

/*
 * Variant Calling Module
 * Joint variant calling using FreeBayes in GVCF mode
 */

process CALL_VARIANTS {
    tag "variant_calling"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    path bam_list
    path ref
    path bed
    val min_base_qual
    val threads
    
    output:
    path "combined.genotyped.vcf.gz"
    
    script:
    """
    # Check bam_list is not empty
    if [ ! -s ${bam_list} ]; then
        echo "ERROR: bam_list.txt is empty!"
        exit 1
    fi
    
    # Run FreeBayes
    freebayes -f ${ref} \\
        -t ${bed} \\
        -L ${bam_list} \\
        --haplotype-length -1 \\
        --min-coverage 10 \\
        --min-base-quality ${min_base_qual} \\
        --gvcf > combined.genotyped.vcf
    
    # Normalize and compress
    bcftools view --threads ${threads} -T ${bed} combined.genotyped.vcf | \\
        bcftools norm -f ${ref} | \\
        bcftools sort -Oz -o combined.genotyped.vcf.gz
    
    # Index VCF
    bcftools index combined.genotyped.vcf.gz
    
    # Clean up intermediate file
    rm combined.genotyped.vcf
    """
}