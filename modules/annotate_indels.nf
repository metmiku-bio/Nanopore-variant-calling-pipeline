#!/usr/bin/env nextflow

/*
 * Indel Annotation Module
 * Annotates indels using snpEff
 */

process ANNOTATE_INDELS {
    tag "annotate_indels"
    publishDir "${params.outdir}/annotations", mode: 'copy'
    
    input:
    path filtered_vcf
    val snpeff_db
    val threads
    
    output:
    path "indels.ann.vcf.gz", emit: annotated_indels_vcf
    path "combined_indels_trans.txt", emit: indels_table
    path "indels.vcf.gz", emit: raw_indels_vcf
    
    script:
    """
    # Extract indels
    bcftools view --threads ${threads} -v indels ${filtered_vcf} -Oz -o indels.vcf.gz
    bcftools index indels.vcf.gz
    
    # Annotate with snpEff

    java -jar  /mnt/storage13/ahri/snpEff/snpEff.jar ${snpeff_db} indels.vcf.gz > indels.ann.vcf
    bgzip -f indels.ann.vcf
    bcftools index indels.ann.vcf.gz
    
    # Convert to tabular format
    (echo -e "SAMPLE\\tCHROM\\tPOS\\tREF\\tALT\\tQUAL\\tGT\\tAD\\tDP\\tANN"; \\
        bcftools query -f "[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%GT\\t%AD\\t%DP\\t%ANN\\n]" \\
        indels.ann.vcf.gz) > combined_indels_trans.txt
    
    """
}