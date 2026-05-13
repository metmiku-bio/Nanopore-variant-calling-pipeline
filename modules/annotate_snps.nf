#!/usr/bin/env nextflow

/*
 * SNP Annotation Module
 * Annotates SNPs using snpEff
 */

process ANNOTATE_SNPS {
    tag "annotate_snps"
    publishDir "${params.outdir}/annotations", mode: 'copy'
    
    input:
    path filtered_vcf
    val snpeff_db
    val threads
    
    output:
    path "snps.ann.vcf.gz", emit: annotated_snps_vcf
    path "combined_snps_trans.txt", emit: snps_table
    path "snps.vcf.gz", emit: raw_snps_vcf
    
    script:
    """
    # Extract SNPs
    bcftools view --threads ${threads} -v snps ${filtered_vcf} -Oz -o snps.vcf.gz
    bcftools index snps.vcf.gz
    
    # Annotate with snpEff
    java -jar  /mnt/storage13/ahri/snpEff/snpEff.jar ${snpeff_db}  -ud 0 snps.vcf.gz > snps.ann.vcf
    bgzip -f snps.ann.vcf
    bcftools index snps.ann.vcf.gz
    
    # Convert to tabular format
    (echo -e "SAMPLE\\tCHROM\\tPOS\\tREF\\tALT\\tQUAL\\tGT\\tAD\\tDP\\tANN"; \\
        bcftools query -f "[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%GT\\t%AD\\t%DP\\t%ANN\\n]" \\
        snps.ann.vcf.gz) > combined_snps_trans.txt
    
    """
}