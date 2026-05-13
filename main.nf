#!/usr/bin/env nextflow

/*
 * Amplicon Nanopore Variant Calling Pipeline
 * -------------------------------------------
 * Nextflow implementation with modular components
 */

// Workflow parameters
params.index_file = ""
params.ref = ""
params.gff = ""
params.bed = ""
params.threads = 10
params.min_base_qual = 20
params.snpeff_db = ""
params.outdir = "./results"

// Import modules
// include { INDEX_REFERENCE } from './modules/index_reference.nf'
include { ALIGN_SAMPLE } from './modules/align_sample.nf'
include { CALL_VARIANTS } from './modules/call_variants.nf'
include { FILTER_VARIANTS } from './modules/filter_variants.nf'
include { ANNOTATE_SNPS } from './modules/annotate_snps.nf'
include { ANNOTATE_INDELS } from './modules/annotate_indels.nf'
include { CALCULATE_COVERAGE } from './modules/calculate_coverage.nf'
include {CREATE_BAM_LIST } from './modules/create_bamlist.nf'

// Main workflow
workflow {
    // Read sample IDs from index file
    Channel.fromPath(params.index_file)
        .splitCsv(header: true , sep : '\t')
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq)
            )
        }
        .set { sample_ch }
    
    sample_ch.view()
    
    // Index reference genome (once)
    // INDEX_REFERENCE(params.ref)
    
    // Process each sample in parallel
    alignment_ch = ALIGN_SAMPLE(sample_ch, params.ref, params.threads)
    
    // Calculate coverage for each sample
    coverage_ch = CALCULATE_COVERAGE(alignment_ch, params.bed)
    
    // Collect BAM files for joint calling
// Collect BAM files for joint calling
    bam_files_ch = alignment_ch
        .map { sample, bam, bai -> bam }
        .collect()

    bam_list_ch = CREATE_BAM_LIST(bam_files_ch)
    bam_files_ch.view()
    // Variant calling
    raw_vcf_ch = CALL_VARIANTS(bam_list_ch, params.ref, params.bed, 
                               params.min_base_qual, params.threads)
    
    // Filter variants
    filtered_vcf_ch = FILTER_VARIANTS(raw_vcf_ch, params.ref, params.bed, params.threads)
    
    // Annotate SNPs and indels
    ANNOTATE_SNPS(filtered_vcf_ch, params.snpeff_db, params.threads)
    ANNOTATE_INDELS(filtered_vcf_ch, params.snpeff_db, params.threads)
    
    // Optional: Create a summary
    Channel.of("Pipeline completed successfully!")
        | view { message -> println message }
}