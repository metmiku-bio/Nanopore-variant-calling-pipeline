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

// Set container images (optional - uncomment if using containers)
// process.container = "biocontainers/fastq2matrix:latest"
// process.container = "biocontainers/freebayes:latest"
// docker.enabled = true

// Import modules
include { INDEX_REFERENCE } from './modules/index_reference.nf'
include { ALIGN_SAMPLE } from './modules/align_sample.nf'
include { CALL_VARIANTS } from './modules/call_variants.nf'
include { FILTER_VARIANTS } from './modules/filter_variants.nf'
include { ANNOTATE_SNPS } from './modules/annotate_snps.nf'
include { ANNOTATE_INDELS } from './modules/annotate_indels.nf'
include { CALCULATE_COVERAGE } from './modules/calculate_coverage.nf'
include { GENERATE_REPORTS } from './modules/generate_reports.nf'

// Main workflow
workflow {
    // Read sample IDs from index file
    Channel.fromPath(params.index_file)
        | splitCsv(header: true)
        | map { row -> row.sample }
        | filter { it != "" }
        | set { sample_ch }
    
    // Index reference genome (once)
    INDEX_REFERENCE(params.ref)
    
    // Process each sample in parallel
    alignment_ch = ALIGN_SAMPLE(sample_ch, params.ref, params.threads)
    
    // Calculate coverage for each sample
    coverage_ch = CALCULATE_COVERAGE(alignment_ch, params.bed, params.threads)
    
    // Collect BAM files for joint calling
    alignment_ch
        | map { sample, bam, bai -> bam }
        | collectFile(name: "bam_list.txt", newLine: true)
        | set { bam_list_ch }
    
    // Variant calling
    raw_vcf_ch = CALL_VARIANTS(bam_list_ch, params.ref, params.bed, 
                               params.min_base_qual, params.threads)
    
    // Filter variants
    filtered_vcf_ch = FILTER_VARIANTS(raw_vcf_ch, params.ref, params.bed, params.threads)
    
    // Separate SNPs and indels
    filtered_vcf_ch
        | map { vcf -> tuple(vcf, "snps") }
        | ANNOTATE_SNPS(params.snpeff_db, params.threads)
        | set { snp_results }
    
    filtered_vcf_ch
        | map { vcf -> tuple(vcf, "indels") }
        | ANNOTATE_INDELS(params.snpeff_db, params.threads)
        | set { indel_results }
    
    // Generate final reports
    GENERATE_REPORTS(snp_results, indel_results, coverage_ch, params.outdir)
}