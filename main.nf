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
params.mitho = "${projectDir}/reference/sequence.fasta"

// Import modules
// include { INDEX_REFERENCE } from './modules/index_reference.nf'

include { ALIGN_SAMPLE ; ALIGN_SAMPLE_ITS } from './modules/align_sample.nf'
include { CALL_VARIANTS ; FREEBAYES_ONLY } from './modules/call_variants.nf'
include { FILTER_VARIANTS } from './modules/filter_variants.nf'
include { ANNOTATE_SNPS } from './modules/annotate_snps.nf'
include { ANNOTATE_INDELS } from './modules/annotate_indels.nf'
include { CALCULATE_COVERAGE  } from './modules/calculate_coverage.nf'
include {CREATE_BAM_LIST } from './modules/create_bamlist.nf'
include {CALL_VARIANTS_ITS} from './modules/call_variant_its.nf'
include {EXTRACTING_READS} from './modules/extracting_its_reads.nf'
include {CONSENSUS_PIPELINE ; CLAIR_ONLY} from './modules/running_clair.nf'
include {CONCAT_FASTAS ; MERGE_FASTAS} from './modules/concatinate_haplotype.nf'
include {PLOT_COVERAGE} from './modules/plot_coverage.nf'

// Main workflow
workflow RESISTANCE_ANALYSIS {
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
    
    
    // Index reference genome (once)
    // INDEX_REFERENCE(params.ref)
    
    // Process each sample in parallel
    alignment_ch = ALIGN_SAMPLE(sample_ch, params.ref, params.threads)
    
    // Calculate coverage for each sample
    coverage_ch = CALCULATE_COVERAGE(alignment_ch, params.bed)
    plotting = PLOT_COVERAGE(coverage_ch.coverage.collect())
    
    // Collect BAM files for joint calling
// Collect BAM files for joint calling
    bam_files_ch = alignment_ch
        .map { sample, bam, bai -> bam }
        .collect()

    bam_list_ch = CREATE_BAM_LIST(bam_files_ch)
    // Variant calling
    raw_vcf_ch = CALL_VARIANTS(bam_list_ch, params.ref, params.bed, 
                               params.min_base_qual, params.threads)
    
    // // variant calling using freebayes
    //     freebayes_output = FREEBAYES_ONLY(
    //     alignment_ch.alignment,
    //     file(params.ref),
    //     file("${params.ref}.fai"),
    //     file(params.bed)
    // )
    // // variant calling using clair 
    // clair_output = CLAIR_ONLY(
    //     alignment_ch.alignment,
    //     file(params.ref),
    //     file(params.bed)
    // ) 
    // Filter variants
    filtered_vcf_ch = FILTER_VARIANTS(raw_vcf_ch, params.ref, params.bed, params.threads)
    
    // Annotate SNPs and indels
    ANNOTATE_SNPS(filtered_vcf_ch.filtered_vcf, params.snpeff_db, params.threads)
    ANNOTATE_INDELS(filtered_vcf_ch.filtered_vcf, params.snpeff_db, params.threads)
    
    // Optional: Create a summary
    Channel.of("Pipeline completed successfully!")
        | view { message -> println message }
}

// workflow ITS_ANALYSIS {
//     // Read sample IDs from index file
//     Channel.fromPath(params.index_file)
//         .splitCsv(header: true , sep : '\t')
//         .map { row ->
//             tuple(
//                 row.sample,
//                 file(row.fastq)
//             )
//         }
//         .set { sample_ch }
    
//     extracted_its = EXTRACTING_READS(sample_ch)

    
//     // Index reference genome (once)
//     // INDEX_REFERENCE(params.ref)
    
//     // Process each sample in parallel
//     alignment_ch = ALIGN_SAMPLE(extracted_its.its_reads, params.ref, params.threads)
//     coverage_ch = CALCULATE_COVERAGE(alignment_ch.alignment ,params.bed )
//     plotting = PLOT_COVERAGE(coverage_ch.coverage.collect())
//     // Calculate coverage for each sample
//     // coverage_ch = CALCULATE_COVERAGE(alignment_ch.alignment, params.bed)
//     // plotting = PLOT_COVERAGE(coverage_ch.coverage.collect())
    
//     // Collect BAM files for joint calling
// // Collect BAM files for joint calling
//     bam_files_ch = alignment_ch.alignment
//         .map { sample, bam, bai -> bam }
//         .collect()

//     bam_list_ch = CREATE_BAM_LIST(bam_files_ch)
//     bam_files_ch.view()
//     // Variant calling and i remove the bed flag from the freebayes
//     raw_vcf_ch = CALL_VARIANTS_ITS(bam_list_ch, params.ref, params.bed, 
//                                params.min_base_qual, params.threads)
//     raw_vcf_ch.view()
    
//     // // Filter variants and but not use of any bed flag here
//     filtered_vcf_ch = FILTER_VARIANTS(raw_vcf_ch, params.ref, params.bed, params.threads)
//     filtered_vcf_ch.filtered_vcf.view()
//     //     clair_output = CONSENSUS_PIPELINE(
//     //     alignment_ch.alignment,
//     //     file(params.ref),
//     //     file(params.bed)
//     // )
//     //     clair_output_only = CLAIR_ONLY(
//     //     alignment_ch.alignment,
//     //     file(params.ref),
//     //     file(params.bed)
//     // ) 
//     /*
//      * Group fasta files by:
//      *   amplicon + haplotype
//      *
//      * Example:
//      *   BC1.cox1.hap1.fa
//      *   BC2.cox1.hap1.fa
//      *
//      * becomes:
//      *   cox1.hap1
//      */

//     // grouped_fastas = clair_output.fasta
//     //     .flatten()
//     //     .map { fasta ->

//     //         def name = fasta.baseName
//     //         def parts = name.tokenize('.')

//     //         def amp = parts[1]
//     //         def hap = parts[2]

//     //         tuple("${amp}.${hap}", fasta)
//     //     }
//     //     .groupTuple()
//     //     .view()

//     // /*
//     //  * Merge grouped fasta files
//     //  */

//     // CONCAT_FASTAS(grouped_fastas)


//     // MERGE_FASTAS(hap1_ch)
//     // grouped_fastas = clair_output.fasta
//     //     .flatten()
//     //     .map { fasta ->

//     //         def name = fasta.baseName
//     //         def parts = name.tokenize('.')

//     //         def amp = parts[1]
//     //         def hap = parts[2]

//     //         tuple("${amp}.${hap}", fasta)
//     //     }
//     //     .groupTuple()

//     // /*
//     //  * Merge grouped fasta files
//     //  */

//     // CONCAT_FASTAS(grouped_fastas)
    
//     // Optional: Create a summary
//     Channel.of("Pipeline completed successfully!")
//         | view { message -> println message }
// }

// workflow COX_ANALYSIS {
//         Channel.fromPath(params.index_file)
//         .splitCsv(header: true , sep : '\t')
//         .map { row ->
//             tuple(
//                 row.sample,
//                 file(row.fastq)
//             )
//         }
//         .set { sample_ch }

    
//     // Index reference genome (once)
//     // INDEX_REFERENCE(params.ref)
    
//     // Process each sample in parallel
//     alignment_ch = ALIGN_SAMPLE(sample_ch, params.mitho, params.threads)
    
//     // Calculate coverage for each sample
//     coverage_ch = CALCULATE_COVERAGE(alignment_ch,params.bed)
//     plotting = PLOT_COVERAGE(coverage_ch.coverage.collect())

    
//     // Collect BAM files for joint calling
// // Collect BAM files for joint calling
//     bam_files_ch = alignment_ch
//         .map { sample, bam, bai -> bam }
//         .collect()

//     bam_list_ch = CREATE_BAM_LIST(bam_files_ch)
//     // Variant calling and i remove the bed flag from the freebayes
//     raw_vcf_ch = CALL_VARIANTS(bam_list_ch, params.mitho, params.bed, 
//                                params.min_base_qual, params.threads)

    
//     // Filter variants and but not use of any bed flag here
//     filtered_vcf_ch = FILTER_VARIANTS(raw_vcf_ch, params.mitho, params.bed, params.threads)
// }

workflow CLAIR {

    sample_ch = Channel
        .fromPath(params.index_file)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq)
            )
        }

    /*
     * Alignment
     */

    alignment_ch = ALIGN_SAMPLE(
        sample_ch,
        params.ref,
        params.threads
    )
    coverage_ch=CALCULATE_COVERAGE(alignment_ch,params.bed)
    plotting = PLOT_COVERAGE(coverage_ch.coverage.collect())
    /*
     * Consensus generation
     */

    clair_output = CONSENSUS_PIPELINE(
        alignment_ch.alignment,
        file(params.ref),
        file(params.bed)
    )

    /*
     * Group fasta files by:
     *   amplicon + haplotype
     *
     * Example:
     *   BC1.cox1.hap1.fa
     *   BC2.cox1.hap1.fa
     *
     * becomes:
     *   cox1.hap1
     */

    grouped_fastas = clair_output.fasta
        .flatten()
        .map { fasta ->

            def name = fasta.baseName
            def parts = name.tokenize('.')

            def amp = parts[1]
            def hap = parts[2]

            tuple("${amp}.${hap}", fasta)
        }
        .groupTuple()

    /*
     * Merge grouped fasta files
     */

    CONCAT_FASTAS(grouped_fastas)
}

workflow {
    
    // Display workflow selection
    log.info "=========================================="
    log.info "WORKFLOW SELECTION: ${params.workflow.toUpperCase()}"
    log.info "=========================================="
    
    // Route to appropriate workflow
    if (params.workflow == "resistance") {
        RESISTANCE_ANALYSIS()
        
    // } else if (params.workflow == "its") {
    //     ITS_ANALYSIS()
    // } else if (params.workflow == "cox") {
    //     COX_ANALYSIS()
    } else if (params.workflow == "clair") {
        CLAIR()
    } else {
        log.error "Invalid workflow selection: ${params.workflow}"
        log.error "Valid options: 'resistance', 'its', 'cox'"
        exit 1
    }
}