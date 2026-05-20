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

include {  ALIGN_SAMPLE_ITS ; CALCULATE_COVERAGE_ITS;PLOT_COVERAGE_ITS;CONSENSUS_PIPELINE_ITS;MAKE_CONSENSUS_ITS; CONCATENATE_HAPLOTYPE_ITS;RAXML_NG_TREE_ITS } from './modules/its_major_pipeline.nf'
include {  ALIGN_SAMPLE_COX ; CALCULATE_COVERAGE_COX;PLOT_COVERAGE_COX;CONSENSUS_PIPELINE_COX;MAKE_CONSENSUS_COX;CONCATENATE_HAPLOTYPE ;RAXML_NG_TREE} from './modules/cox_major_pipeline.nf'
include { ALIGN_SAMPLE ; CALCULATE_COVERAGE ;PLOT_COVERAGE;CALL_VARIANTS;ANNOTATE_SNPS;ANNOTATE_INDELS ;CREATE_BAM_LIST ;FILTER_VARIANTS } from './modules/resistance_analysis.nf'
// include { FILTER_VARIANTS } from './modules/filter_variants.nf'
// include { ANNOTATE_SNPS } from './modules/annotate_snps.nf'
// include { ANNOTATE_INDELS } from './modules/annotate_indels.nf'
// include { CALCULATE_COVERAGE  } from './modules/calculate_coverage.nf'
// include {CREATE_BAM_LIST } from './modules/create_bamlist.nf'
// include {CALL_VARIANTS_ITS} from './modules/call_variant_its.nf'
// include {EXTRACTING_READS} from './modules/extracting_its_reads.nf'
// // include {CONSENSUS_PIPELINE_ITS  ;MAKE_CONSENSUS ;MAKE_CONSENSUS_COX} from './modules/running_clair.nf'
// include {CONCAT_FASTAS ; MERGE_FASTAS} from './modules/concatinate_haplotype.nf'
// include {PLOT_COVERAGE} from './modules/plot_coverage.nf'

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



workflow CLAIR {

    /*
     * --------------------------------
     * LOAD SAMPLES
     * --------------------------------
     */

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
     * ============================================
     * ITS PIPELINE
     * ============================================
     */

    its_align = ALIGN_SAMPLE_ITS(
        sample_ch,
        file(params.its_ref),
        params.threads
    )

    its_coverage = CALCULATE_COVERAGE_ITS(
        its_align.alignment,
        file(params.its_bed)
    )

    PLOT_COVERAGE_ITS(
        its_coverage.coverage.collect()
    )

    its_clair = CONSENSUS_PIPELINE_ITS(
        its_align.alignment,
        file(params.its_ref),
        file(params.its_bed)
    )

    its_consensus = MAKE_CONSENSUS_ITS(
        its_clair.phased_file,
        file(params.its_ref),
        file(params.its_bed)
    )
    hap1_ch = its_consensus.hap1.collect()
    hap2_ch = its_consensus.hap2.collect()
    concatenate_its=CONCATENATE_HAPLOTYPE_ITS(hap1_ch,hap2_ch)
    tree_input_ch_its = concatenate_its.hap1_mafft.mix(concatenate_its.hap2_mafft)
    RAXML_NG_TREE_ITS(tree_input_ch_its)

    /*
     * ============================================
     * COX PIPELINE
     * ============================================
     */

    cox_align = ALIGN_SAMPLE_COX(
        sample_ch,
        file(params.cox_ref),
        params.threads
    )

    cox_coverage = CALCULATE_COVERAGE_COX(
        cox_align.alignment,
        file(params.cox_bed)
    )

    PLOT_COVERAGE_COX(
        cox_coverage.coverage.collect()
    )

    cox_clair = CONSENSUS_PIPELINE_COX(
        cox_align.alignment,
        file(params.cox_ref),
        file(params.cox_bed)
    )

    cox_consensus = MAKE_CONSENSUS_COX(
        cox_clair.phased_file,
        file(params.cox_ref),
        file(params.cox_bed)
    )
    hap1_ch = cox_consensus.hap1.collect()
    hap2_ch = cox_consensus.hap2.collect()
    concatenate = CONCATENATE_HAPLOTYPE(hap1_ch,hap2_ch)

    tree_input_ch = concatenate.hap1_mafft.mix(concatenate.hap2_mafft)
    raxml_hap1 = RAXML_NG_TREE(tree_input_ch)
    // raxml_hap2= RAXML_NG_TREE(concatenate.hap2_mafft)
}
workflow {
    
    // Display workflow selection
    log.info "=========================================="
    log.info "WORKFLOW SELECTION: ${params.workflow.toUpperCase()}"
    log.info "=========================================="
    
    // Route to appropriate workflow
    if (params.workflow == "resistance") {
        RESISTANCE_ANALYSIS()
        
    // } else if (params.workflow == "test") {
    //     TEST()
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