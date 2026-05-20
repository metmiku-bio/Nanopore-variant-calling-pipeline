process ALIGN_SAMPLE_COX {
    tag "${sample}"
    publishDir "${params.outdir}/cox/alignment/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq) 
    path(ref) 
    val(threads)
    
    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: alignment
    
    script:
    def prefix = sample
    """
    # Alignment
    minimap2 -x map-ont --MD -t ${threads} -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:nanopore' -a ${ref} ${fastq} | \
        samtools sort -@ ${threads} -o ${sample}.sorted.bam -
    
    # Index the BAM file (separate command, not piped)
    samtools index ${sample}.sorted.bam
    # Create a symlink in the work directory to the published location
    # This ensures the path is accessible
    """
}

process CALCULATE_COVERAGE_COX {
    tag "flagstat: ${sample}"
    publishDir "${params.outdir}/cox/flagstat", mode: 'copy'
    
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

process PLOT_COVERAGE_COX {
    conda "/mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/envs/r_plot.yml"

    publishDir "${params.outdir}/cox/plots", mode: 'copy'

    input:
    path coverage_files

    output:
     path "*.pdf" 


    script:
    """
    Rscript ${projectDir}/script/coverage_plot.R ${coverage_files.join(' ')}
    """
}

process CONSENSUS_PIPELINE_COX {

    publishDir "${params.outdir}/cox/consensus", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(bed)

output:
tuple val(sample),
      path("sample_consensus_amplicon_fastas/${sample}/${sample}.phased.vcf.gz"),
      path("sample_consensus_amplicon_fastas/${sample}/${sample}.phased.vcf.gz.tbi"),
      path("sample_consensus_amplicon_fastas/${sample}/${sample}.mask.bed"),
      emit: phased_file
    script:
    """
    bash ${projectDir}/script/update_clair.sh \
        --bam ${bam} \
        -r ${ref} \
        -g ${bed}
    """
}

process MAKE_CONSENSUS_COX {

    publishDir "${params.outdir}/cox/consensus", mode: 'copy'

    input:
    tuple val(sample),
          path(vcf),
          path(tbi),
          path(mask)

    path(ref)
    path(bed)

    output:
    path("${sample}.*.hap1.fa"), emit: hap1
    path("${sample}.*.hap2.fa"), emit: hap2
    script:
    """
    bash ${projectDir}/script/make_consensus_cox.sh \
        ${sample} \
        ${ref} \
        ${bed} \
        ${vcf} \
        ${mask}
    """
}

process CONCATENATE_HAPLOTYPE {

    conda "${projectDir}/envs/mafft.yml"

    publishDir "${params.outdir}/cox/concatenate", mode: 'copy'

    input:
    path(hap1_files)
    path(hap2_files)

    output:
    path("concatenate_hap1_cox_region.fa"), emit: hap1
    path("concatenate_hap2_cox_region.fa"), emit: hap2
    path ("concatenate_hap1_cox_region.afa"),emit:hap1_mafft
    path ("concatenate_hap2_cox_region.afa"),emit:hap2_mafft

    script:
    """
    cat ${hap1_files.join(' ')} > concatenate_hap1_cox_region.fa

    cat ${hap2_files.join(' ')} > concatenate_hap2_cox_region.fa

    mafft --auto concatenate_hap1_cox_region.fa > concatenate_hap1_cox_region.afa
    mafft --auto concatenate_hap2_cox_region.fa > concatenate_hap2_cox_region.afa

    """
}

process RAXML_NG_TREE {

    conda "${projectDir}/envs/raxml.yml"

    publishDir "${params.outdir}/cox/phylogeny", mode: 'copy'

    input:
    path(alignment)

    output:
    path("*.raxml.bestTree"), emit: best_tree
    path("*.raxml.support"), emit: support_tree
    path("*.raxml.log"), emit: log
    path("*.raxml.*"), emit: all_results

    script:

    def prefix = alignment.baseName

    """
    raxml-ng --all \
        --msa ${alignment} \
        --msa-format FASTA \
        --model GTR+G \
        --prefix ${prefix} \
        --seed 826482 \
        --bs-metric tbe \
        --tree rand{1} \
        --bs-trees 1000
    """
}