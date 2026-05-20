process CONSENSUS_PIPELINE_ITS {

    publishDir "${params.outdir}/consensus", mode: 'copy'

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


process MAKE_CONSENSUS {

    publishDir "${params.outdir}/consensus", mode: 'copy'

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
    path("${sample}.*.fa"), emit: fasta

    script:
    """
    bash ${projectDir}/script/make_consensus.sh \
        ${sample} \
        ${ref} \
        ${bed} \
        ${vcf} \
        ${mask}
    """
}

process MAKE_CONSENSUS_COX {

    publishDir "${params.outdir}/consensus", mode: 'copy'

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
    path("${sample}.*.fa"), emit: fasta

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