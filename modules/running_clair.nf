process CONSENSUS_PIPELINE {

    publishDir "${params.outdir}/consensus", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(bed)

    output:
    path("sample_consensus_amplicon_fastas/**/*.fa"), emit: fasta
    path("sample_consensus_amplicon_fastas/**/*.merge_output.vcf.gz"), emit: merge_output

    script:
    """
    bash ${projectDir}/script/new_clair.sh \
        --bam ${bam} \
        -r ${ref} \
        -g ${bed}
    """
}

process CLAIR_ONLY {
    publishDir "${params.outdir}/clair_only", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(bed)


    output:
    path("sample_consensus_amplicon_fastas/**/*.consensus.fa"),
        emit: fasta,
        optional: true

    path("sample_consensus_amplicon_fastas/**/*.merge_output.vcf.gz"),
        emit: merge_output

    script:
    """
    bash ${projectDir}/script/clair_only.sh \
        --bam ${bam} \
        -r ${ref} \
        -g ${bed}
    
    """
}