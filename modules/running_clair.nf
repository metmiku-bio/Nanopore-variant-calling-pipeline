process CONSENSUS_PIPELINE {

    publishDir "${params.outdir}/consensus", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(bed)


    output:
    path("sample_consensus_amplicon_fastas/**/*.fa"), emit: fasta
    path("sample_consensus_amplicon_fastas/**"), emit: all_results

    script:
    """
    bash ${projectDir}/script/new_clair.sh \
        --bam ${bam} \
        -r ${ref} \
        -g ${bed}
    """
}