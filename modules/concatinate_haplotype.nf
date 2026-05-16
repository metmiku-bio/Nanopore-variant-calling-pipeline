process CONCAT_FASTAS {

    publishDir "${params.outdir}/combined_fastas", mode: 'copy'

    input:
    tuple val(group), path(fastas)

    output:
    path("combined.${group}.fasta")

    script:
    """
    cat ${fastas.join(' ')} > combined.${group}.fasta
    """
}
