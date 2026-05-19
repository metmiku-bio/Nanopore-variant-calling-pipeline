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

process MERGE_FASTAS {

    tag "${hap_type}"

    publishDir "${params.outdir}/haplotypes", mode: 'copy'

    input:
    tuple val(hap_type), path(fasta_files)

    output:
    path("combined.${hap_type}.fasta")

    script:
    """
    cat ${fasta_files} > combined.${hap_type}.fasta
    """
}