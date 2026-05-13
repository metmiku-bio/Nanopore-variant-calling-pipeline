process CREATE_BAM_LIST {

    publishDir "${params.outdir}/bam_lists", mode: 'copy'

    input:
    val bam_files

    output:
    path "bam_list.txt"

    script:
    def bam_lines = bam_files.collect {
        it.toAbsolutePath().toString()
    }.join('\n')

    """
    printf '%s\n' "${bam_lines}" > bam_list.txt
    """
}