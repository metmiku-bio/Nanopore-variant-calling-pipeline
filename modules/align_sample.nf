process ALIGN_SAMPLE {
    tag "${sample}"
    publishDir "${params.outdir}/resistance/alignment/${sample}", mode: 'copy'
    
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