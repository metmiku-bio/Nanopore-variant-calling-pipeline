process ALIGN_SAMPLE {
    tag "${sample}"
    publishDir "${params.outdir}/alignment/${sample}", mode: 'copy'
    
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

process ALIGN_SAMPLE_ITS {

    tag "${sample}"
    publishDir "${params.outdir}/alignment_its/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq)
    path(ref)
    val(threads)

    output:
    tuple val(sample),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai"),
          emit: alignment

    path("${sample}_region.bed"), emit: bed

    script:
    """
    minimap2 -x map-ont --MD -t ${threads} \
        -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:nanopore" \
        -a ${ref} ${fastq} | \
    samtools sort -@ ${threads} -o ${sample}.tmp.sorted.bam -

    # INDEX FIRST
    samtools index ${sample}.tmp.sorted.bam

    samtools view -@ ${threads} \
        -b \
        -F 260 \
        -q 10 \
        ${sample}.tmp.sorted.bam \
        NW_023405169.1 \
        > ${sample}.tmp.filtered.bam

    samtools sort -@ ${threads} \
        -o ${sample}.sorted.bam \
        ${sample}.tmp.filtered.bam

    samtools index ${sample}.sorted.bam

    bedtools bamtobed \
        -i ${sample}.sorted.bam \
        > ${sample}_region.bed
    """
}