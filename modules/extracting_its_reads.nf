process EXTRACTING_READS {

    publishDir "${params.outdir}/its_reads/", mode: 'copy'

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("*.fastq.gz"), emit: its_reads    

    script:
    """
    python ${projectDir}/script/demux_nanopore_amplicon.py \
        --fastq ${fastq} \
        --barcodes ${params.barcodes} \
        --max-mismatch ${params.mismatch} \
        --edge-size ${params.edge_size} \
        --log-prefix ${sample}
    rm -r  *.unassigned.fastq
    gzip -k *.fastq
    """
    
}