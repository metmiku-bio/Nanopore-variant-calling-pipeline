#!/usr/bin/env nextflow

/*
 * Index Reference Module
 * Creates sequence dictionary and FASTA index
 */

process INDEX_REFERENCE {
    tag "reference: ${ref}"
    publishDir "${params.outdir}/reference", mode: 'copy'
    
    input:
    path ref
    
    output:
    path "${ref}.fai"
    // path "${ref/\.fa/}.dict"
    
    script:
    """
    # Create FASTA index
    samtools faidx ${ref}
    
    """
}