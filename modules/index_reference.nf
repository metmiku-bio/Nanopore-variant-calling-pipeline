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
    path "${ref/\.fa/}.dict"
    
    script:
    """
    # Create FASTA index
    samtools faidx ${ref}
    
    # Create sequence dictionary (python approach using fastq2matrix)
    python << 'EOF'
import subprocess as sp
import sys
from fastq2matrix import create_seq_dict
try:
    create_seq_dict("${ref}")
except Exception as e:
    print(f"Error creating sequence dictionary: {e}", file=sys.stderr)
    sys.exit(1)
EOF
    
    # Alternative if fastq2matrix is not available:
    # samtools dict ${ref} > ${ref/\.fa/}.dict
    """
}