#!/usr/bin/env nextflow

/*
 * Coverage Calculation Module
 * Calculates per-base and per-amplicon coverage
 */

process CALCULATE_COVERAGE {
    tag "coverage: ${sample}"
    publishDir "${params.outdir}/coverage/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path bed
    val threads
    
    output:
    tuple val(sample), path("${sample}_coverage_summary.txt"), path("${sample}.per-base.bed.gz")
    
    script:
    """
    # Calculate coverage statistics
    mosdepth -x -b ${bed} ${sample} ${bam}
    
    # Calculate mean coverage per amplicon
    bedtools coverage -a ${bed} -b ${bam} -mean > ${sample}_coverage_mean.txt
    
    # Create summary file
    cat > ${sample}_coverage_summary.txt << EOL
Sample: ${sample}
Total reads: \$(samtools view -c ${bam})
Mean coverage: \$(awk '{sum+=\$NF} END {print sum/NR}' ${sample}_coverage_mean.txt)
EOL
    
    # Rename per-base file for consistency
    mv ${sample}.per-base.bed.gz ${sample}.per-base.bed.gz.tmp || true
    """
}