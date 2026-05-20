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


process CALCULATE_COVERAGE {
    tag "flagstat: ${sample}"
    publishDir "${params.outdir}/flagstat", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai) 
    path(bed)
    
    output:
    path "${sample}_flagstat.txt", emit: flagstat
    path "${sample}_coverage_mean.txt" , emit:coverage
    
    script:
    """
    # Generate flagstat statistics
    samtools flagstat ${bam} > ${sample}_flagstat.txt
    bedtools coverage -a ${bed} -b ${bam} -mean > ${sample}_coverage_mean.txt
    
    # Display basic stats (optional, for logging)
    echo "=== Flagstat for ${sample} ==="
    cat ${sample}_flagstat.txt
    """
}

process PLOT_COVERAGE {
    conda "/mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/envs/r_plot.yml"

    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path coverage_files

    output:
     path "*.pdf" 


    script:
    """
    Rscript ${projectDir}/script/coverage_plot.R ${coverage_files.join(' ')}
    """
}

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

process CALL_VARIANTS {
    conda "${projectDir}/envs/freebayes.yml"
    tag "variant_calling"
    publishDir "${params.outdir}/variants/raw", mode: 'copy'
    
    input:
    path bam_list
    path ref
    path bed
    val min_base_qual
    val threads
    
    output:
    path "combined.genotyped.vcf.gz", emit: raw_vcf
    
    script:
    """
    # Check if BAM list exists
    if [ ! -f ${bam_list} ]; then
        echo "ERROR: BAM list file ${bam_list} not found"
        exit 1
    fi
    
    # Convert relative paths to absolute paths and remove empty lines
    cat ${bam_list} | while read line; do
        if [ -n "\$line" ]; then
            # Resolve to absolute path if relative
            if [[ "\$line" != /* ]]; then
                realpath "\$line" 2>/dev/null || readlink -f "\$line" 2>/dev/null || echo "\$line"
            else
                echo "\$line"
            fi
        fi
    done | grep -v '^\$' > bam_list_resolved.txt
    
    # Verify BAM files exist
    while read bam_file; do
        if [ ! -f "\$bam_file" ]; then
            echo "ERROR: BAM file not found: \$bam_file"
            exit 1
        fi
    done < bam_list_resolved.txt
    
    echo "Found \$(wc -l < bam_list_resolved.txt) BAM files"
    
    # Run FreeBayes
    freebayes -f ${ref} \\
        -t ${bed} \\
        -L bam_list_resolved.txt \\
        --haplotype-length -1 \\
        --min-coverage 10 \\
        --min-base-quality ${min_base_qual} \\
        --gvcf > combined.genotyped.vcf
    
    # Normalize and compress
    bcftools view --threads ${threads} -T ${bed} combined.genotyped.vcf | \\
        bcftools norm -f ${ref} -Oz -o combined.genotyped.vcf.gz
    
    bcftools view combined.genotyped.vcf -Oz -o tmp.vcf.gz

    bcftools sort tmp.vcf.gz -Oz -o combined.genotyped.vcf.gz

    bcftools index combined.genotyped.vcf.gz
    """
}


process FILTER_VARIANTS {
    tag "filtering"
    publishDir "${params.outdir}/variants/filtered", mode: 'copy'
    
    input:
    path raw_vcf
    path ref
    path bed
    val threads
    
    output:
    path "filtered.vcf.gz", emit: filtered_vcf
    path "filtering_stats.txt", emit: stats
    
    script:
    """
    # Apply filters: DP > 10 in at least one sample and QUAL > 30
    bcftools filter -i 'FMT/DP>10' -S . ${raw_vcf} | \\
        bcftools view --threads ${threads} -i 'QUAL>30' | \\
        bcftools sort -Oz -o filtered.vcf.gz
    
    # Index filtered VCF
    bcftools index filtered.vcf.gz
    
    # Generate filtering stats
    bcftools stats filtered.vcf.gz > filtering_stats.txt
    """
}

process ANNOTATE_SNPS {

    tag "annotate_snps"

    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path filtered_vcf
    val snpeff_db
    val threads

    output:
    path "snps.ann.vcf.gz", emit: annotated_snps_vcf
    path "missense.vcf.gz", emit: missense_vcf
    path "missense_table.txt", emit: missense_table
    path "combined_snps_trans.txt", emit: snps_table
    path "snps.vcf.gz", emit: raw_snps_vcf

    script:
    """
    # Extract SNPs
    bcftools view \
        --threads ${threads} \
        -v snps \
        ${filtered_vcf} \
        -Oz -o snps.vcf.gz

    bcftools index snps.vcf.gz

    # Annotate with snpEff
    java -jar /mnt/storage13/ahri/snpEff/snpEff.jar \
        ${snpeff_db} \
        -ud 0 \
        snps.vcf.gz > snps.ann.vcf

    bgzip -f snps.ann.vcf
    bcftools index snps.ann.vcf.gz

    # Extract only missense variants
    bcftools view \
        -i 'ANN~"missense_variant"' \
        snps.ann.vcf.gz \
        -Oz -o missense.vcf.gz

    bcftools index missense.vcf.gz

    # Full SNP table
    (
        echo -e "SAMPLE\\tCHROM\\tPOS\\tREF\\tALT\\tQUAL\\tGT\\tAD\\tDP\\tANN"

        bcftools query \
            -f "[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%GT\\t%AD\\t%DP\\t%ANN\\n]" \
            snps.ann.vcf.gz

    ) > combined_snps_trans.txt

    # Missense-only table
    (
        echo -e "SAMPLE\\tCHROM\\tPOS\\tREF\\tALT\\tQUAL\\tGT\\tAD\\tDP\\tANN"

        bcftools query \
            -f "[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%GT\\t%AD\\t%DP\\t%ANN\\n]" \
            missense.vcf.gz

    ) > missense_table.txt
    """
}



process ANNOTATE_INDELS {
    tag "annotate_indels"
    publishDir "${params.outdir}/annotations", mode: 'copy'
    
    input:
    path filtered_vcf
    val snpeff_db
    val threads
    
    output:
    path "indels.ann.vcf.gz", emit: annotated_indels_vcf
    path "combined_indels_trans.txt", emit: indels_table
    path "indels.vcf.gz", emit: raw_indels_vcf
    
    script:
    """
    # Extract indels
    bcftools view --threads ${threads} -v indels ${filtered_vcf} -Oz -o indels.vcf.gz
    bcftools index indels.vcf.gz
    
    # Annotate with snpEff

    java -jar  /mnt/storage13/ahri/snpEff/snpEff.jar ${snpeff_db} indels.vcf.gz > indels.ann.vcf
    bgzip -f indels.ann.vcf
    bcftools index indels.ann.vcf.gz
    
    # Convert to tabular format
    (echo -e "SAMPLE\\tCHROM\\tPOS\\tREF\\tALT\\tQUAL\\tGT\\tAD\\tDP\\tANN"; \\
        bcftools query -f "[%SAMPLE\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%GT\\t%AD\\t%DP\\t%ANN\\n]" \\
        indels.ann.vcf.gz) > combined_indels_trans.txt
    
    """
}