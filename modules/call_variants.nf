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

process FREEBAYES_ONLY {
    conda "${projectDir}/envs/freebayes.yml"
    tag "variant_calling using freebayes for individual samples"
    publishDir "${params.outdir}/freebayes/", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(fai)
    path(bed)

    
    output:
    path "${sample}.freebayes.filtered.vcf.gz", emit: raw_vcf



    script:
    """
# Variant calling
    # Variant calling
    freebayes \
        -f ${ref} \
        -t ${bed} \
        ${bam} \
        --haplotype-length -1 \
        --min-coverage 10 \
        --min-base-quality 20 \
        --min-alternate-fraction 0.1 \
        --ploidy 2 \
        > ${sample}.freebayes.vcf

    # Normalize
    bcftools norm \
        -f ${ref} \
        -m -both \
        ${sample}.freebayes.vcf \
        -Ou \
    | bcftools sort \
        -Oz -o ${sample}.freebayes.norm.sorted.vcf.gz

    # Filter
    bcftools filter -i 'FMT/DP>10' -S .  ${sample}.freebayes.norm.sorted.vcf.gz | \\
        bcftools view --threads 10 -i 'QUAL>20' | \\
            bcftools sort -Oz -o ${sample}.freebayes.filtered.vcf.gz

    # Index filtered VCF
    bcftools index -f ${sample}.freebayes.filtered.vcf.gz
    """
}