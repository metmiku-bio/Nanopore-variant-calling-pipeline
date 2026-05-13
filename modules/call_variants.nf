process CALL_VARIANTS {
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