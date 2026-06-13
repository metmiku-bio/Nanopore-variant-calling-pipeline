#!/usr/bin/env bash
set -euo pipefail

sample=$1
ref=$2
bed=$3
vcf=$4
mask=$5


region="NC_028223.1:1536-2031"
name="cox1_region"
    for H in 1 2; do

        samtools faidx "$ref" "$region" | \
        bcftools consensus \
            -H $H \
            -m "$mask" \
            "$vcf" \
        | sed "1s/^>.*/>${sample}_${name}_hap${H}/" \
        > "${sample}.${name}.hap${H}.fa"

    done

