#!/usr/bin/env bash
set -euo pipefail

sample=$1
ref=$2
bed=$3
vcf=$4
mask=$5


region="NW_023405050.1:4-511"
name="its_region"
    for H in 1 2; do

        samtools faidx "$ref" "$region" | \
        bcftools consensus \
            -H $H \
            -m "$mask" \
            "$vcf" \
        | sed "1s/^>.*/>${sample}_${name}_hap${H}/" \
        > "${sample}.${name}.hap${H}.fa"

    done

