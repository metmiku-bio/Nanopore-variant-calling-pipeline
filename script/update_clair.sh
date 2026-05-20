#!/usr/bin/env bash
set -euxo pipefail

# ============================================================
# ONT Amplicon Consensus Pipeline
#
# Steps:
#   1. Clair3 variant calling
#   2. WhatsHap phasing
#   3. Low-depth masking
#   4. No-call masking
#   5. Consensus generation
#
# Output:
#   sample_consensus_amplicon_fastas/<sample>/
#
# ============================================================

# ------------------------------------------------------------
# Default parameters
# ------------------------------------------------------------

BAM=""
REFERENCE=""
REGIONS_BED=""

THREADS=${THREADS:-8}
MASK_DP=${MASK_DP:-5}

CLAIR3_MODEL=${CLAIR3_MODEL:-/mnt/storage13/ahri/miniforge3/envs/clair3_env/bin/models/r1041_e82_400bps_sup_v500}

CLAIR3_BIN=${CLAIR3_BIN:-/mnt/storage13/ahri/miniforge3/envs/muhaps-clair3/bin/run_clair3.sh}

OUTDIR="sample_consensus_amplicon_fastas"

# ------------------------------------------------------------
# Usage
# ------------------------------------------------------------

usage() {
    echo
    echo "Usage:"
    echo "  $0 --bam sample.sorted.bam -r reference.fa -g regions.bed"
    echo
    exit 1
}

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bam)
            BAM="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -g|--regions)
            REGIONS_BED="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            usage
            ;;
    esac
done

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

[[ -z "$BAM" ]] && usage
[[ -z "$REFERENCE" ]] && usage
[[ -z "$REGIONS_BED" ]] && usage

[[ ! -f "$BAM" ]] && { echo "ERROR: BAM not found"; exit 1; }
[[ ! -f "$REFERENCE" ]] && { echo "ERROR: Reference not found"; exit 1; }
[[ ! -f "$REGIONS_BED" ]] && { echo "ERROR: BED not found"; exit 1; }

# ------------------------------------------------------------
# Fix BED line endings
# ------------------------------------------------------------

sed -i 's/\r$//' "$REGIONS_BED"

# ------------------------------------------------------------
# Prepare reference indexes
# ------------------------------------------------------------

echo "📚 Preparing reference indexes"

[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

DICT="${REFERENCE%.*}.dict"

if [[ ! -f "$DICT" ]]; then
    gatk CreateSequenceDictionary \
        -R "$REFERENCE"
fi

# ------------------------------------------------------------
# Sample naming
# ------------------------------------------------------------

sample=$(basename "$BAM" .sorted.bam)

sample_dir="${OUTDIR}/${sample}"

mkdir -p "$sample_dir"

echo "🧬 Processing sample: $sample"

# ------------------------------------------------------------
# Ensure BAM index exists
# ------------------------------------------------------------

[[ -f "${BAM}.bai" ]] || samtools index "$BAM"

# ============================================================
# 1. Clair3 Variant Calling
# ============================================================

echo "🔬 Running Clair3"

"$CLAIR3_BIN" \
    --bam_fn="$BAM" \
    --ref_fn="$REFERENCE" \
    --threads="$THREADS" \
    --platform=ont \
    --output="${sample_dir}/${sample}.clair3" \
    --model_path="$CLAIR3_MODEL" \
    --include_all_ctgs \
    --bed_fn="$REGIONS_BED" \
    --sample_name="$sample"

# ============================================================
# 2. Rename Clair3 outputs
# ============================================================

echo "📝 Renaming Clair3 outputs"

clair_dir="${sample_dir}/${sample}.clair3"

for f in merge_output pileup full_alignment; do

    src_vcf=$(find "$clair_dir" -name "${f}.vcf.gz" | head -1 || true)

    if [[ -n "${src_vcf:-}" && -f "$src_vcf" ]]; then

        cp "$src_vcf" \
           "${sample_dir}/${sample}.${f}.vcf.gz"

        if [[ -f "${src_vcf}.tbi" ]]; then
            cp "${src_vcf}.tbi" \
               "${sample_dir}/${sample}.${f}.vcf.gz.tbi"
        fi
    fi
done

# ============================================================
# 3. Select final VCF
# ============================================================

in_vcf="${sample_dir}/${sample}.merge_output.vcf.gz"

if [[ ! -s "$in_vcf" ]]; then
    echo "ERROR: Missing merge_output VCF"
    exit 1
fi

tabix -f "$in_vcf" 2>/dev/null || \
bcftools index -t -f "$in_vcf"

# ============================================================
# 4. Filter multiallelic / symbolic sites
# ============================================================

echo "🧹 Filtering variants"

filtered_vcf="${sample_dir}/${sample}.filtered.vcf.gz"

bcftools view \
    -m2 -M2 \
    -e 'ALT="*" || ALT="<NON_REF>"' \
    "$in_vcf" \
    -Oz -o "$filtered_vcf"

bcftools index -t -f "$filtered_vcf"

# ============================================================
# 5. Determine sample name
# ============================================================

vcf_sample=$(bcftools query -l "$filtered_vcf" | head -1 || true)

sample_arg=()

if [[ -n "${vcf_sample:-}" ]]; then
    sample_arg=(-s "$vcf_sample")
fi

# ============================================================
# 6. WhatsHap phasing
# ============================================================

echo "🧩 Running WhatsHap"

phased_vcf="${sample_dir}/${sample}.phased.vcf.gz"

if [[ -n "${vcf_sample:-}" ]]; then

    whatshap phase \
        --reference "$REFERENCE" \
        --indels \
        --sample "$vcf_sample" \
        --ignore-read-groups \
        -o "$phased_vcf" \
        "$filtered_vcf" \
        "$BAM" \
    || {

        echo "⚠️ WhatsHap failed — using unphased VCF"

        cp "$filtered_vcf" "$phased_vcf"
    }

else

    echo "⚠️ Could not determine VCF sample"

    cp "$filtered_vcf" "$phased_vcf"
fi

tabix -f "$phased_vcf" 2>/dev/null || \
bcftools index -t -f "$phased_vcf"

# ============================================================
# 7. Generate low-depth mask
# ============================================================

echo "📉 Generating low-depth mask"

{
    T="$MASK_DP"

    samtools depth -a "$BAM" \
    | awk -v T="$T" 'BEGIN{OFS="\t"}
    {
        if ($3<T) {
            if (s=="") {
                chr=$1
                s=$2-1
                e=$2
            }
            else if ($1==chr && $2==e+1) {
                e=$2
            }
            else {
                print chr,s,e
                chr=$1
                s=$2-1
                e=$2
            }
        }
        else if (s!="") {
            print chr,s,e
            s=""
        }
    }
    END{
        if (s!="")
            print chr,s,e
    }'

} > "${sample_dir}/${sample}.depthmask.bed"

# ============================================================
# 8. Generate no-call mask
# ============================================================

echo "🚫 Generating no-call mask"

if [[ -n "${vcf_sample:-}" ]]; then

    bcftools query \
        -s "$vcf_sample" \
        -i 'GT=".|." || GT="./."' \
        -f '%CHROM\t%POS0\t%POS\n' \
        "$phased_vcf" \
    > "${sample_dir}/${sample}.nocall.mask.bed"

else

    touch "${sample_dir}/${sample}.nocall.mask.bed"
fi

# ============================================================
# 9. Merge masks
# ============================================================

echo "🧬 Merging masks"

cat \
    "${sample_dir}/${sample}.depthmask.bed" \
    "${sample_dir}/${sample}.nocall.mask.bed" \
| sort -k1,1 -k2,2n \
| bedtools merge \
> "${sample_dir}/${sample}.mask.bed"

# ============================================================
# 10. Build haplotype consensus
# ============================================================

