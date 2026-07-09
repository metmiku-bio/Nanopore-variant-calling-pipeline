# Amplicon Nanopore Variant Calling Pipeline

This is a modular Nextflow pipeline designed for amplicon-based variant calling and consensus generation using Oxford Nanopore sequencing data.

## Features

- **Read Alignment:** Maps Nanopore reads to a reference genome.
- **Coverage Calculation & Plotting:** Calculates depth of coverage across amplicons and generates pdf plots/reports (via R ).
- **Variant Calling & Filtering:** Performs variant calling using freebayes for the resistance and clair workflows .
- **Annotation:** Annotates SNPs and Indels using SnpEff.
- **Consensus Generation (Clair3):** Utilizes Clair3 for accurate variant calling and consensus sequence generation, separating alleles by haplotype.

## Prerequisites - still underdevelopment the profile is now working on conda
- **Clair3:** installing clair3 is mandatory for running the diversity markers then give the path of the tools and  version of the model  in the script/update_clair.sh 
- **Nextflow:** Core execution engine.
- **Docker / Singularity:** For running the custom pipeline environment. - still under development but can be used with conda environment since integrating clair3 makes it difficult

### Docker Environment

The pipeline relies on a comprehensive Docker image containing all dependencies. It is built upon the official Clair3 GPU image and adds essential tools.

To build the image:
```bash
docker build -t amplicon-pipeline:latest .
```

*Tools included in the image:* `clair3`, `samtools`, `bcftools`, `bedtools`, `freebayes`, `seqkit`, `R` (with `ggplot2`, `dplyr`, `plotly`, `DT`, `rmarkdown`), and `quarto`.

## Input Requirements

1. **Index File (`--index_file`):** A tab-separated values (TSV) file with a header mapping sample IDs to their corresponding FASTQ files.
   ```tsv
   sample	fastq
   Sample_A	/path/to/Sample_A.fastq.gz
   Sample_B	/path/to/Sample_B.fastq.gz
   ```
2. **Reference Genome (`--ref`):** Path to the reference FASTA file.
3. **Target BED File (`--bed`):** Path to the BED file specifying the amplicon coordinates.

## Usage

The pipeline execution is driven by the `--workflow` parameter. Currently, two main workflows are supported:

### 1. Resistance Workflow
This workflow aligns reads, calculates coverage, performs joint variant calling, filters the variants, and annotates them.

```bash
nextflow run ../main.nf     --index_file /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/test_data/sample_file_path.tsv     --ref /mnt/storage13/ahri/Anopheles_stephensi/reference/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna     --gff /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/reference/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.gff     --bed /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/bed/resistance_region.bed     --snpeff_db An_stephen     --outdir ./     -resume     --workflow resistance     -profile conda
```

### 2. Clair Consensus Workflow
This workflow aligns reads, calculates coverage, runs Clair3 to phase and call variants, and extracts grouped consensus FASTA files based on amplicons and haplotypes.

```bash
nextflow run main.nf     --index_file /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/test_data/sample_file_path.tsv     --its_ref /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/reference/its_reference_NW_023405050.1:151506-152137.fasta   --cox_ref /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/reference/sequence.fasta  --gff /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/reference/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.gff        --snpeff_db An_stephen     --outdir /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/its_updated_reference_result      --workflow clair  -profile conda --its_bed /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/bed/its_NW_023405050.1:151506-152137.bed  --cox_bed /mnt/storage13/ahri/Anopheles_stephensi/Nanopore-variant-calling-pipeline/bed/cox.bed
```

## Key Parameters

| Parameter | Description | Default |
| :--- | :--- | :--- |
| `--workflow` | Select the pipeline sub-workflow to run (`resistance` or `clair`). | *Required* |
| `--index_file` | Path to the TSV file mapping samples to read files. | `''` |
| `--ref` | Path to the reference FASTA file. | `''` |
| `--bed` | Path to the BED file with target regions. | `''` |
| `--gff` | Path to the GFF annotation file. | `''` |
| `--outdir` | Output directory for pipeline results. | `./results` |
| `--threads` | Number of CPU threads allocated for parallel processes. | `10` |
| `--min_base_qual`| Minimum base quality threshold for variant calling. | `20` |
| `--snpeff_db` | Name of the SnpEff database to use for annotation. | `''` |

## Pipeline Architecture

- `main.nf`: The main Nextflow script defining the workflows.
- `nextflow.config`: Contains default parameter values, execution profiles, and resource requests (CPUs, memory).
- `modules/`: A directory containing the individual, reusable process definitions (e.g., alignment, variant calling, plotting). The main modules are the cox_major_pipeline.nf , its_major_pipeline.nf,and resistance_analysis.nf