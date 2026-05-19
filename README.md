# Amplicon Nanopore Variant Calling Pipeline

This is a modular Nextflow pipeline designed for amplicon-based variant calling and consensus generation using Oxford Nanopore sequencing data.

## Features

- **Read Alignment:** Maps Nanopore reads to a reference genome.
- **Coverage Calculation & Plotting:** Calculates depth of coverage across amplicons and generates HTML plots/reports (via R and Quarto).
- **Variant Calling & Filtering:** Performs joint variant calling and subsequent filtering.
- **Annotation:** Annotates SNPs and Indels using SnpEff.
- **Consensus Generation (Clair3):** Utilizes Clair3 for accurate variant calling and consensus sequence generation, separating alleles by haplotype.

## Prerequisites

- **Nextflow:** Core execution engine.
- **Docker / Singularity:** For running the custom pipeline environment.

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
nextflow run main.nf \
    --workflow resistance \
    --index_file path/to/index.tsv \
    --ref path/to/reference.fasta \
    --bed path/to/amplicons.bed \
    --snpeff_db "your_snpeff_db" \
    --outdir ./results
```

### 2. Clair Consensus Workflow
This workflow aligns reads, calculates coverage, runs Clair3 to phase and call variants, and extracts grouped consensus FASTA files based on amplicons and haplotypes.

```bash
nextflow run main.nf \
    --workflow clair \
    --index_file path/to/index.tsv \
    --ref path/to/reference.fasta \
    --bed path/to/amplicons.bed \
    --outdir ./results
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
- `modules/`: A directory containing the individual, reusable process definitions (e.g., alignment, variant calling, plotting).
