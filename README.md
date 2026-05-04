# Amplicon Nanopore Variant Calling Pipeline

A Nextflow pipeline for amplicon-based variant calling from Nanopore sequencing data.

## Pipeline Overview

This pipeline takes Nanopore FASTQ files and a sample index, aligns the reads to a reference genome, calls and filters variants (SNPs and Indels), annotates them, calculates coverage, and generates final reports.

**Pipeline Steps:**
1. **Reference Indexing**: Prepares the reference genome (`INDEX_REFERENCE`)
2. **Alignment**: Minimap2 alignment of Nanopore reads and sorting with samtools (`ALIGN_SAMPLE`)
3. **Variant Calling**: Calls raw variants (`CALL_VARIANTS`)
4. **Variant Filtering**: Filters variants based on minimum base quality and other metrics (`FILTER_VARIANTS`)
5. **Annotation**: Annotates SNPs and Indels using SnpEff (`ANNOTATE_SNPS`, `ANNOTATE_INDELS`)
6. **Coverage Calculation**: Determines sequencing coverage across amplicon regions (`CALCULATE_COVERAGE`)
7. **Reporting**: Generates summary statistics and reports (`GENERATE_REPORTS`)

## Prerequisites
- [Nextflow](https://www.nextflow.io/) (>= 20.04.0)
- Container engine (Docker/Singularity/Conda) if using standard bioinformatics containers. 

## Usage

### 1. Prepare your Sample Index
Create a CSV file containing your sample information. The header must include a column named `sample`. The pipeline will automatically parse this file and filter out empty entries.

### 2. Run the Pipeline
Execute the pipeline from the command line, providing all the necessary paths and parameters:

```bash
nextflow run main.nf \
    --index_file path/to/samples.csv \
    --ref path/to/reference.fasta \
    --gff path/to/annotation.gff \
    --bed path/to/amplicons.bed \
    --snpeff_db my_snpeff_database \
    --outdir ./results \
    -resume
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--index_file` | `''` | CSV file containing sample IDs |
| `--ref` | `''` | Path to the reference FASTA genome |
| `--gff` | `''` | Path to the GFF annotation file |
| `--bed` | `''` | Path to the BED file defining amplicon target regions |
| `--threads` | `10` | Number of threads to use for parallelizable processes (e.g., Minimap2) |
| `--min_base_qual`| `20` | Minimum base quality threshold for variant calling/filtering |
| `--snpeff_db` | `''` | SnpEff database name for variant annotation |
| `--outdir` | `./results` | Output directory for pipeline results |

## Resource Management

The pipeline handles retries and resource scaling dynamically:
- Default tasks receive 8 GB RAM and standard time allocations.
- If a task fails and is retried, Nextflow doubles the memory and time allocation automatically.
- Tasks labeled `high_memory` start with 32 GB RAM.

## Output Structure

Results are stored in the directory defined by `--outdir` (default: `./results`). Some expected outputs include:
- `bam/`: Sorted BAM files, index `.bai` files, and `samtools flagstat` outputs.
- `pipeline_info/`: Nextflow execution reports (timeline, resource usage trace, DAG).
- Additional folders will contain the VCFs, coverage tables, and generated reports.
