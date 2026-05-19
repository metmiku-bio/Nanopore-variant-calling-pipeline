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