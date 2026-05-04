#!/usr/bin/env nextflow

/*
 * Report Generation Module
 * Creates final summary reports
 */

process GENERATE_REPORTS {
    tag "generate_reports"
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    input:
    path snp_results
    path indel_results
    path coverage_ch
    val outdir
    
    output:
    path "summary_report.html"
    path "variant_stats.txt"
    
    script:
    """
    # Create summary statistics
    cat > summary_stats.txt << EOL
=== Variant Calling Summary ===
SNP Count: \$(zcat ${snp_results[0]} 2>/dev/null | grep -v "^#" | wc -l)
Indel Count: \$(zcat ${indel_results[0]} 2>/dev/null | grep -v "^#" | wc -l)
EOL
    
    # Generate HTML report
    cat > summary_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Amplicon Variant Calling Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        .section { margin-bottom: 30px; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
    </style>
</head>
<body>
    <h1>Amplicon Nanopore Variant Calling Pipeline Report</h1>
    <div class="section">
        <h2>Pipeline Parameters</h2>
        <ul>
            <li>Reference: ${params.ref}</li>
            <li>Target BED: ${params.bed}</li>
            <li>Min Base Quality: ${params.min_base_qual}</li>
            <li>SnpEff Database: ${params.snpeff_db}</li>
        </ul>
    </div>
    <div class="section">
        <h2>Variant Statistics</h2>
        <pre>$(cat summary_stats.txt)</pre>
    </div>
</body>
</html>
EOF
    
    # Create final variant summary
    echo "Sample\tSNP_Count\tIndel_Count" > variant_stats.txt
    """
}