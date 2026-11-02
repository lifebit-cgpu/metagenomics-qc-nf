process QC_SUMMARY {
    publishDir "${params.outdir}/qc_summary", mode: 'copy'

    input:
    path script
    path samples
    path fastp_json,     stageAs: 'fastp/*'
    path host_stats,     stageAs: 'host/*'
    path kraken_reports, stageAs: 'kraken/*'

    output:
    path 'qc_summary.csv', emit: csv
    path 'top_taxa.csv',   emit: taxa
    path '*_mqc.tsv',      emit: mqc

    script:
    """
    python3 ${script} \\
        --samples ${samples} \\
        --fastp_dir fastp --host_dir host --kraken_dir kraken \\
        --min_reads ${params.min_reads} \\
        --max_host_fraction ${params.max_host_fraction} \\
        --min_retention ${params.min_retention} \\
        --max_contaminant_fraction ${params.max_contaminant_fraction} \\
        --contaminant_taxids '${params.contaminant_taxids}' \\
        --top_n ${params.top_n}
    """
}
