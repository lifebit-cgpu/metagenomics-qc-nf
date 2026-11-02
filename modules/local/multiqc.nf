process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path files, stageAs: '?/*'
    path config

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc -f -c ${config} -n multiqc_report.html .
    """
}
