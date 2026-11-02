process KRAKEN2 {
    tag "${meta.id}"
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: 'in?/*')
    path db

    output:
    tuple val(meta), path('*.kraken2.report.txt'), emit: report
    tuple val(meta), path('*.kraken2.log'),        emit: log

    script:
    """
    kraken2 --db ${db} --threads ${task.cpus} --paired --gzip-compressed \\
        --report ${meta.id}.kraken2.report.txt --output /dev/null \\
        ${reads[0]} ${reads[1]} 2> ${meta.id}.kraken2.log
    """
}
