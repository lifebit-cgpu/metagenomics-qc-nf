process FASTP {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.{json,html,log}'

    input:
    tuple val(meta), path(reads, stageAs: 'in?/*')

    output:
    tuple val(meta), path('*_trimmed_R{1,2}.fastq.gz'), emit: reads
    tuple val(meta), path('*.fastp.json'),             emit: json
    tuple val(meta), path('*.fastp.html'),             emit: html
    tuple val(meta), path('*.fastp.log'),              emit: log

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${meta.id}_trimmed_R1.fastq.gz -O ${meta.id}_trimmed_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
        --length_required ${params.fastp_min_length} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
        2> ${meta.id}.fastp.log
    """
}
