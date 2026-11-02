process FASTQC {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: 'in?/*')

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip'),  emit: zip

    script:
    """
    ln -sf ${reads[0]} ${meta.id}_R1.fastq.gz
    ln -sf ${reads[1]} ${meta.id}_R2.fastq.gz
    fastqc --threads ${task.cpus} --quiet ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz
    """
}
