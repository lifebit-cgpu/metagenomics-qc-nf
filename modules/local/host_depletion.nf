process HOST_DEPLETION {
    tag "${meta.id}"
    publishDir "${params.outdir}/host_depletion", mode: 'copy', pattern: '*.{log,tsv}'

    input:
    tuple val(meta), path(reads, stageAs: 'in?/*')
    path index

    output:
    tuple val(meta), path('*_nonhost_R{1,2}.fastq.gz'), emit: reads
    tuple val(meta), path('*.bowtie2.log'),            emit: log
    tuple val(meta), path('*.host_stats.tsv'),         emit: stats

    script:
    """
    IDX=\$(find -L ${index} -name '*.rev.1.bt2*' | head -n 1 | sed -E 's/\\.rev\\.1\\.bt2l?\$//')
    if [ -z "\$IDX" ]; then echo "No bowtie2 index found in ${index}" >&2; exit 1; fi

    bowtie2 -p ${task.cpus} -x "\$IDX" \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --un-conc-gz ${meta.id}_nonhost_R%.fastq.gz \\
        2> ${meta.id}.bowtie2.log > /dev/null

    touch ${meta.id}_nonhost_R1.fastq.gz ${meta.id}_nonhost_R2.fastq.gz
    IN=\$(( \$(gzip -dc ${reads[0]} | wc -l) / 4 ))
    OUT=\$(( \$(gzip -dc ${meta.id}_nonhost_R1.fastq.gz | wc -l) / 4 ))
    printf 'sample\\tinput_pairs\\tnonhost_pairs\\n%s\\t%s\\t%s\\n' ${meta.id} \$IN \$OUT > ${meta.id}.host_stats.tsv
    """
}
