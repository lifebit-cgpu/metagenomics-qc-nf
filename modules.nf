// Processes for metagenomics-qc-nf

process MAKE_TEST_DATA {
    tag 'synthetic'
    label 'process_single'
    publishDir "${params.outdir}/test_data", mode: params.publish_mode

    input:
    path script

    output:
    path 'samplesheet.csv',    emit: samplesheet
    path 'reads/*.fastq.gz',   emit: reads
    path 'host_synthetic.fa',  emit: host
    path 'panel_synthetic.fa', emit: panel
    path 'ground_truth.csv',   emit: truth

    script:
    """
    python3 ${script} --outdir . --seed ${params.test_seed} \\
        --min-raw-pairs ${params.min_raw_pairs} \\
        --warn-host-frac ${params.warn_host_frac} \\
        --warn-classified-frac ${params.warn_classified_frac}
    """
}

process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc_raw", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.zip',  emit: zip
    path '*_fastqc.html', emit: html

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """
}

process FASTP {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/fastp", mode: params.publish_mode, pattern: '*.{json,html,log}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_{1,2}.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html
    path "${meta.id}.fastp.log",  emit: log

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${meta.id}_trimmed_1.fastq.gz -O ${meta.id}_trimmed_2.fastq.gz \\
        --detect_adapter_for_pe \\
        --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \\
        --qualified_quality_phred 15 \\
        --length_required ${params.min_length} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
        2> ${meta.id}.fastp.log
    """
}

process BOWTIE2_BUILD {
    tag "${name}"
    label 'process_medium'

    input:
    tuple val(name), path(fasta)

    output:
    path "${name}_bowtie2_index", emit: index

    script:
    """
    mkdir ${name}_bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} ${name}_bowtie2_index/${name} > ${name}.bowtie2_build.log
    """
}

process HOST_DEPLETE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/host_depletion", mode: params.publish_mode, pattern: '*.{log,tsv}'
    publishDir "${params.outdir}/host_depleted_reads", mode: params.publish_mode, pattern: '*.fastq.gz', enabled: params.save_host_depleted

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}_nohost_{1,2}.fastq.gz"), emit: reads
    path "${meta.id}.host.bowtie2.log", emit: log
    path "${meta.id}.host_counts.tsv",  emit: counts

    script:
    // Index prefix = the single *.rev.1.bt2[l] file in the index directory, minus its suffix
    """
    INDEX=\$(ls ${index}/*.rev.1.bt2* | sed -E 's/\\.rev\\.1\\.bt2l?\$//')
    bowtie2 -p ${task.cpus} --sensitive -x \$INDEX \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        --un-conc-gz ${meta.id}_nohost_%.fastq.gz \\
        -S /dev/null 2> ${meta.id}.host.bowtie2.log
    INPUT=\$(grep -m1 'reads; of these' ${meta.id}.host.bowtie2.log | awk '{print \$1}')
    if [ -s ${meta.id}_nohost_1.fastq.gz ]; then
        NONHOST=\$(( \$(zcat ${meta.id}_nohost_1.fastq.gz | wc -l) / 4 ))
    else
        NONHOST=0
        touch ${meta.id}_nohost_1.fastq.gz ${meta.id}_nohost_2.fastq.gz
    fi
    printf 'sample\\tinput_pairs\\tnonhost_pairs\\n%s\\t%s\\t%s\\n' ${meta.id} \$INPUT \$NONHOST > ${meta.id}.host_counts.tsv
    """
}

process CLASSIFY_PANEL {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/classification", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    path index

    output:
    path "${meta.id}.panel_counts.tsv",  emit: counts
    path "${meta.id}.panel.bowtie2.log", emit: log

    script:
    // Read 1 of each (host-depleted) pair is aligned to the reference panel; FASTA headers are "<taxon>|<contig>"
    """
    INDEX=\$(ls ${index}/*.rev.1.bt2* | sed -E 's/\\.rev\\.1\\.bt2l?\$//')
    set -o pipefail
    bowtie2 -p ${task.cpus} --sensitive -x \$INDEX -U ${reads[0]} \\
        --no-unal --no-hd 2> ${meta.id}.panel.bowtie2.log \\
      | cut -f 3 | sort | uniq -c \\
      | awk 'BEGIN{OFS="\\t"; print "reference","reads"} {print \$2,\$1}' > ${meta.id}.panel_counts.tsv
    """
}

process KRAKEN2 {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/classification", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    path db

    output:
    path "${meta.id}.kraken2.report.txt", emit: report

    script:
    """
    kraken2 --db ${db} --threads ${task.cpus} --paired \\
        --report ${meta.id}.kraken2.report.txt --output /dev/null \\
        ${reads[0]} ${reads[1]}
    """

    stub:
    """
    printf ' 20.00\\t200\\t200\\tU\\t0\\tunclassified\\n 80.00\\t800\\t0\\tR\\t1\\troot\\n 50.00\\t500\\t500\\tS\\t816\\t    Bacteroides stub\\n 30.00\\t300\\t300\\tS\\t853\\t    Faecalibacterium stub\\n' > ${meta.id}.kraken2.report.txt
    """
}

process QC_SUMMARY {
    label 'process_single'
    publishDir "${params.outdir}/qc_summary", mode: params.publish_mode

    input:
    path script
    path samples
    path fastp_json
    path host_counts
    path classification
    path truth

    output:
    path 'qc_summary.csv',        emit: csv
    path 'top_taxa.csv',          emit: taxa
    path 'qc_summary_mqc.tsv',    emit: mqc
    path 'validation_report.csv', emit: validation, optional: true

    script:
    """
    python3 ${script} --samples ${samples} \\
        --min-raw-pairs ${params.min_raw_pairs} \\
        --min-retained-frac ${params.min_retained_frac} \\
        --warn-host-frac ${params.warn_host_frac} \\
        --warn-classified-frac ${params.warn_classified_frac} \\
        --top-n ${params.top_n_taxa}
    """
}

process MULTIQC {
    label 'process_single'
    publishDir "${params.outdir}/multiqc", mode: params.publish_mode

    input:
    path '*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc --force --title 'Metagenomics QC' --filename multiqc_report.html .
    """
}
