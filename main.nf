#!/usr/bin/env nextflow
/*
 * metagenomics-qc-nf
 * QC for vendor-delivered (e.g. CosmosID) paired-end shotgun metagenomics FASTQ.
 *   FastQC (raw) -> fastp trimming -> host-read detection and depletion (only when metagenomic == TRUE)
 *   -> classified fraction and top taxa (Kraken2 database, or a bowtie2 reference panel)
 *   -> per-sample QC table with PASS/WARN/FAIL flags -> MultiQC report
 */
nextflow.enable.dsl = 2

include { MAKE_TEST_DATA; FASTQC; FASTP; HOST_DEPLETE; CLASSIFY_PANEL; KRAKEN2; QC_SUMMARY; MULTIQC } from './modules.nf'
include { BOWTIE2_BUILD as BUILD_HOST_INDEX; BOWTIE2_BUILD as BUILD_PANEL_INDEX } from './modules.nf'

def isTrue(v) {
    return v != null && v.toString().trim().toUpperCase() in ['TRUE', 'T', '1', 'YES', 'Y']
}

workflow {
    def input_desc = params.generate_test_data ? "synthetic test data (seed ${params.test_seed})" : "${params.input}"
    def host_desc  = params.host_index ?: (params.host_fasta ?: (params.generate_test_data ? 'synthetic host genome' : 'none'))
    def class_desc = params.kraken2_db ? "kraken2 db ${params.kraken2_db}" : (params.classifier_fasta ? "bowtie2 panel ${params.classifier_fasta}" : (params.generate_test_data ? 'bowtie2 panel (synthetic)' : 'none'))
    log.info "metagenomics-qc-nf ${workflow.manifest.version} | input: ${input_desc} | host: ${host_desc} | classifier: ${class_desc} | outdir: ${params.outdir}"

    qc_script = file("${projectDir}/bin/qc_summary.py", checkIfExists: true)

    if (params.generate_test_data) {
        MAKE_TEST_DATA(file("${projectDir}/bin/make_test_data.py", checkIfExists: true))
        ch_files = MAKE_TEST_DATA.out.reads.flatten().map { f -> [f.name, f] }
        ch_samples = MAKE_TEST_DATA.out.samplesheet
            .splitCsv(header: true)
            .map { r -> [r.fastq_1, r] }
            .join(ch_files)
            .map { n1, r, f1 -> [r.fastq_2, r, f1] }
            .join(ch_files)
            .map { n2, r, f1, f2 -> [[id: r.sample, metagenomic: isTrue(r.metagenomic)], [f1, f2]] }
        ch_host_fasta  = MAKE_TEST_DATA.out.host
        ch_panel_fasta = MAKE_TEST_DATA.out.panel
        ch_truth       = MAKE_TEST_DATA.out.truth
    } else {
        if (!params.input) {
            error "Provide --input samplesheet.csv (columns: sample,fastq_1,fastq_2,metagenomic) or run with -profile test"
        }
        ch_samples = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { r -> [[id: r.sample, metagenomic: isTrue(r.metagenomic)], [file(r.fastq_1, checkIfExists: true), file(r.fastq_2, checkIfExists: true)]] }
        ch_host_fasta  = params.host_fasta ? Channel.value(file(params.host_fasta, checkIfExists: true)) : Channel.empty()
        ch_panel_fasta = params.classifier_fasta ? Channel.value(file(params.classifier_fasta, checkIfExists: true)) : Channel.empty()
        ch_truth       = Channel.value(file("${projectDir}/assets/NO_FILE", checkIfExists: true))
    }

    FASTQC(ch_samples)
    FASTP(ch_samples)

    // Host-read detection and depletion, only for samples with metagenomic == TRUE
    if (params.host_index) {
        ch_host_index = Channel.value(file(params.host_index, type: 'dir', checkIfExists: true))
    } else if (params.host_fasta || params.generate_test_data) {
        BUILD_HOST_INDEX(ch_host_fasta.map { f -> ['host', f] })
        ch_host_index = BUILD_HOST_INDEX.out.index
    } else {
        error "Provide --host_index (directory holding a bowtie2 index) or --host_fasta for host-read detection"
    }

    ch_trimmed = FASTP.out.reads.branch { meta, reads ->
        host: meta.metagenomic
        skip: true
    }
    HOST_DEPLETE(ch_trimmed.host, ch_host_index)
    ch_clean = HOST_DEPLETE.out.reads.mix(ch_trimmed.skip)

    // Classified fraction and top taxa
    ch_class = Channel.empty()
    if (params.kraken2_db) {
        KRAKEN2(ch_clean, Channel.value(file(params.kraken2_db, type: 'dir', checkIfExists: true)))
        ch_class = KRAKEN2.out.report
    } else if (params.classifier_fasta || params.generate_test_data) {
        BUILD_PANEL_INDEX(ch_panel_fasta.map { f -> ['panel', f] })
        CLASSIFY_PANEL(ch_clean, BUILD_PANEL_INDEX.out.index)
        ch_class = CLASSIFY_PANEL.out.counts.mix(CLASSIFY_PANEL.out.log)
    } else {
        log.warn "No --kraken2_db or --classifier_fasta given: classified fraction and top taxa will be NA"
    }

    ch_sample_list = ch_samples
        .map { meta, reads -> "${meta.id},${meta.metagenomic}".toString() }
        .collectFile(name: 'samples.csv', newLine: true, sort: true)

    QC_SUMMARY(
        qc_script,
        ch_sample_list,
        FASTP.out.json.collect(),
        HOST_DEPLETE.out.counts.collect().ifEmpty([]),
        ch_class.collect().ifEmpty([]),
        ch_truth
    )

    MULTIQC(
        FASTQC.out.zip
            .mix(FASTP.out.json, HOST_DEPLETE.out.log, QC_SUMMARY.out.mqc)
            .collect()
    )
}

workflow.onComplete {
    log.info "metagenomics-qc-nf ${workflow.success ? 'completed successfully' : 'FAILED'}; results in ${params.outdir}"
}
