#!/usr/bin/env nextflow
/*
 * metagenomics-qc-nf
 * QC for vendor-delivered (e.g. CosmosID) paired-end shotgun metagenomics FASTQ.
 *   FastQC (raw) -> fastp (trimming) -> bowtie2 host depletion (only metagenomic=TRUE)
 *   -> Kraken2 (classified fraction, top taxa, contamination) -> QC summary -> MultiQC
 */
nextflow.enable.dsl = 2

include { MAKE_TEST_DATA } from './modules/local/make_test_data'
include { FASTQC         } from './modules/local/fastqc'
include { FASTP          } from './modules/local/fastp'
include { BOWTIE2_BUILD  } from './modules/local/bowtie2_build'
include { HOST_DEPLETION } from './modules/local/host_depletion'
include { KRAKEN2_BUILD  } from './modules/local/kraken2_build'
include { KRAKEN2        } from './modules/local/kraken2'
include { QC_SUMMARY     } from './modules/local/qc_summary'
include { MULTIQC        } from './modules/local/multiqc'

def toBool(v) {
    return v != null && v.toString().trim().toUpperCase() in ['TRUE', 'T', 'YES', 'Y', '1']
}

def resolvePath(p, base) {
    def s = p.toString().trim()
    if (s ==~ /^(s3|gs|az|https?|ftp):\/\/.*/ || s.startsWith('/')) {
        return file(s, checkIfExists: true)
    }
    return file("${base}/${s}", checkIfExists: true)
}

workflow {
    log.info """
    metagenomics-qc-nf  v${workflow.manifest.version}
    ==========================================
    test_data   : ${params.test_data}
    input       : ${params.input}
    host_fasta  : ${params.host_fasta}
    bowtie2_idx : ${params.bowtie2_index}
    kraken2_db  : ${params.kraken2_db}
    outdir      : ${params.outdir}
    """.stripIndent()

    gen_script = file("${projectDir}/scripts/make_synthetic_data.py", checkIfExists: true)
    qc_script  = file("${projectDir}/scripts/qc_summary.py", checkIfExists: true)
    mqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)

    def hostOn   = params.test_data || params.host_fasta || params.bowtie2_index
    def krakenOn = !params.skip_kraken2 && (params.test_data || params.kraken2_db)

    ch_host_fasta  = Channel.empty()
    ch_kraken_lib  = Channel.empty()
    ch_kraken_tax  = Channel.empty()

    if (params.test_data) {
        // ---- synthetic, deterministic test data generated inside the pipeline ----
        MAKE_TEST_DATA(gen_script, params.test_seed)
        ch_meta  = MAKE_TEST_DATA.out.samplesheet
                        .splitCsv(header: true)
                        .map { r -> [ r.sample, toBool(r.metagenomic) ] }
        ch_files = MAKE_TEST_DATA.out.reads
                        .flatten()
                        .map { f -> [ f.name.replaceAll(/_R[12]\.fastq\.gz$/, ''), f ] }
                        .groupTuple(sort: true)
        ch_reads = ch_meta.join(ch_files)
                        .map { s, mg, fq -> [ [ id: s, metagenomic: mg ], fq ] }
        ch_host_fasta = MAKE_TEST_DATA.out.host_fasta
        ch_kraken_lib = MAKE_TEST_DATA.out.kraken_library
        ch_kraken_tax = MAKE_TEST_DATA.out.taxonomy
    } else {
        if (!params.input) {
            error "Provide --input <samplesheet.csv> (columns: sample,fastq_1,fastq_2,metagenomic) or run with -profile test"
        }
        def base = file(params.input).parent
        ch_reads = Channel.fromPath(params.input, checkIfExists: true)
                        .splitCsv(header: true)
                        .map { r ->
                            if (!r.sample || !r.fastq_1 || !r.fastq_2) {
                                error "Samplesheet row is missing sample/fastq_1/fastq_2: ${r}"
                            }
                            [ [ id: r.sample.trim(), metagenomic: toBool(r.metagenomic) ],
                              [ resolvePath(r.fastq_1, base), resolvePath(r.fastq_2, base) ] ]
                        }
        if (params.host_fasta) {
            ch_host_fasta = Channel.value(file(params.host_fasta, checkIfExists: true))
        }
        if (!hostOn) {
            log.warn "No --host_fasta or --bowtie2_index given: host depletion is skipped for ALL samples"
        }
        if (!krakenOn) {
            log.warn "No --kraken2_db given (or --skip_kraken2): classified fraction and top taxa are not computed"
        }
    }

    // ---------------- QC of raw reads + trimming ----------------
    FASTQC(ch_reads)
    FASTP(ch_reads)

    // ---------------- host index ----------------
    if (params.bowtie2_index) {
        ch_index = Channel.value(file(params.bowtie2_index, checkIfExists: true))
    } else {
        BOWTIE2_BUILD(ch_host_fasta)
        ch_index = BOWTIE2_BUILD.out.index.first()
    }

    // ---------------- host depletion controlled by the metagenomic flag ----------------
    ch_trim = FASTP.out.reads.branch { meta, r ->
        deplete: hostOn && meta.metagenomic
        skip: true
    }
    HOST_DEPLETION(ch_trim.deplete, ch_index)
    ch_clean = HOST_DEPLETION.out.reads.mix(ch_trim.skip)

    // ---------------- taxonomic classification ----------------
    ch_kraken_reports = Channel.empty()
    if (krakenOn) {
        if (params.kraken2_db) {
            ch_kraken_db = Channel.value(file(params.kraken2_db, checkIfExists: true))
        } else {
            KRAKEN2_BUILD(ch_kraken_lib, ch_kraken_tax)
            ch_kraken_db = KRAKEN2_BUILD.out.db.first()
        }
        KRAKEN2(ch_clean, ch_kraken_db)
        ch_kraken_reports = KRAKEN2.out.report.map { it[1] }
    }

    // ---------------- per-sample QC summary ----------------
    ch_samples = ch_reads
                    .map { m, r -> "${m.id},${m.metagenomic}" }
                    .collectFile(name: 'samples.csv', newLine: true, sort: true)

    QC_SUMMARY(
        qc_script,
        ch_samples,
        FASTP.out.json.map { it[1] }.collect(),
        HOST_DEPLETION.out.stats.map { it[1] }.collect().ifEmpty([]),
        ch_kraken_reports.collect().ifEmpty([])
    )

    // ---------------- MultiQC ----------------
    ch_mqc = FASTQC.out.zip.map { it[1] }
                .mix(FASTP.out.json.map { it[1] })
                .mix(HOST_DEPLETION.out.log.map { it[1] })
                .mix(ch_kraken_reports)
                .mix(QC_SUMMARY.out.mqc)
                .flatten()
                .collect()
    MULTIQC(ch_mqc, mqc_config)
}

workflow.onComplete {
    log.info "Pipeline finished: ${workflow.success ? 'SUCCESS' : 'FAILED'} | duration ${workflow.duration} | results in ${params.outdir}"
}
