process MAKE_TEST_DATA {
    tag "seed ${seed}"
    publishDir "${params.outdir}/test_data", mode: 'copy'

    input:
    path script
    val seed

    output:
    path 'samplesheet.csv',     emit: samplesheet
    path 'reads/*.fastq.gz',    emit: reads
    path 'host_genome.fa',      emit: host_fasta
    path 'kraken_library.fa',   emit: kraken_library
    path 'taxonomy',            emit: taxonomy
    path 'ground_truth*.csv',   emit: truth

    script:
    """
    python3 ${script} --outdir . --seed ${seed}
    """
}
