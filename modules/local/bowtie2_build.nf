process BOWTIE2_BUILD {
    tag "${fasta.name}"

    input:
    path fasta

    output:
    path 'bt2_index', emit: index

    script:
    """
    mkdir -p bt2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bt2_index/host > bowtie2_build.log 2>&1
    """
}
