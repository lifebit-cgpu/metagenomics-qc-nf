process KRAKEN2_BUILD {
    tag "${library.name}"

    input:
    path library
    path taxonomy

    output:
    path 'kraken2_db', emit: db

    script:
    """
    mkdir -p kraken2_db/taxonomy
    cp ${taxonomy}/nodes.dmp ${taxonomy}/names.dmp kraken2_db/taxonomy/
    kraken2-build --add-to-library ${library} --db kraken2_db --no-masking
    kraken2-build --build --db kraken2_db --threads ${task.cpus}
    ls -l kraken2_db
    """
}
