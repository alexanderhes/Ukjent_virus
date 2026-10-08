/*
 * EXTRACT_HOST_FASTA + MAKE_HOST_BLASTDB
 *
 * Build the host BLAST database used by FILTER_HOST_CONTIGS when it does not
 * exist yet (main.nf checks params.host_blastdb at start-up). The host
 * sequences are reconstructed from the bowtie2 host index with bowtie2-inspect,
 * so the BLAST database always matches the read-level host reference and no
 * separate FASTA has to be configured.
 *
 * MAKE_HOST_BLASTDB writes into the directory of params.host_blastdb via
 * storeDir, so the database is built once and reused by later runs. The user
 * running the pipeline needs write access to that directory.
 */

process EXTRACT_HOST_FASTA {
    tag "${index}"
    label 'process_medium'

    input:
    path(bt2_index)   // all *.bt2 host index files staged into work dir

    output:
    path("host_from_index.fa.gz"), emit: fasta

    script:
    index = file(params.host_index).name
    """
    bowtie2-inspect ${index} | gzip -1 > host_from_index.fa.gz
    """
}

process MAKE_HOST_BLASTDB {
    tag "${db_name}"
    label 'process_medium'

    storeDir "${file(params.host_blastdb).parent}"

    input:
    path(fasta_gz)

    output:
    path("${db_name}.*"), emit: db

    script:
    db_name = file(params.host_blastdb).name
    """
    gzip -dc ${fasta_gz} \\
        | makeblastdb -in - -dbtype nucl -out ${db_name} \\
            -title "Host: ${file(params.host_index).name} (bowtie2-inspect of the host index)"
    """
}
