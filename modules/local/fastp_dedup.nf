/*
 * FASTP_DEDUP
 *
 * Pass 2 — deduplication only (no re-trimming).
 * Runs on already-trimmed reads so that adapter/quality differences
 * no longer mask identical sequences. This reveals the full duplicate
 * population (~40-50%) compared to ~5% visible in raw reads.
 */

process FASTP_DEDUP {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/fastp_dedup/${meta.id}", mode: 'copy',
        pattern: "*.{json,html,log}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_dedup_R1.fastq.gz"),
                     path("${meta.id}_dedup_R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_fastp_dedup.json"),  emit: json
    path "${meta.id}_fastp_dedup.html",                    emit: html
    path "${meta.id}_fastp_dedup.log",                     emit: log

    script:
    """
    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${meta.id}_dedup_R1.fastq.gz \\
        -O ${meta.id}_dedup_R2.fastq.gz \\
        --dedup \\
        --thread ${task.cpus} \\
        --json ${meta.id}_fastp_dedup.json \\
        --html ${meta.id}_fastp_dedup.html \\
        2> ${meta.id}_fastp_dedup.log
    """
}
