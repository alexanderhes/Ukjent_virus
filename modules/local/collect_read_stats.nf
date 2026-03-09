/*
 * COLLECT_READ_STATS
 *
 * Parses the raw-read count (from HOST_FILTER) and the two fastp JSON
 * files (trim + dedup) to produce a single-row TSV per sample with:
 *
 *   sample_ID | raw_reads | host_filtered_reads | trimmed_reads | dedup_reads
 *   host_removal_pct | trim_removed_pct | dup_rate_pct
 *
 * raw_reads            = R1 count × 2 (pre-host-filter, true sequencer output)
 * host_filtered_reads  = paired reads entering FASTP_TRIM (post-re-pairing)
 * trimmed_reads        = reads passing FASTP_TRIM quality/complexity filters
 * dedup_reads          = reads after FASTP_DEDUP (= EsViritu denominator)
 *
 * Logic lives in bin/collect_read_stats.py.
 */

process COLLECT_READ_STATS {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/read_stats/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(raw_reads_txt), path(fastp_trim_json), path(fastp_dedup_json)

    output:
    tuple val(meta), path("${meta.id}_read_stats.tsv"), emit: tsv

    script:
    """
    collect_read_stats.py \\
        ${meta.id} \\
        ${raw_reads_txt} \\
        ${fastp_trim_json} \\
        ${fastp_dedup_json} \\
        ${meta.id}_read_stats.tsv
    """
}
