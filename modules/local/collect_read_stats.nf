/*
 * COLLECT_READ_STATS
 *
 * Parses the two fastp JSON files (trim + dedup) and the HOST_FILTER read
 * count to produce a single-row TSV per sample with:
 *
 *   sample_ID | raw_reads | trimmed_reads | dedup_reads | host_filtered_reads
 *   trim_removed_pct | dup_rate_pct | host_removal_pct
 *
 * Pipeline order is FASTP_TRIM -> FASTP_DEDUP -> HOST_FILTER:
 * raw_reads            = reads entering FASTP_TRIM (R1 + R2, true sequencer output)
 * trimmed_reads        = reads passing FASTP_TRIM quality/complexity filters
 * dedup_reads          = reads after FASTP_DEDUP (entering HOST_FILTER)
 * host_filtered_reads  = reads left after HOST_FILTER (= EsViritu/SPAdes input)
 * Percentages are relative to each step's input (raw, trimmed, dedup reads).
 *
 * Logic lives in bin/collect_read_stats.py.
 */

process COLLECT_READ_STATS {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/read_stats/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(fastp_trim_json), path(fastp_dedup_json), path(hostfilt_count)

    output:
    tuple val(meta), path("${meta.id}_read_stats.tsv"), emit: tsv

    script:
    """
    collect_read_stats.py \\
        ${meta.id} \\
        ${fastp_trim_json} \\
        ${fastp_dedup_json} \\
        ${hostfilt_count} \\
        ${meta.id}_read_stats.tsv
    """
}
