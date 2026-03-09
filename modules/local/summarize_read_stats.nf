/*
 * SUMMARIZE_READ_STATS
 *
 * Runs after all COLLECT_READ_STATS processes and SUMMARIZE_ESV have
 * finished. Produces two outputs in esviritu_batch/esv_summary/:
 *
 *   esv_staged.read_stats.tsv
 *     Batch-level read-count funnel: one row per sample with
 *     raw_reads, host_filtered_reads, trimmed_reads, dedup_reads,
 *     and the three derived removal percentages.
 *
 *   esv_staged.detected_virus.info.enriched.tsv
 *     The EsViritu per-accession info table joined with raw_reads.
 *     RPM/RPKMR are computed downstream in make_overview_table.R at
 *     assembly-set level (correctly handling segmented virus lengths).
 */

process SUMMARIZE_READ_STATS {
    label 'process_low'

    publishDir "${params.outdir}/esviritu_batch/esv_summary", mode: 'copy'

    input:
    path(read_stats_files)   // collected list of per-sample _read_stats.tsv
    path(info_tsv)           // esv_staged.detected_virus.info.tsv from SUMMARIZE_ESV

    output:
    path "esv_staged.read_stats.tsv",                        emit: read_stats
    path "esv_staged.detected_virus.info.enriched.tsv",      emit: info_enriched

    script:
    """
    Rscript summarize_read_stats.R
    """
}
