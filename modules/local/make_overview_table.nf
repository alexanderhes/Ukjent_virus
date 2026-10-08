/*
 * MAKE_OVERVIEW_TABLE
 *
 * Produces a comprehensive per-sample and batch-level overview TSV combining
 * three sections for each detected virus (one row per assembly-set):
 *
 *   Part 1 — Read funnel:
 *     raw_reads → host_filtered_reads → trimmed_reads → dedup_reads,
 *     with percentage removal at each step.
 *
 *   Part 2 — EsViritu results:
 *     Breadth of coverage (%), RPM, RPKMR, RPKMF, ANI, Pi, read count.
 *     Correctly handles segmented viruses (one row per assembly-set).
 *
 *   Part 3 — Assembly + BLAST validation (NA when --validate not used):
 *     assembly_status, n_contigs, longest_contig_bp,
 *     blast_coverage (single % for non-segmented, per-segment string for segmented),
 *     blast_identity_pct.
 *
 *   Part 4 — Verdict (verdict_* params):
 *     esv_verdict (mapping evidence), blast_verdict (de novo + BLAST evidence,
 *     from contig_aln_bp / contig_aln_pct), combined verdict, and flags
 *     (general, esv_flags, blast_flags).
 *
 * The R logic lives in bin/make_overview_table.R; scalar params are passed
 * as command-line arguments using Rscript's trailingOnly mechanism.
 *
 * Outputs:
 *   esv_staged.overview.tsv  — all samples combined
 *   {sample_id}_overview.tsv — one file per sample (also in overview/ dir)
 */

process MAKE_OVERVIEW_TABLE {
    label 'process_low'

    publishDir "${params.outdir}/esviritu_batch/esv_summary", mode: 'copy',
        pattern: "esv_staged.overview.tsv"
    publishDir "${params.outdir}/overview", mode: 'copy',
        pattern: "esv_staged.overview.tsv",
        saveAs: { "all_samples_overview.tsv" }
    publishDir "${params.outdir}/overview", mode: 'copy',
        pattern: "*_overview.tsv"

    input:
    path(read_stats_tsv)
    path(assembly_summary_tsv)
    path(info_enriched_tsv)
    path(validation_tsvs)       // collected validation summary TSVs, or a dummy file when --validate is off
    path(has_contigs_files)     // per-species has_contigs.txt markers from BLASTN_VALIDATE, or dummy
    path(ref_lengths_files)     // per-species reference sequence length TSVs from BLASTN_VALIDATE, or dummy
    path(db_metadata_tsv)       // staged EsViritu DB metadata TSV for BLAST-only taxonomy enrichment
    val(run_validate)
    val(validate_min_reads)
    val(assembly_taxon_level)

    output:
    path "esv_staged.overview.tsv",  emit: batch_overview
    path "*_overview.tsv",           emit: sample_overviews

    script:
    """
    make_overview_table.R ${run_validate} ${validate_min_reads} ${assembly_taxon_level} \
        ${params.verdict_esv_min_reads} ${params.verdict_esv_min_breadth_pct} \
        ${params.verdict_contig_min_aln_bp} ${params.verdict_contig_min_aln_pct} \
        ${params.verdict_artefact_max_aln_bp} ${params.verdict_artefact_max_aln_pct} \
        ${params.verdict_esv_divergent_identity_pct} ${params.verdict_blast_divergent_identity_pct} \
        ${params.verdict_recurrent_min_samples}
    """
}

