/*
 * VISUALIZE_VALIDATION
 *
 * Reads the per-sample validation BLAST summary TSV and produces a multi-page
 * PDF: one page per detected viral species, showing each assembled contig as a
 * horizontal bar spanning its aligned region on the reference genome.
 * Contigs are ordered by alignment length (longest first).
 */

process VISUALIZE_VALIDATION {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    tuple val(meta), path(summary_tsv)
    path(ref_lengths_files)     // all per-species ref length TSVs (staged, used to set genome x-axis limits)

    output:
    tuple val(meta), path("${meta.id}_validation_contigs.pdf"), emit: pdf

    script:
    """
    visualize_validation.R \\
        ${summary_tsv} \\
        ${meta.id}_validation_contigs.pdf \\
        .
    """
}
