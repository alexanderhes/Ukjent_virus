/*
 * VISUALIZE_VALIDATION
 *
 * Reads the per-sample validation BLAST summary TSV and produces two
 * multi-page PDFs, one page per detected viral species, showing each assembled
 * contig as a horizontal bar spanning its aligned region on the reference
 * genome (longest alignment first):
 *   {id}_validation_contigs.pdf   : contigs that pass the artefact thresholds,
 *                                   i.e. the contigs SCAFFOLD_GROUPS scaffolds
 *   {id}_validation_artefacts.pdf : likely artefact contigs only (BLAST hits
 *                                   cover < params.verdict_artefact_max_aln_pct %
 *                                   or < params.verdict_artefact_max_aln_bp bp
 *                                   of the contig), labelled with their aligned share
 */

process VISUALIZE_VALIDATION {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    tuple val(meta), path(summary_tsv)
    path(ref_lengths_files)     // all per-species ref length TSVs (staged, used to set genome x-axis limits)

    output:
    tuple val(meta), path("${meta.id}_validation_contigs.pdf"),   emit: pdf
    tuple val(meta), path("${meta.id}_validation_artefacts.pdf"), emit: artefacts_pdf

    script:
    """
    visualize_validation.R contigs \\
        ${summary_tsv} \\
        . \\
        ${meta.id}_validation_contigs.pdf \\
        ${meta.id}_validation_artefacts.pdf \\
        ${params.verdict_artefact_max_aln_pct} \\
        ${params.verdict_artefact_max_aln_bp}
    """
}
