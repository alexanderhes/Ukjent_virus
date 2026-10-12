/*
 * VISUALIZE_GENOMES
 *
 * Plots the polished genomes of one sample (POLISH_SCAFFOLD) on the references
 * they were built against (logic in bin/visualize_validation.R, genomes mode):
 *   {id}_validation_genomes.pdf : one page per species, laid out like the
 *                                 contig plot; full-length genomes with bases
 *                                 coloured by source (contig base kept,
 *                                 read-called inside a contig, read-called
 *                                 outside the contigs, N) at their reference
 *                                 positions;
 *                                 contig overhangs beyond the reference ends
 *   {id}_genome_coverage.pdf    : one page per genome (per segment): read
 *                                 depth (dashed line at params.polish_min_depth),
 *                                 consensus by source with N in black, and the
 *                                 contigs placed by SCAFFOLD_MINIMAP2
 */

process VISUALIZE_GENOMES {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    tuple val(meta), path(genome_stats), path(genome_positions), path(contig_aln)
    path(ref_lengths_files)     // all per-species ref length TSVs (segment labels, x-axis limits)

    output:
    tuple val(meta), path("${meta.id}_validation_genomes.pdf"), emit: pdf
    tuple val(meta), path("${meta.id}_genome_coverage.pdf"),    emit: coverage_pdf

    script:
    """
    visualize_validation.R genomes \\
        ${genome_stats} \\
        ${genome_positions} \\
        ${contig_aln} \\
        . \\
        ${meta.id}_validation_genomes.pdf \\
        ${meta.id}_genome_coverage.pdf \\
        ${params.polish_min_depth}
    """
}
