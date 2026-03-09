/*
 * SUMMARIZE_VALIDATION
 *
 * Concatenates all per-species BLAST result TSVs for a single sample into one
 * summary file. The header is taken from the first file; subsequent files have
 * their header row stripped before appending.
 *
 * Output: {id}_validation_summary.tsv published to validation/{id}/
 */

process SUMMARIZE_VALIDATION {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/validation/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(blast_tsvs)

    output:
    tuple val(meta), path("${meta.id}_validation_summary.tsv"), emit: summary

    script:
    """
    # Write header from the first available TSV
    first_file=\$(ls ${meta.id}_*_blastn.tsv | head -1)
    head -1 "\${first_file}" > ${meta.id}_validation_summary.tsv

    # Append data rows from all TSVs (skip header line in each)
    for f in ${meta.id}_*_blastn.tsv; do
        tail -n +2 "\${f}" >> ${meta.id}_validation_summary.tsv
    done
    """
}
