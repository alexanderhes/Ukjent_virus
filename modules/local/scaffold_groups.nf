/*
 * SCAFFOLD_GROUPS
 *
 * Groups each sample's BLAST-assigned contigs for reference-guided scaffolding:
 * one group per taxon, or per segment for segmented taxa, each with the
 * reference accession its contigs will be ordered against (highest total
 * bitscore in the group). Logic lives in bin/prepare_scaffold_groups.R.
 *
 * Artefact contigs are left out: only contigs whose BLAST hits cover
 * >= params.verdict_artefact_max_aln_pct % of the contig and
 * >= params.verdict_artefact_max_aln_bp bp are scaffolded (the overview's
 * artefact_suspected thresholds). This keeps long rRNA contigs with a short
 * match at the end of a viral reference out of the scaffolders.
 *
 * Outputs:
 *   {id}_scaffold_groups.tsv  : sample safe_taxon segment group_id ref_accession
 *                               ref_length n_contigs contig_ids
 *                               (header only when no contig passes)
 *   {id}_scaffold_contigs.tsv : per contig: sample safe_taxon contig contig_len
 *                               aln_bp aln_pct scaffolded (published with the
 *                               genomes by POLISH_SCAFFOLD)
 */

process SCAFFOLD_GROUPS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(blast_tsv), path(ref_lengths)

    output:
    tuple val(meta), path("${meta.id}_scaffold_groups.tsv"),  emit: groups
    tuple val(meta), path("${meta.id}_scaffold_contigs.tsv"), emit: contigs

    script:
    """
    prepare_scaffold_groups.R ${meta.id} ${blast_tsv} ${ref_lengths} \\
        ${meta.id}_scaffold_groups.tsv ${meta.id}_scaffold_contigs.tsv \\
        ${params.verdict_artefact_max_aln_pct} ${params.verdict_artefact_max_aln_bp}
    """
}
