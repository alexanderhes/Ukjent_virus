/*
 * SPLIT_VIRAL_READS
 *
 * Extracts reads from the third-pass EsViritu BAM that map to any of the
 * reference accessions belonging to a given viral species. Multiple accessions
 * per species (different strains / segments) are pooled into one paired FASTQ.
 *
 * Input meta must carry:
 *   meta.id          — sample identifier
 *   meta.safe_species — filesystem-safe species name (spaces → underscores,
 *                       leading "s__" stripped)
 */

process SPLIT_VIRAL_READS {
    tag "${meta.id}:${meta.safe_species}"
    label 'process_medium'

    publishDir "${params.outdir}/validation/${meta.id}/reads", mode: 'copy'

    input:
    tuple val(meta), val(accessions), path(bam), path(bai)

    output:
    tuple val(meta),
          path("${meta.id}_${meta.safe_species}_R1.fastq.gz"),
          path("${meta.id}_${meta.safe_species}_R2.fastq.gz"), emit: reads

    script:
    // Build space-separated region list for samtools view
    def regions = accessions.join(' ')
    """
    # Extract reads mapping to any accession in this species group.
    # -F 2304 excludes secondary (0x100) and supplementary (0x800) alignments
    # so each read name appears exactly once in the name-sorted BAM and the
    # paired FASTQ output is clean.
    samtools view -b -F 2304 ${bam} ${regions} \\
        | samtools sort -n -o extracted_sorted.bam

    # Convert to paired FASTQ; discard singletons (reads whose mate is absent)
    # and any residual non-primary alignments.
    samtools fastq \\
        -F 2304 \\
        -1 ${meta.id}_${meta.safe_species}_R1.fastq.gz \\
        -2 ${meta.id}_${meta.safe_species}_R2.fastq.gz \\
        -0 /dev/null \\
        -s /dev/null \\
        extracted_sorted.bam
    """
}
