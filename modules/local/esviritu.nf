/*
 * ESVIRITU
 *
 * Sensitive, specific virus detection and quantification against the
 * EsViritu curated database of ~20k dereplicated virus genomes.
 *
 * Internal QC (-q False) and deduplication (--dedup False) are disabled —
 * both are handled by upstream HOST_FILTER, FASTP_TRIM, and FASTP_DEDUP.
 * No --filter_dir is passed — host filtering is already complete.
 *
 * --keep True retains intermediate BAM/alignment files in --temp for
 * use by future downstream modules (species-level assembly etc.).
 */

process ESVIRITU {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/esviritu/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path(esviritu_db)   // staged DB directory

    output:
    tuple val(meta),
          path("esv_out/${meta.id}.detected_virus.info.tsv"),
          path("esv_out/${meta.id}.detected_virus.assembly_summary.tsv"),
          path("esv_out/${meta.id}.tax_profile.tsv"),
          path("esv_out/${meta.id}_final_consensus.fasta"),
          path("esv_out/${meta.id}_EsViritu_reactable.html"),   emit: results
    // Glob: captures all files EsViritu writes, used for batch summary
    path "esv_out/*",                                          emit: tsv_files
    // Tuple emit carrying meta — used by validation sub-workflow for splitCsv
    tuple val(meta), path("esv_out/${meta.id}.detected_virus.info.tsv"), emit: tsv_info_meta
    path "${meta.id}_temp",                                    emit: temp, optional: true
    // Third-pass BAM/BAI — retained when keep is active; required for validation
    tuple val(meta), path("${meta.id}_temp/${meta.id}.third.filt.sorted.bam"),     emit: third_bam, optional: true
    tuple val(meta), path("${meta.id}_temp/${meta.id}.third.filt.sorted.bam.bai"), emit: third_bai, optional: true

    script:
    // Force --keep True whenever validate is enabled so BAMs are available
    def keep_flag = (params.esviritu_keep || params.validate) ? '--keep True' : ''
    """
    mkdir -p esv_out

    EsViritu \\
        -r ${r1} ${r2} \\
        -s ${meta.id} \\
        -o esv_out \\
        -p paired \\
        -q False \\
        --dedup False \\
        -t ${task.cpus} \\
        --db ${esviritu_db} \\
        --temp ${meta.id}_temp \\
        ${keep_flag}
    """
}
