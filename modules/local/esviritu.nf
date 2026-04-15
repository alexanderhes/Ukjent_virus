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
 * --keep True retains intermediate BAM/alignment files in --temp.
 * Enable via params.esviritu_keep when BAMs are needed for other purposes.
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
    // --keep True retains intermediate temp BAMs; driven solely by esviritu_keep param
    def keep_flag = params.esviritu_keep ? '--keep True' : ''
    """
    mkdir -p esv_out

    read_count=\$(seqkit stats -T -j 1 ${r1} 2>/dev/null | awk 'NR==2{print \$4+0}')

    if [ "\${read_count}" -gt 0 ]; then
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
            ${keep_flag} || true

        # EsViritu may exit 0 but produce no output when no reads map to the DB
        if [ ! -f "esv_out/${meta.id}.detected_virus.info.tsv" ]; then
            echo "WARNING: ${meta.id} — EsViritu produced no output (no reads mapped to DB) — writing placeholder files"
            printf 'sample_ID\tName\tdescription\tLength\tSegment\tAccession\tAssembly\tAsm_length\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tRPKMF\tread_count\tcovered_bases\tmean_coverage\tavg_read_identity\tPi\tfiltered_reads_in_sample\n' \
                > esv_out/${meta.id}.detected_virus.info.tsv
            printf 'sample_ID\tfiltered_reads_in_sample\tAssembly\tAsm_length\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tread_count\tcovered_bases\tavg_read_identity\tAccession\tSegment\tRPKMF\n' \
                > esv_out/${meta.id}.detected_virus.assembly_summary.tsv
            printf 'sample_ID\tfiltered_reads_in_sample\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tread_count\tRPKMF\tavg_read_identity\tassembly_list\n' \
                > esv_out/${meta.id}.tax_profile.tsv
            touch esv_out/${meta.id}_final_consensus.fasta
            touch esv_out/${meta.id}_EsViritu_reactable.html
        fi
    else
        # No reads after host filtering — produce header-only placeholder outputs
        # so all downstream processes receive valid (empty) files.
        echo "WARNING: ${meta.id} has 0 reads after host filtering — skipping EsViritu"

        printf 'sample_ID\tName\tdescription\tLength\tSegment\tAccession\tAssembly\tAsm_length\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tRPKMF\tread_count\tcovered_bases\tmean_coverage\tavg_read_identity\tPi\tfiltered_reads_in_sample\n' \
            > esv_out/${meta.id}.detected_virus.info.tsv

        printf 'sample_ID\tfiltered_reads_in_sample\tAssembly\tAsm_length\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tread_count\tcovered_bases\tavg_read_identity\tAccession\tSegment\tRPKMF\n' \
            > esv_out/${meta.id}.detected_virus.assembly_summary.tsv

        printf 'sample_ID\tfiltered_reads_in_sample\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tread_count\tRPKMF\tavg_read_identity\tassembly_list\n' \
            > esv_out/${meta.id}.tax_profile.tsv

        touch esv_out/${meta.id}_final_consensus.fasta
        touch esv_out/${meta.id}_EsViritu_reactable.html
    fi
    """
}
