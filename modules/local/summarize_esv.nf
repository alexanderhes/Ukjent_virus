/*
 * SUMMARIZE_ESV
 *
 * Runs after ALL per-sample ESVIRITU processes have finished (enforced by
 * .collect() in main.nf). Stages all per-sample output files in a flat
 * directory, then calls summarize_esv_runs to produce the batch-level
 * summary tables and interactive HTML report.
 */

process SUMMARIZE_ESV {
    label 'process_low'

    publishDir "${params.outdir}/esviritu_batch", mode: 'copy'

    input:
    path(sample_files)   // collected list of all output files from all ESVIRITU runs

    output:
    path "esv_summary/*",                                                        emit: summary
    path "esv_summary/esv_staged.detected_virus.info.tsv",                      emit: info_tsv
    path "esv_summary/esv_staged.detected_virus.assembly_summary.tsv",          emit: assembly_summary_tsv

    script:
    """
    # Stage all per-sample files into one flat directory.
    # Nextflow has already staged every file in the process workDir.
    mkdir -p esv_staged

    for f in ${sample_files}; do
        cp "\$f" esv_staged/ || echo "WARNING: failed to copy \$f" >&2
    done

    # Run batch summarisation against the staging directory
    summarize_esv_runs esv_staged --outdir esv_summary
    """
}
