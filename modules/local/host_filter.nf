/*
 * HOST_FILTER
 *
 * Remove host (human T2T + PhiX) reads using bowtie2.
 *
 * Default (sensitive_host_filter = false):
 *   Single paired-end run with --sensitive-local --no-discordant --no-mixed.
 *   Unmapped pairs are written directly via --un-conc-gz. Faster (3–5×) and
 *   produces synchronised R1/R2 output without a re-pairing step.
 *
 * High-sensitivity mode (sensitive_host_filter = true):
 *   bowtie2 is run independently on R1 and R2 (unpaired, --very-sensitive-local)
 *   so each read is judged on its own merits (Constantinides et al.,
 *   Bioinformatics 2023). seqkit pair then re-synchronises the output.
 *   Use when minimising host read carry-through matters more than speed.
 */

process HOST_FILTER {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/host_filtered/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path(bt2_index)   // all *.bt2 index files staged into work dir

    output:
    tuple val(meta), path("${meta.id}_hostfilt_R1.fastq.gz"),
                     path("${meta.id}_hostfilt_R2.fastq.gz"), emit: reads
    path "${meta.id}_hostfilt.log",                           emit: log
    tuple val(meta), path("${meta.id}_raw_reads.txt"),        emit: raw_read_count

    script:
    // Derive the bowtie2 prefix from the staged file names (e.g. host_index)
    def index = file(params.host_index).name
    if (params.sensitive_host_filter) {
        """
        # ── Count raw reads (R1 count × 2 for paired) ───────────────────────────
        seqkit stats -T -j ${task.cpus} ${r1} \
            | awk 'NR==2{print \$4 * 2}' > ${meta.id}_raw_reads.txt

        # ── R1 — unmapped reads written directly via --un-gz ────────────────────
        bowtie2 \\
            --very-sensitive-local \\
            -x ${index} \\
            -U ${r1} \\
            -p ${task.cpus} \\
            --un-gz ${meta.id}_unfilt_R1.fastq.gz \\
            --no-unal \\
            -S /dev/null \\
            2>> ${meta.id}_hostfilt.log

        # ── R2 — unmapped reads written directly via --un-gz ────────────────────
        bowtie2 \\
            --very-sensitive-local \\
            -x ${index} \\
            -U ${r2} \\
            -p ${task.cpus} \\
            --un-gz ${meta.id}_unfilt_R2.fastq.gz \\
            --no-unal \\
            -S /dev/null \\
            2>> ${meta.id}_hostfilt.log

        # ── Re-pair: synchronise R1 and R2 by read name ─────────────────────────
        # Independent per-read filtering means a read may be present in one file
        # but host-filtered away in the other. seqkit pair matches by read name
        # and discards singletons, producing two guaranteed-synchronised FASTQs.
        # Guard: if either file is empty (all reads were host-filtered) skip pairing
        # to avoid seqkit pair crashing on an empty gzip (EOF error).
        r1_count=\$(seqkit stats -T ${meta.id}_unfilt_R1.fastq.gz 2>/dev/null | awk 'NR==2{print \$4+0}')
        r2_count=\$(seqkit stats -T ${meta.id}_unfilt_R2.fastq.gz 2>/dev/null | awk 'NR==2{print \$4+0}')

        if [ "\${r1_count}" -gt 0 ] && [ "\${r2_count}" -gt 0 ]; then
            seqkit pair \\
                -1 ${meta.id}_unfilt_R1.fastq.gz \\
                -2 ${meta.id}_unfilt_R2.fastq.gz \\
                -O paired_out

            mv paired_out/${meta.id}_unfilt_R1.fastq.gz ${meta.id}_hostfilt_R1.fastq.gz
            mv paired_out/${meta.id}_unfilt_R2.fastq.gz ${meta.id}_hostfilt_R2.fastq.gz
            rm -rf paired_out
        else
            # All reads were host-filtered; rename as-is so output files are present
            # and downstream processes receive empty-but-valid FASTQs
            mv ${meta.id}_unfilt_R1.fastq.gz ${meta.id}_hostfilt_R1.fastq.gz
            mv ${meta.id}_unfilt_R2.fastq.gz ${meta.id}_hostfilt_R2.fastq.gz
        fi
        """
    } else {
        """
        # ── Count raw reads (R1 count × 2 for paired) ───────────────────────────
        seqkit stats -T -j ${task.cpus} ${r1} \
            | awk 'NR==2{print \$4 * 2}' > ${meta.id}_raw_reads.txt

        # ── Single paired-end run; unmapped pairs written via --un-conc-gz ──────
        # --no-discordant and --no-mixed prevent bowtie2 from attempting solo-read
        # alignments when a pair fails, keeping --un-conc-gz as the sole output path.
        bowtie2 \\
            --sensitive-local \\
            -x ${index} \\
            -1 ${r1} \\
            -2 ${r2} \\
            -p ${task.cpus} \\
            --no-discordant \\
            --no-mixed \\
            --un-conc-gz ${meta.id}_unfilt.fastq.gz \\
            --no-unal \\
            -S /dev/null \\
            2> ${meta.id}_hostfilt.log

        mv ${meta.id}_unfilt.fastq.gz.1 ${meta.id}_hostfilt_R1.fastq.gz
        mv ${meta.id}_unfilt.fastq.gz.2 ${meta.id}_hostfilt_R2.fastq.gz
        """
    }
}
