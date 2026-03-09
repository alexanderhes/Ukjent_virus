/*
 * HOST_FILTER
 *
 * Remove host (human T2T + PhiX) reads using bowtie2 in unpaired mode with
 * --very-sensitive-local. Each read is judged independently, regardless of
 * mate pairedness, which maximises host-removal sensitivity (Constantinides
 * et al., Bioinformatics 2023). Unmapped reads (-f 4) are retained.
 *
 * Outputs are ENA-ready: host-depleted R1 and R2 FASTQs.
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
    seqkit pair \\
        -1 ${meta.id}_unfilt_R1.fastq.gz \\
        -2 ${meta.id}_unfilt_R2.fastq.gz \\
        -O paired_out

    mv paired_out/${meta.id}_unfilt_R1.fastq.gz ${meta.id}_hostfilt_R1.fastq.gz
    mv paired_out/${meta.id}_unfilt_R2.fastq.gz ${meta.id}_hostfilt_R2.fastq.gz

    rm -rf paired_out
    """
}
