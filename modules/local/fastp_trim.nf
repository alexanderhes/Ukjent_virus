/*
 * FASTP_TRIM
 *
 * Pass 1 — quality trimming only (NO deduplication).
 * Adapter auto-detection, quality filtering, length filtering,
 * low-complexity filtering, and poly-G/X trimming.
 *
 * Deduplication is intentionally deferred to FASTP_DEDUP so that it
 * operates on already-trimmed reads, exposing the full ~40-50% duplicate
 * population rather than the ~5% visible in raw reads.
 */

process FASTP_TRIM {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/fastp_trim/${meta.id}", mode: 'copy',
        pattern: "*.{json,html,log}"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_trim_R1.fastq.gz"),
                     path("${meta.id}_trim_R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_fastp_trim.json"),  emit: json
    path "${meta.id}_fastp_trim.html",                    emit: html
    path "${meta.id}_fastp_trim.log",                     emit: log

    script:
    """
    read_count=\$(seqkit stats -T -j 1 ${r1} 2>/dev/null | awk 'NR==2{print \$4+0}')

    if [ "\${read_count}" -gt 0 ]; then
        fastp \\
            -i ${r1} \\
            -I ${r2} \\
            -o ${meta.id}_trim_R1.fastq.gz \\
            -O ${meta.id}_trim_R2.fastq.gz \\
            --detect_adapter_for_pe \\
            -q ${params.quality_cutoff} \\
            --length_required ${params.min_length} \\
            --low_complexity_filter \\
            --complexity_threshold ${params.complexity_threshold} \\
            --trim_poly_g \\
            --trim_poly_x \\
            --thread ${task.cpus} \\
            --json ${meta.id}_fastp_trim.json \\
            --html ${meta.id}_fastp_trim.html \\
            2> ${meta.id}_fastp_trim.log
    else
        echo "WARNING: ${meta.id} has 0 reads — skipping FASTP_TRIM"
        cp ${r1} ${meta.id}_trim_R1.fastq.gz
        cp ${r2} ${meta.id}_trim_R2.fastq.gz
        printf '{"summary":{"before_filtering":{"total_reads":0},"after_filtering":{"total_reads":0}}}' \\
            > ${meta.id}_fastp_trim.json
        touch ${meta.id}_fastp_trim.html ${meta.id}_fastp_trim.log
    fi
    """
}
