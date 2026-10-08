/*
 * FILTER_HOST_CONTIGS
 *
 * Removes host (human) contigs from the SPAdes assembly before BLAST.
 *
 * Purpose: remove human rRNA that assembles into contigs. Reads from the
 * multi-copy rRNA arrays (chr13/14/15/21/22) can pass HOST_FILTER because their
 * mates align to different repeat copies and the pair is never concordant.
 * Such contigs then BLAST-hit viral references whose GenBank records carry
 * human rRNA at their ends, producing recurring artefact detections and false
 * genomes in the scaffolding steps.
 *
 * Method: megablast of the contigs against a BLAST database of the host
 * reference (params.host_blastdb; built with makeblastdb from the same FASTA
 * as the bowtie2 host index). For each contig the query bases covered by host
 * HSPs with identity >= params.host_contig_min_identity are merged; a contig is
 * removed when they cover >= params.host_contig_max_pct % of its length.
 * Repeats (rRNA arrays) give many HSPs per contig; -max_target_seqs and
 * -max_hsps bound that without affecting the covered fraction.
 *
 * Outputs:
 *   {id}_query.nohost.fasta  : contigs kept (input to BLAST and scaffolding)
 *   {id}_host_contigs.tsv    : per-contig report: contig, length, host_cov_bp,
 *                              host_cov_pct, best_host_hit, best_identity, removed
 *   {id}_host_blast.tsv      : raw host BLAST hits (outfmt 6 + qlen)
 */

process FILTER_HOST_CONTIGS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/validation/${meta.id}/assembly", mode: 'copy'

    input:
    tuple val(meta), path(query)
    path(blastdb_files)   // all host BLAST database files staged into work dir

    output:
    tuple val(meta), path("${meta.id}_query.nohost.fasta"), emit: query
    tuple val(meta), path("${meta.id}_host_contigs.tsv"),   emit: report
    path "${meta.id}_host_blast.tsv",                        emit: blast

    script:
    def db           = file(params.host_blastdb).name
    def max_pct      = params.host_contig_max_pct
    def min_identity = params.host_contig_min_identity
    """
    echo -e "contig\\tlength\\thost_cov_bp\\thost_cov_pct\\tbest_host_hit\\tbest_identity\\tremoved" > ${meta.id}_host_contigs.tsv
    : > ${meta.id}_host_blast.tsv

    if [ \$(grep -c "^>" ${query} 2>/dev/null || echo 0) -eq 0 ]; then
        cp ${query} ${meta.id}_query.nohost.fasta
        exit 0
    fi

    blastn \\
        -task megablast \\
        -query ${query} \\
        -db ${db} \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \\
        -evalue 1e-10 \\
        -max_target_seqs 10 \\
        -max_hsps 50 \\
        -num_threads ${task.cpus} \\
        -out ${meta.id}_host_blast.tsv

    # Host-covered query bases per contig: merge HSP intervals (identity filter),
    # and keep the best-scoring host hit for the report.
    awk -F'\\t' -v OFS='\\t' -v min_id=${min_identity} '\$3 >= min_id { print \$1, \$7, \$8 }' ${meta.id}_host_blast.tsv \\
        | sort -k1,1 -k2,2n \\
        | awk -F'\\t' -v OFS='\\t' '
            function flush() { if (cur != "") cov[cur] += e - s + 1 }
            \$1 != cur { flush(); cur = \$1; s = \$2; e = \$3; next }
            \$2 > e    { cov[cur] += e - s + 1; s = \$2; e = \$3; next }
                       { if (\$3 > e) e = \$3 }
            END        { flush(); for (c in cov) print c, cov[c] }
        ' > host_cov.tsv
    sort -t \$'\\t' -k1,1 -k12,12gr ${meta.id}_host_blast.tsv | awk -F'\\t' -v OFS='\\t' '!seen[\$1]++ { print \$1, \$2, \$3 }' > best_hit.tsv

    seqkit fx2tab -n -i -l ${query} \\
        | awk -F'\\t' -v OFS='\\t' -v max_pct=${max_pct} '
            FILENAME == ARGV[1] { cov[\$1] = \$2; next }
            FILENAME == ARGV[2] { hit[\$1] = \$2; idn[\$1] = \$3; next }
            {
                id = \$1; len = \$2; c = (id in cov) ? cov[id] : 0
                pct = (len > 0) ? c / len * 100 : 0
                print id, len, c, sprintf("%.1f", pct),
                      (id in hit) ? hit[id] : "NA", (id in idn) ? idn[id] : "NA",
                      (pct >= max_pct) ? "true" : "false"
            }
        ' host_cov.tsv best_hit.tsv - >> ${meta.id}_host_contigs.tsv

    awk -F'\\t' 'NR > 1 && \$7 == "true" { print \$1 }' ${meta.id}_host_contigs.tsv > removed_ids.txt
    if [ -s removed_ids.txt ]; then
        seqkit grep -v -f removed_ids.txt ${query} > ${meta.id}_query.nohost.fasta
    else
        cp ${query} ${meta.id}_query.nohost.fasta
    fi

    echo "INFO: ${meta.id}: removed \$(wc -l < removed_ids.txt) of \$(tail -n +2 ${meta.id}_host_contigs.tsv | wc -l) contigs as host (>= ${max_pct}% of length in host BLAST hits >= ${min_identity}% identity)"
    """
}
