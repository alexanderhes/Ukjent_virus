/*
 * POLISH_SCAFFOLD
 *
 * Turns the per-group templates from SCAFFOLD_MINIMAP2 into read-supported
 * genomes, once per sample:
 *   1. All templates of the sample go into one bowtie2 index, so reads choose
 *      between taxa instead of being forced onto one template.
 *   2. Round 1 maps all host-filtered reads (bowtie2 --local); later rounds map
 *      only the pairs that mapped in round 1. Each round calls a new template
 *      with samtools consensus (insertions applied, deletions removed).
 *      Positions with depth < params.polish_min_depth become N, so only
 *      read-supported bases reach the output -- reference bases used to fill
 *      template gaps and ends never do.
 *   3. Leading/trailing Ns are trimmed and headers renamed to
 *      {id}|{taxon}|{segment}|{ref_accession}.
 *
 * Outputs, published together to validation/{id}/genomes/:
 *   {id}_{group_id}.fasta   : one polished genome per taxon (or segment)
 *   {id}_genomes.fasta      : all genomes of the sample in one multi-FASTA
 *   {id}_scaffold_contigs.tsv : which contigs were used (from SCAFFOLD_GROUPS)
 *   {id}_genome_stats.tsv   : per-group statistics
 *     sample safe_taxon segment group_id ref_accession ref_length n_contigs
 *     n_contigs_used contig_bp contig_cov_pct genome_len n_count N_pct
 *     ref_cov_pct mean_depth
 *   ref_cov_pct = called (non-N) bases / reference length, capped at 100.
 *   mean_depth is the final round's mean read depth over the called positions
 *   (depth >= params.polish_min_depth), i.e. over the bases of the genome.
 *   Groups without a template or without mapped reads get genome_len 0.
 */

process POLISH_SCAFFOLD {
    tag "${meta.id}"
    label 'process_medium'

    // Publish only the folder (the stats file inside it is emitted separately).
    publishDir "${params.outdir}/validation/${meta.id}", mode: 'copy',
        pattern: "*_genomes_out",
        saveAs: { fn -> 'genomes' }

    input:
    tuple val(meta), path(templates_dir), path(templates_tsv), path(contig_report), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_genomes_out"),                       emit: genomes
    tuple val(meta), path("${meta.id}_genomes_out/${meta.id}_genome_stats.tsv"), emit: stats

    script:
    def rounds    = params.polish_rounds
    def min_depth = params.polish_min_depth
    """
    genomes=${meta.id}_genomes_out
    mkdir -p "\${genomes}"
    cp ${contig_report} "\${genomes}/"

    if ls ${templates_dir}/*.fa > /dev/null 2>&1; then
        cat ${templates_dir}/*.fa | seqkit seq -u > template.fa
    else
        : > template.fa
    fi

    : > final_depth.tsv
    : > polished.fa

    if [ -s template.fa ]; then
        reads1=${r1}
        reads2=${r2}
        for round in \$(seq 1 ${rounds}); do
            bowtie2-build --threads ${task.cpus} -q template.fa idx_r\${round}
            bowtie2 --local -p ${task.cpus} -x idx_r\${round} -1 "\${reads1}" -2 "\${reads2}" \\
                --no-unal -S round\${round}.sam 2> bowtie2_round\${round}.log
            echo "INFO: ${meta.id} polish round \${round}: \$(tail -1 bowtie2_round\${round}.log)"
            samtools sort -@ ${task.cpus} -o round.bam round\${round}.sam
            samtools index round.bam

            if [ "\${round}" -eq 1 ] && [ ${rounds} -gt 1 ]; then
                samtools collate -u -O round1.sam \\
                    | samtools fastq -n -1 mapped_R1.fastq.gz -2 mapped_R2.fastq.gz \\
                        -0 /dev/null -s /dev/null
                reads1=mapped_R1.fastq.gz
                reads2=mapped_R2.fastq.gz
            fi
            rm -f round\${round}.sam

            samtools depth -aa round.bam > final_depth.tsv
            samtools consensus -a --show-ins yes --show-del no -m simple \\
                -d ${min_depth} round.bam -o consensus.fa
            seqkit seq -u consensus.fa > template.fa
            if [ ! -s template.fa ]; then
                echo "WARNING: no reads mapped to any template of ${meta.id} in round \${round}"
                break
            fi
        done
        cp template.fa polished.fa
    fi

    # Trim end Ns; write one FASTA per group and its length / N count.
    seqkit fx2tab polished.fa | awk -F'\\t' '{ s = \$2; sub(/^N+/, "", s); sub(/N+\$/, "", s); print \$1 "\\t" s }' \\
        > polished.tab

    awk -F'\\t' -v OFS='\\t' -v sample="${meta.id}" -v outdir="\${genomes}" -v min_depth=${min_depth} '
        FILENAME == ARGV[1] { seq[\$1] = \$2; next }
        FILENAME == ARGV[2] { if (\$3 >= min_depth) { dsum[\$1] += \$3; dn[\$1]++ }; next }
        FNR == 1 {
            print "sample", "safe_taxon", "segment", "group_id", "ref_accession",
                  "ref_length", "n_contigs", "n_contigs_used", "contig_bp", "contig_cov_pct",
                  "genome_len", "n_count", "N_pct", "ref_cov_pct", "mean_depth"
            next
        }
        {
            taxon = \$2; segment = \$3; gid = \$4; ref = \$5; ref_len = \$6
            n_contigs = \$7; n_used = \$9; contig_bp = \$10
            s = (gid in seq) ? seq[gid] : ""
            len = length(s)
            tmp = s; n_count = gsub(/N/, "", tmp)
            called = len - n_count
            if (len > 0) {
                file = outdir "/" sample "_" gid ".fasta"
                print ">" sample "|" taxon "|" segment "|" ref > file
                for (i = 1; i <= len; i += 60) print substr(s, i, 60) > file
                close(file)
            }
            ref_cov = (ref_len > 0) ? called / ref_len * 100 : 0
            if (ref_cov > 100) ref_cov = 100
            print sample, taxon, segment, gid, ref, ref_len, n_contigs, n_used, contig_bp,
                  (ref_len > 0 ? sprintf("%.1f", contig_bp / ref_len * 100) : "NA"),
                  len, n_count,
                  (len > 0 ? sprintf("%.1f", n_count / len * 100) : "NA"),
                  sprintf("%.1f", ref_cov),
                  (gid in dn ? sprintf("%.1f", dsum[gid] / dn[gid]) : 0)
        }
    ' polished.tab final_depth.tsv ${templates_tsv} > \${genomes}/${meta.id}_genome_stats.tsv

    # All genomes of the sample in one multi-FASTA
    if ls \${genomes}/${meta.id}_*.fasta > /dev/null 2>&1; then
        cat \${genomes}/${meta.id}_*.fasta > ${meta.id}_genomes.fasta
        mv ${meta.id}_genomes.fasta \${genomes}/
    fi
    """
}
