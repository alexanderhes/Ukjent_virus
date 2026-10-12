/*
 * POLISH_SCAFFOLD
 *
 * Turns the per-group templates from SCAFFOLD_MINIMAP2 into genomes, once per
 * sample. Template case marks the source of each base: uppercase = contig
 * (including contig overhangs beyond the reference ends), lowercase = reference
 * fill between and around the contigs. Each template base also carries its
 * reference position ({group_id}.pos), which is kept through every round.
 *   1. All templates of the sample go into one bowtie2 index, so reads choose
 *      between taxa instead of being forced onto one template.
 *   2. Round 1 maps all host-filtered reads (bowtie2 --local); the pairs that
 *      mapped are kept for later rounds. Each round
 *      counts the read bases at every template position (samtools mpileup ->
 *      polish_merge.py counts; no quality filters, overlapping mates both
 *      counted, deletions included in the depth) and merges the reads into the
 *      template (polish_merge.py merge). Where the reads decide, the genome gets
 *      the most frequent read base when >= params.polish_min_depth reads cover
 *      the position (a deletion if deletions are most frequent; the template
 *      base wins ties), and N below that. No 75% threshold and no IUPAC codes:
 *      a mixed position gets its majority base. Read insertions carried by the
 *      majority of the reads are added (samtools consensus -c 0.5 --mark-ins).
 *   3. params.polish_contig_bases sets where the reads decide:
 *        'fallback' (default): everywhere except contig positions with fewer
 *          than polish_min_depth reads, which keep the contig base. SPAdes
 *          contigs can carry a minority allele (metaSPAdes picks a path through
 *          a variant bubble, not the majority), so where enough reads cover a
 *          contig, the reads decide;
 *        'keep': outside the contigs only (every contig base kept);
 *        'reads': everywhere (N below polish_min_depth, also in contigs).
 *      Read insertions between two kept contig bases are ignored. Reference
 *      bases used to fill gaps and ends never reach the genome as such: they
 *      become read-called bases or N (also when no reads map). Read-called
 *      bases are written in lowercase (they count as fill in a later round).
 *      With params.polish_rounds = 0 no reads are used: contig bases only,
 *      every fill position N.
 *   4. bin/polish_merge.py finalize writes each genome twice: full length (end
 *      Ns kept, i.e. the reference length plus any overhangs) and trimmed
 *      (leading/trailing Ns removed). Headers: {id}|{taxon}|{segment}|{ref_accession}.
 *
 * Outputs, published together to validation/{id}/genomes/:
 *   {id}_{group_id}.fasta           : full-length genome per taxon (or segment)
 *   {id}_{group_id}_trimmed.fasta   : the same without end Ns
 *   {id}_genomes.fasta / {id}_genomes_trimmed.fasta : all genomes of the sample
 *   {id}_scaffold_contigs.tsv       : which contigs were used (from SCAFFOLD_GROUPS)
 *   {id}_genome_stats.tsv           : per-group statistics
 *     sample safe_taxon segment group_id ref_accession ref_length n_contigs
 *     n_contigs_used contig_bp contig_cov_pct genome_len n_count N_pct
 *     ref_cov_pct mean_depth trimmed_len n_end_5p n_end_3p n_internal
 *     ref_called_bp overhang_5p overhang_3p
 *   genome_len / n_count / N_pct describe the full-length genome; trimmed_len
 *   the trimmed one; n_end_5p/3p are the Ns trimmed off, n_internal the Ns left.
 *   ref_called_bp = reference positions with a called (non-N) base;
 *   ref_cov_pct = ref_called_bp / ref_length. Overhang and inserted bases
 *   are not counted. mean_depth = the last round's mean read depth (base
 *   counts) over the positions with depth >= params.polish_min_depth.
 *   Groups without a template or without any called base get genome_len 0.
 *
 * For the genome plots (VISUALIZE_GENOMES), not published:
 *   {id}_genome_positions.tsv : group_id pos depth base refpos in_contig -- per
 *                               base of the full-length genomes: the read depth
 *                               that decided it (base counts of the last
 *                               polishing round; 0 with polish_rounds = 0), the
 *                               base (uppercase = contig base kept, lowercase =
 *                               read-called), its reference position and whether
 *                               it lies on a contig position of the unpolished
 *                               template.
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
    tuple val(meta), path("${meta.id}_genomes_out"),                              emit: genomes
    tuple val(meta), path("${meta.id}_genomes_out/${meta.id}_genome_stats.tsv"), emit: stats
    tuple val(meta), path("${meta.id}_genome_positions.tsv"),                    emit: plot_data

    script:
    def rounds    = params.polish_rounds
    def min_depth = params.polish_min_depth
    def contig_bases = params.polish_contig_bases
    if (!(contig_bases in ['keep', 'reads', 'fallback'])) {
        error "polish_contig_bases must be keep, reads or fallback (got '${contig_bases}')"
    }
    """
    genomes=${meta.id}_genomes_out
    mkdir -p "\${genomes}"
    cp ${contig_report} "\${genomes}/"

    # Cased template with reference positions: name <tab> sequence <tab> positions
    : > template.tab
    for fa in ${templates_dir}/*.fa; do
        [ -e "\${fa}" ] || continue
        printf "%s\\t%s\\t%s\\n" "\$(basename "\${fa}" .fa)" "\$(seqkit seq -s -w 0 "\${fa}")" "\$(cat "\${fa%.fa}.pos")" >> template.tab
    done
    cp template.tab template0.tab   # contig positions before polishing (plots)

    # Base counts and template of the last polishing round (depth track and
    # mean_depth): empty without reads or with polish_rounds = 0.
    : > counts.tsv
    cp template.tab template_called.tab
    reads1=${r1}
    reads2=${r2}

    if [ -s template.tab ]; then
        for round in \$(seq 1 ${rounds}); do
            awk -F'\\t' '{ print ">" \$1; print toupper(\$2) }' template.tab > template.fa
            bowtie2-build --threads ${task.cpus} -q template.fa idx_r\${round}
            bowtie2 --local -p ${task.cpus} -x idx_r\${round} -1 "\${reads1}" -2 "\${reads2}" \\
                --no-unal -S round\${round}.sam 2> bowtie2_round\${round}.log
            echo "INFO: ${meta.id} polish round \${round}: \$(tail -1 bowtie2_round\${round}.log)"
            samtools sort -@ ${task.cpus} -o round.bam round\${round}.sam
            samtools index round.bam

            # Later rounds use only the pairs that mapped in round 1.
            if [ "\${round}" -eq 1 ]; then
                samtools collate -u -O round1.sam \\
                    | samtools fastq -n -1 mapped_R1.fastq.gz -2 mapped_R2.fastq.gz \\
                        -0 /dev/null -s /dev/null
                reads1=mapped_R1.fastq.gz
                reads2=mapped_R2.fastq.gz
            fi
            rm -f round\${round}.sam

            # Base counts per template position: the read calls. Counted like
            # samtools consensus simple mode: no quality filters, overlapping
            # mates both counted, non-proper pairs included, no depth cap.
            samtools faidx template.fa
            samtools mpileup -aa -A -B -x -Q 0 -q 0 -d 0 -f template.fa round.bam 2> /dev/null \\
                | polish_merge.py counts > counts.tsv
            # Consensus only for the read insertions (majority of the reads).
            samtools consensus -a --show-ins yes --mark-ins --show-del yes -m simple -c 0.5 \\
                -d ${min_depth} round.bam -o consensus.fa
            # No reads: every read-decided position becomes N.
            cp template.tab template_called.tab
            polish_merge.py merge --template template.tab --cons consensus.fa --out merged.tab \\
                --contig-bases ${contig_bases} --counts counts.tsv --min-depth ${min_depth}
            mv merged.tab template.tab
            if [ ! -s consensus.fa ]; then
                echo "WARNING: no reads mapped to any template of ${meta.id} in round \${round}"
                break
            fi
        done
        # polish_rounds = 0: contigs only, every fill position becomes N.
        if [ ${rounds} -eq 0 ]; then
            : > consensus.fa
            : > counts.tsv
            polish_merge.py merge --template template.tab --cons consensus.fa --counts counts.tsv \\
                --out merged.tab --contig-bases keep
            mv merged.tab template.tab
        fi
    fi

    polish_merge.py finalize --template template.tab --groups ${templates_tsv} --depth counts.tsv \\
        --sample ${meta.id} --outdir "\${genomes}" --final-tab final.tab --min-depth ${min_depth}

    # All genomes of the sample in one multi-FASTA (full length and trimmed)
    if ls "\${genomes}"/${meta.id}_*_trimmed.fasta > /dev/null 2>&1; then
        for f in "\${genomes}"/${meta.id}_*.fasta; do
            case "\${f}" in *_trimmed.fasta) cat "\${f}" >> ${meta.id}_genomes_trimmed.fasta ;; *) cat "\${f}" >> ${meta.id}_genomes.fasta ;; esac
        done
        mv ${meta.id}_genomes.fasta ${meta.id}_genomes_trimmed.fasta "\${genomes}/"
    fi

    # Plot data: each genome base with the read depth that decided it (the last
    # round's base counts), its reference position and whether it is on a contig.
    polish_merge.py positions --final-tab final.tab --template0 template0.tab \\
        --template-called template_called.tab --counts counts.tsv --out ${meta.id}_genome_positions.tsv
    """
}
