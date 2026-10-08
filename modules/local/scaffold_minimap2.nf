/*
 * SCAFFOLD_MINIMAP2
 *
 * Reference-guided combination of the contigs of each scaffolding group
 * (one taxon, or one segment of a segmented taxon; see SCAFFOLD_GROUPS):
 *   1. Extracts the group's reference accession from the EsViritu .fna and
 *      the group's contigs from the sample's SPAdes contigs.
 *   2. Aligns the contigs to the reference with minimap2 -x map-ont. The
 *      assembly preset (asm20) fails on contigs ~80-86% identical to their
 *      reference (e.g. Puumala: L contig unaligned, M/S only 68-74% aligned),
 *      and unaligned contig parts never reach the consensus. map-ont is tuned
 *      for ~85-90% identity and aligns such contigs end to end. Shorter seeds
 *      (-k 11 -w 5 instead of map-ont's k=15) are needed in hypervariable
 *      stretches: without them minimap2 finds no seeds in the last ~400 bp of
 *      a Puumala S contig and soft-clips it.
 *   3. Builds a contig consensus in reference coordinates with samtools
 *      consensus; overlapping contigs merge, uncovered positions become N.
 *   4. Fills the N positions with the reference base in lowercase to make the
 *      polishing template, so reads can map across the gaps between contigs.
 *      Reference bases never reach the final genome: POLISH_SCAFFOLD keeps only
 *      read-supported bases.
 *
 * Outputs (consumed by POLISH_SCAFFOLD):
 *   {id}_minimap2_templates/{group_id}.fa  : one template per group with aligned contigs
 *   {id}_minimap2_templates.tsv            : groups TSV plus n_contigs_used and
 *                                            contig_bp (reference bases covered by contigs)
 */

process SCAFFOLD_MINIMAP2 {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(contigs), path(groups)
    path(db_dir)

    output:
    tuple val(meta), path("${meta.id}_minimap2_templates"), path("${meta.id}_minimap2_templates.tsv"), emit: templates

    script:
    """
    fna=\$(find -L ${db_dir} -name "*.fna" | head -1)
    outdir=${meta.id}_minimap2_templates
    mkdir -p "\${outdir}"

    printf "%s\\tn_contigs_used\\tcontig_bp\\n" "\$(head -1 ${groups})" > ${meta.id}_minimap2_templates.tsv

    # Read via \\037: tab is IFS whitespace, so empty fields (segment) would collapse.
    tail -n +2 ${groups} | tr '\\t' '\\037' | while IFS=\$'\\037' read -r sample taxon segment group_id ref_acc ref_len n_contigs contig_ids; do
        seqkit grep -p "\${ref_acc}" "\${fna}" > ref.fa
        echo "\${contig_ids}" | tr ',' '\\n' > contig_ids.txt
        seqkit grep -f contig_ids.txt ${contigs} > group_contigs.fa

        minimap2 -ax map-ont -k 11 -w 5 --secondary=no -t ${task.cpus} ref.fa group_contigs.fa \\
            | samtools sort -o group.bam -
        samtools index group.bam

        n_used=\$(samtools view -F 0x904 group.bam | cut -f1 | sort -u | wc -l)
        contig_bp=0

        if [ "\${n_used}" -gt 0 ]; then
            samtools consensus -a --show-ins no --show-del yes -m simple -d 1 -c 0.5 \\
                group.bam -o contig_cons.fa

            # Consensus is in reference coordinates (insertions hidden), so position i
            # of the consensus is position i of the reference. N -> lowercase reference
            # base; deletion markers (*) are dropped.
            paste <(seqkit seq -s -w 0 contig_cons.fa) <(seqkit seq -s -w 0 ref.fa) \\
                | awk -F'\\t' -v name="\${group_id}" '
                    {
                        c = toupper(\$1); r = tolower(\$2); out = ""; covered = 0
                        for (i = 1; i <= length(c); i++) {
                            b = substr(c, i, 1)
                            if (b == "N") b = substr(r, i, 1)
                            else covered++
                            if (b != "*") out = out b
                        }
                        print ">" name > "/dev/stdout"
                        print out > "/dev/stdout"
                        print covered > "covered.txt"
                    }' > "\${outdir}/\${group_id}.fa"
            contig_bp=\$(cat covered.txt)
        else
            echo "WARNING: no contig of \${group_id} aligned to \${ref_acc} -- no template"
        fi

        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "\${sample}" "\${taxon}" "\${segment}" "\${group_id}" "\${ref_acc}" \\
            "\${ref_len}" "\${n_contigs}" "\${contig_ids}" "\${n_used}" "\${contig_bp}" \\
            >> ${meta.id}_minimap2_templates.tsv
    done
    """
}
