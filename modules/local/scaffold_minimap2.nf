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
 *   3. Builds a contig consensus along the reference with samtools consensus;
 *      overlapping contigs merge, uncovered positions become N. Contig
 *      insertions relative to the reference are kept (--mark-ins) and contig
 *      deletions dropped, so the template holds the contig sequence itself.
 *   4. bin/build_scaffold_template.py fills the N positions with the reference
 *      base in lowercase to make the polishing template, so reads can map
 *      across the gaps between contigs. Case marks the source: uppercase =
 *      contig base, kept as is in the final genome; lowercase = reference fill,
 *      which POLISH_SCAFFOLD replaces with read-called bases (or N). Reference
 *      bases never reach the genome. It also adds contig overhangs beyond the
 *      reference ends (references often lack UTR ends) when the alignment
 *      reaches that end, the overhang has no host BLAST hit, and after removing
 *      a terminal homopolymer it is <= params.scaffold_max_overhang bp; and it
 *      records the reference position of every template base.
 *
 * Outputs (consumed by POLISH_SCAFFOLD):
 *   {id}_minimap2_templates/{group_id}.fa  : one template per group with aligned contigs
 *   {id}_minimap2_templates/{group_id}.pos : reference position of each template base
 *                                            (comma-separated; see build_scaffold_template.py)
 *   {id}_minimap2_templates.tsv            : groups TSV plus n_contigs_used, contig_bp
 *                                            (reference bases covered by contigs),
 *                                            overhang_5p and overhang_3p (bp added)
 * and, for the genome coverage plot (VISUALIZE_GENOMES):
 *   {id}_contig_ref.tsv : group_id contig ref_accession ref_start ref_end strand
 *                         one row per primary/supplementary contig alignment
 *                         (1-based, inclusive reference coordinates)
 */

process SCAFFOLD_MINIMAP2 {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(contigs), path(groups), path(host_hits)   // host_hits: FILTER_HOST_CONTIGS BLAST hits, or [] when the filter is off
    path(db_dir)

    output:
    tuple val(meta), path("${meta.id}_minimap2_templates"), path("${meta.id}_minimap2_templates.tsv"), emit: templates
    tuple val(meta), path("${meta.id}_contig_ref.tsv"),                                              emit: contig_aln

    script:
    def host_arg = host_hits ? "--host-blast ${host_hits}" : ''
    """
    fna=\$(find -L ${db_dir} -name "*.fna" | head -1)
    outdir=${meta.id}_minimap2_templates
    mkdir -p "\${outdir}"

    printf "%s\\tn_contigs_used\\tcontig_bp\\toverhang_5p\\toverhang_3p\\n" "\$(head -1 ${groups})" > ${meta.id}_minimap2_templates.tsv
    printf "group_id\\tcontig\\tref_accession\\tref_start\\tref_end\\tstrand\\n" > ${meta.id}_contig_ref.tsv

    # Read via \\037: tab is IFS whitespace, so empty fields (segment) would collapse.
    tail -n +2 ${groups} | tr '\\t' '\\037' | while IFS=\$'\\037' read -r sample taxon segment group_id ref_acc ref_len n_contigs contig_ids; do
        seqkit grep -p "\${ref_acc}" "\${fna}" > ref.fa
        echo "\${contig_ids}" | tr ',' '\\n' > contig_ids.txt
        seqkit grep -f contig_ids.txt ${contigs} > group_contigs.fa

        minimap2 -ax map-ont -k 11 -w 5 --secondary=no -t ${task.cpus} ref.fa group_contigs.fa \\
            | samtools sort -o group.bam -
        samtools index group.bam

        n_used=\$(samtools view -F 0x904 group.bam | cut -f1 | sort -u | wc -l)

        # Reference span of each contig alignment: M, D, N, = and X consume the reference.
        samtools view -F 0x104 group.bam | awk -F'\\t' -v OFS='\\t' -v gid="\${group_id}" '{
            ref_len = 0; cig = \$6
            while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
                n = substr(cig, 1, RLENGTH - 1); op = substr(cig, RLENGTH, 1)
                if (op ~ /[MDN=X]/) ref_len += n
                cig = substr(cig, RLENGTH + 1)
            }
            print gid, \$1, \$3, \$4, \$4 + ref_len - 1, (int(\$2 / 16) % 2 ? "-" : "+")
        }' >> ${meta.id}_contig_ref.tsv

        contig_bp=0
        overhang_5p=0
        overhang_3p=0

        if [ "\${n_used}" -gt 0 ]; then
            samtools consensus -a --show-ins yes --mark-ins --show-del yes -m simple -d 1 -c 0.5 \\
                group.bam -o contig_cons.fa

            samtools view -F 0x104 group.bam > group.sam
            build_scaffold_template.py \\
                --ref ref.fa --cons contig_cons.fa --sam group.sam --group "\${group_id}" \\
                --out-fa "\${outdir}/\${group_id}.fa" --out-pos "\${outdir}/\${group_id}.pos" \\
                --max-overhang ${params.scaffold_max_overhang} \\
                --host-min-identity ${params.host_contig_min_identity} ${host_arg} > template_stats.txt
            read -r contig_bp overhang_5p overhang_3p < template_stats.txt
        else
            echo "WARNING: no contig of \${group_id} aligned to \${ref_acc} -- no template"
        fi

        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "\${sample}" "\${taxon}" "\${segment}" "\${group_id}" "\${ref_acc}" \\
            "\${ref_len}" "\${n_contigs}" "\${contig_ids}" "\${n_used}" "\${contig_bp}" \\
            "\${overhang_5p}" "\${overhang_3p}" \\
            >> ${meta.id}_minimap2_templates.tsv
    done
    """
}
