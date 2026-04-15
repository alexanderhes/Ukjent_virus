/*
 * BLASTN_VALIDATE
 *
 * Assembles all quality-filtered reads per sample with SPAdes (upstream), then:
 *   1. Builds an accession → safe_taxon lookup from the full EsViritu metadata TSV
 *      using the configured assembly_taxon_level (subspecies / species).
 *   2. Derives the set of detected safe_taxa from the per-sample
 *      detected_virus.info.tsv produced by ESVIRITU.
 *   3. Builds a BLAST nucleotide DB from the full EsViritu .fna.
 *   4. BLASTs the whole-sample contigs against that DB.
 *   5. Assigns each contig exclusively to the detected safe_taxon that
 *      accumulates the highest total BLAST bitscore across all its hits.
 *      Hits to accessions outside the detected safe_taxa are discarded.
 *   6. Writes per-sample output files compatible with the downstream R scripts:
 *      - {id}_blastn.tsv      : hits with Sample and Species (safe_taxon) columns
 *      - {id}_ref_lengths.tsv : lengths of all reference sequences for detected taxa
 *      - {id}_has_contigs.txt : one row per detected safe_taxon (true/false)
 *
 * Output format columns (blastn.tsv):
 *   Sample  Species  Scaffold_ID  Matched_Reference  Identity_%  Align_Len
 *   Query_Len  Mismatches  Gap_Opens  Q_Start  Q_End  S_Start  S_End
 *   E-value  Bit_Score  Cov_%
 */

process BLASTN_VALIDATE {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/validation/${meta.id}/blast", mode: 'copy'

    input:
    tuple val(meta), path(query), path(detected_info_tsv)
    path(db_dir)    // staged EsViritu DB directory (contains .fna and metadata TSV)

    output:
    tuple val(meta), path("${meta.id}_blastn.tsv"),     emit: blast_results
    path "${meta.id}_ref_lengths.tsv",                  emit: ref_lengths
    path "${meta.id}_has_contigs.txt",                  emit: has_contigs

    script:
    def taxon_level = params.assembly_taxon_level
    """
    fna=\$(find -L ${db_dir} -name "*.fna"      | head -1)
    meta_tsv=\$(find -L ${db_dir} -name "*.tsv" | head -1)

    # ── Step 1: Build full accession → safe_taxon lookup from the metadata TSV ──
    # Metadata TSV columns: 1=Accession  11=species  12=subspecies
    # Safe transform: strip ^[st]__ prefix, then replace [^A-Za-z0-9._-] with _
    # IMPORTANT: this transform MUST stay in sync with to_safe() in
    # bin/make_overview_table.R and the Groovy replaceAll in main.nf.
    awk -F'\\t' -v level="${taxon_level}" '
        BEGIN { OFS="\\t" }
        NR==1 { next }
        {
            acc        = \$1
            species    = \$11
            subspecies = \$12
            taxon = (level == "subspecies" && subspecies != "" && subspecies != "NA") \\
                    ? subspecies : species
            sub(/^[st]__/, "", taxon)
            gsub(/[^A-Za-z0-9._-]/, "_", taxon)
            if (taxon != "") print acc, taxon
        }
    ' "\${meta_tsv}" > acc_safe_taxon_full.tsv

    # ── Step 2: Derive detected safe_taxa from the ESVIRITU info TSV ────────────
    awk -F'\\t' -v level="${taxon_level}" '
        BEGIN { OFS="\\t" }
        NR==1 {
            for (i=1; i<=NF; i++) {
                if (\$i == "species")    spcol  = i
                if (\$i == "subspecies") sspcol = i
            }
            next
        }
        {
            species    = \$spcol
            subspecies = \$sspcol
            taxon = (level == "subspecies" && subspecies != "" && subspecies != "NA") \\
                    ? subspecies : species
            sub(/^[st]__/, "", taxon)
            gsub(/[^A-Za-z0-9._-]/, "_", taxon)
            if (taxon != "") print taxon
        }
    ' "${detected_info_tsv}" | sort -u > detected_safe_taxa.txt

    # ── Step 3: Filter lookup to detected safe_taxa only ────────────────────────
    awk -F'\\t' 'NR==FNR { t[\$1]=1; next } (\$2 in t) { print }' \\
        detected_safe_taxa.txt acc_safe_taxon_full.tsv > acc_safe_taxon_filtered.tsv

    # ── Step 4: Extract reference sequences for detected taxa; write ref_lengths ─
    awk -F'\\t' '{ print \$1 }' acc_safe_taxon_filtered.tsv > ref_accessions.txt
    seqkit grep -f ref_accessions.txt "\${fna}" > species_refs.fasta

    echo -e "Accession\\tLength" > ${meta.id}_ref_lengths.tsv
    awk '/^>/{if(len>0) print name"\\t"len; name=substr(\$0,2); gsub(/ .*/,"",name); len=0} \\
         !/^>/{len+=length(\$0)} \\
         END{if(len>0) print name"\\t"len}' \\
        species_refs.fasta >> ${meta.id}_ref_lengths.tsv

    # ── Write TSV header (output exists even if we exit early) ──────────────────
    echo -e "Sample\\tSpecies\\tScaffold_ID\\tMatched_Reference\\tIdentity_%\\tAlign_Len\\tQuery_Len\\tMismatches\\tGap_Opens\\tQ_Start\\tQ_End\\tS_Start\\tS_End\\tE-value\\tBit_Score\\tCov_%" \\
        > ${meta.id}_blastn.tsv

    # Guard: empty query FASTA → has_contigs all false, skip BLAST
    if [ \$(grep -c "^>" ${query} 2>/dev/null || echo 0) -eq 0 ]; then
        awk -v sample="${meta.id}" '{ print sample"\\t"\$1"\\tfalse" }' \\
            detected_safe_taxa.txt > ${meta.id}_has_contigs.txt
        echo "WARNING: empty query FASTA for ${meta.id} -- skipping BLAST"
        exit 0
    fi

    # Guard: no viruses detected → empty has_contigs, skip BLAST
    if [ ! -s detected_safe_taxa.txt ]; then
        touch ${meta.id}_has_contigs.txt
        echo "WARNING: no viruses detected for ${meta.id} -- skipping BLAST"
        exit 0
    fi

    # ── Step 5: Build BLAST DB from the full .fna ────────────────────────────────
    makeblastdb \\
        -in "\${fna}" \\
        -dbtype nucl \\
        -out full_blastdb \\
        -title "${meta.id}_fulldb"

    # ── Step 6: BLAST whole-sample contigs against the full DB ──────────────────
    # outfmt 6: col1=qseqid col2=sseqid col3=pident col4=length col5=qlen
    #           col6=mismatch col7=gapopen col8=qstart col9=qend col10=sstart
    #           col11=send col12=evalue col13=bitscore col14=qcovs
    blastn \\
        -query ${query} \\
        -db full_blastdb \\
        -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \\
        -evalue 1e-5 \\
        -num_threads ${task.cpus} \\
        -out blast_raw.tsv

    # ── Step 7: Annotate hits with safe_taxon; filter to detected taxa only ──────
    # Appends safe_taxon (col 15) to each blast row; rows whose Matched_Reference
    # is absent from acc_safe_taxon_filtered are dropped (non-detected viruses).
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{t[\$1]=\$2; next} \\
         { if (\$2 in t) print \$0"\\t"t[\$2] }' \\
        acc_safe_taxon_filtered.tsv blast_raw.tsv > blast_with_taxon.tsv

    # ── Step 8: Assign each contig exclusively to its best-match safe_taxon ──────
    # Sum bitscore (col 13) per (contig, safe_taxon); pick taxon with highest sum.
    awk -F'\\t' '
        BEGIN { OFS="\\t" }
        {
            contig = \$1; taxon = \$15; bs = \$13+0
            score[contig, taxon] += bs
        }
        END {
            for (key in score) {
                n = split(key, parts, SUBSEP)
                contig = parts[1]; taxon = parts[2]
                if (!(contig in best_score) || score[key] > best_score[contig]) {
                    best_score[contig] = score[key]
                    best_taxon[contig] = taxon
                }
            }
            for (contig in best_taxon) print contig, best_taxon[contig]
        }
    ' blast_with_taxon.tsv > contig_assignments.tsv

    # ── Step 9: Reconstruct BLAST TSV — keep only hits matching assigned taxon ───
    awk -F'\\t' -v sample="${meta.id}" '
        BEGIN { OFS="\\t" }
        NR==FNR { best[\$1]=\$2; next }
        {
            contig = \$1; taxon = \$15
            if (contig in best && best[contig] == taxon)
                print sample, taxon, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14
        }
    ' contig_assignments.tsv blast_with_taxon.tsv >> ${meta.id}_blastn.tsv

    # ── Step 10: Write has_contigs.txt — one row per detected safe_taxon ─────────
    awk -F'\\t' '{ print \$2 }' contig_assignments.tsv | sort -u > assigned_taxa.txt
    awk -F'\\t' -v sample="${meta.id}" '
        NR==FNR { assigned[\$1]=1; next }
        { taxon=\$1; print sample"\\t"taxon"\\t"(taxon in assigned ? "true" : "false") }
    ' assigned_taxa.txt detected_safe_taxa.txt > ${meta.id}_has_contigs.txt
    """
}
