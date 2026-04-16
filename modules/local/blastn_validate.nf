/*
 * BLASTN_VALIDATE
 *
 * Assembles all quality-filtered reads per sample with SPAdes (upstream), then:
 *   1. Builds an accession -> safe_taxon lookup from the full EsViritu metadata TSV
 *      using the configured assembly_taxon_level (subspecies / species).
 *   2. Builds a BLAST nucleotide DB from the full EsViritu .fna.
 *   3. BLASTs the whole-sample contigs against that DB.
 *   4. Assigns each contig exclusively to the DB taxon that accumulates the
 *      highest total BLAST bitscore across all its hits.
 *   5. Writes per-sample output files compatible with the downstream R scripts:
 *      - {id}_blastn.tsv      : hits with Sample and Species (safe_taxon) columns
 *      - {id}_ref_lengths.tsv : lengths of reference sequences retained in final hits
 *      - {id}_has_contigs.txt : sample assembly sentinel plus one row per BLAST taxon
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
    tuple val(meta), path(query)
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

    # Step 1: Build full accession -> safe_taxon lookup from the metadata TSV.
    awk -F'\t' -v level="${taxon_level}" '
        BEGIN { OFS="\t" }
        NR==1 { next }
        {
            acc        = \$1
            species    = \$11
            subspecies = \$12
            taxon = (level == "subspecies" && subspecies != "" && subspecies != "NA") \
                    ? subspecies : species
            sub(/^[st]__/, "", taxon)
            gsub(/[^A-Za-z0-9._-]/, "_", taxon)
            if (taxon != "") print acc, taxon
        }
    ' "\${meta_tsv}" > acc_safe_taxon_full.tsv

    echo -e "Sample\tSpecies\tScaffold_ID\tMatched_Reference\tIdentity_%\tAlign_Len\tQuery_Len\tMismatches\tGap_Opens\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score\tCov_%" \
        > ${meta.id}_blastn.tsv
    echo -e "Accession\tLength" > ${meta.id}_ref_lengths.tsv

    # Empty query FASTA means SPAdes produced no usable contigs.
    if [ \$(grep -c "^>" ${query} 2>/dev/null || echo 0) -eq 0 ]; then
        echo -e "${meta.id}\t__sample_assembly__\tfalse" > ${meta.id}_has_contigs.txt
        echo "WARNING: empty query FASTA for ${meta.id} -- skipping BLAST"
        exit 0
    fi

    # Step 2: Build BLAST DB from the full .fna.
    makeblastdb \
        -in "\${fna}" \
        -dbtype nucl \
        -out full_blastdb \
        -title "${meta.id}_fulldb"

    # Step 3: BLAST whole-sample contigs against the full DB.
    blastn \
        -query ${query} \
        -db full_blastdb \
        -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
        -evalue 1e-5 \
        -num_threads ${task.cpus} \
        -out blast_raw.tsv

    # Step 4: Annotate hits with safe_taxon from the full metadata lookup.
    awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{t[\$1]=\$2; next} \
         { if (\$2 in t) print \$0"\t"t[\$2] }' \
        acc_safe_taxon_full.tsv blast_raw.tsv > blast_with_taxon.tsv

    # Step 5: Assign each contig exclusively to its best-match safe_taxon.
    awk -F'\t' '
        BEGIN { OFS="\t" }
        {
            contig = \$1; taxon = \$15; bs = \$13+0
            score[contig, taxon] += bs
        }
        END {
            for (key in score) {
                split(key, parts, SUBSEP)
                contig = parts[1]; taxon = parts[2]
                if (!(contig in best_score) || score[key] > best_score[contig]) {
                    best_score[contig] = score[key]
                    best_taxon[contig] = taxon
                }
            }
            for (contig in best_taxon) print contig, best_taxon[contig]
        }
    ' blast_with_taxon.tsv > contig_assignments.tsv

    # Step 6: Keep only hits matching the assigned taxon for each contig.
    awk -F'\t' -v sample="${meta.id}" '
        BEGIN { OFS="\t" }
        NR==FNR { best[\$1]=\$2; next }
        {
            contig = \$1; taxon = \$15
            if (contig in best && best[contig] == taxon)
                print sample, taxon, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14
        }
    ' contig_assignments.tsv blast_with_taxon.tsv >> ${meta.id}_blastn.tsv

    # Step 7: Capture reference lengths for accessions retained in final hits.
    tail -n +2 ${meta.id}_blastn.tsv | cut -f4 | sort -u > ref_accessions.txt

    if [ -s ref_accessions.txt ]; then
        seqkit grep -f ref_accessions.txt "\${fna}" > species_refs.fasta
        awk '/^>/{if(len>0) print name"\t"len; name=substr(\$0,2); gsub(/ .*/,"",name); len=0} \
             !/^>/{len+=length(\$0)} \
             END{if(len>0) print name"\t"len}' \
            species_refs.fasta >> ${meta.id}_ref_lengths.tsv
    fi

    # Step 8: Write sample-level assembly status plus BLAST-assigned taxa.
    {
        echo -e "${meta.id}\t__sample_assembly__\ttrue"
        awk -F'\t' -v sample="${meta.id}" '{ print sample"\t"\$2"\ttrue" }' contig_assignments.tsv | sort -u
    } > ${meta.id}_has_contigs.txt
    """
}
