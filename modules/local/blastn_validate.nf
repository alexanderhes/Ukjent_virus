/*
 * BLASTN_VALIDATE
 *
 * Assembles all quality-filtered reads per sample with SPAdes (upstream), then:
 *   1. Builds an accession -> safe_taxon lookup from the full EsViritu metadata TSV
 *      using the configured assembly_taxon_level (subspecies / species).
 *   2. Builds a BLAST nucleotide DB from the full EsViritu .fna.
 *   3. BLASTs the whole-sample contigs against that DB.
 *   4. Assigns each contig exclusively to the DB accession that accumulates the
 *      highest total BLAST bitscore across all its hits.
 *   5. Maps the winning accession to taxonomy using the metadata TSV.
 *   6. Writes per-sample output files compatible with the downstream R scripts:
 *      - {id}_blastn.tsv      : hits with Sample and Species (safe_taxon) columns
 *      - {id}_ref_lengths.tsv : length and raw segment label of reference sequences
 *                               retained in final hits (Accession, Length, Segment)
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
    # At subspecies level the label is "species subspecies", because bare
    # subspecies labels (e.g. serotype "1") are shared across unrelated species.
    # Must stay in sync with pick_taxon_label() in bin/make_overview_table.R.
    # Also writes accession -> raw segment label (empty when not segmented).
    awk -F'\t' -v level="${taxon_level}" '
        BEGIN { OFS="\t" }
        NR==1 { next }
        {
            acc        = \$1
            segment    = \$4
            species    = \$11
            subspecies = \$12
            gsub(/"/, "", segment)
            if (segment == "NA") segment = ""
            print acc, segment > "acc_segment.tsv"
            sub(/^[a-z]__/, "", species)
            sub(/^[a-z]__/, "", subspecies)
            taxon = species
            if (level == "subspecies" && subspecies != "" && subspecies != "NA")
                taxon = (species != "" && index(subspecies, species)) ? subspecies : species " " subspecies
            gsub(/[^A-Za-z0-9._-]/, "_", taxon)
            if (taxon != "") print acc, taxon
        }
    ' "\${meta_tsv}" > acc_safe_taxon_full.tsv

    echo -e "Sample\tSpecies\tScaffold_ID\tMatched_Reference\tIdentity_%\tAlign_Len\tQuery_Len\tMismatches\tGap_Opens\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score\tCov_%" \
        > ${meta.id}_blastn.tsv
    echo -e "Accession\tLength\tSegment" > ${meta.id}_ref_lengths.tsv

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

    # Step 4: Assign each contig exclusively to its best-match accession.
    awk -F'\t' '
        BEGIN { OFS="\t" }
        {
            contig = \$1; accession = \$2; bs = \$13+0
            score[contig, accession] += bs
        }
        END {
            for (key in score) {
                split(key, parts, SUBSEP)
                contig = parts[1]; accession = parts[2]
                if (!(contig in best_score) || score[key] > best_score[contig]) {
                    best_score[contig] = score[key]
                    best_accession[contig] = accession
                }
            }
            for (contig in best_accession) print contig, best_accession[contig]
        }
    ' blast_raw.tsv > contig_best_accessions.tsv

    # Step 5: Keep only hits matching the assigned accession for each contig,
    # then annotate those hits with the winning accession's taxonomy.
    awk -F'\t' -v sample="${meta.id}" '
        BEGIN { OFS="\t" }
        FILENAME == ARGV[1] { best[\$1]=\$2; next }
        FILENAME == ARGV[2] { taxon[\$1]=\$2; next }
        {
            contig = \$1; accession = \$2
            if (contig in best && best[contig] == accession && accession in taxon)
                print sample, taxon[accession], \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14
        }
    ' contig_best_accessions.tsv acc_safe_taxon_full.tsv blast_raw.tsv >> ${meta.id}_blastn.tsv

    # Step 6: Capture reference lengths and segment labels for accessions retained in final hits.
    tail -n +2 ${meta.id}_blastn.tsv | cut -f4 | sort -u > ref_accessions.txt

    if [ -s ref_accessions.txt ]; then
        seqkit grep -f ref_accessions.txt "\${fna}" > species_refs.fasta
        awk '/^>/{if(len>0) print name"\t"len; name=substr(\$0,2); gsub(/ .*/,"",name); len=0} \
             !/^>/{len+=length(\$0)} \
             END{if(len>0) print name"\t"len}' \
            species_refs.fasta \
        | awk -F'\t' 'BEGIN { OFS="\t" } \
                      FNR==NR { seg[\$1]=\$2; next } \
                      { print \$1, \$2, ((\$1 in seg) ? seg[\$1] : "") }' \
            acc_segment.tsv - >> ${meta.id}_ref_lengths.tsv
    fi

    # Step 7: Write sample-level assembly status plus BLAST-assigned taxa.
    {
        echo -e "${meta.id}\t__sample_assembly__\ttrue"
        tail -n +2 ${meta.id}_blastn.tsv | cut -f2 | sort -u | awk -v sample="${meta.id}" '{ print sample"\t"\$1"\ttrue" }'
    } > ${meta.id}_has_contigs.txt
    """
}
