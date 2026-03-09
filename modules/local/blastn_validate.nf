/*
 * BLASTN_VALIDATE
 *
 * For each detected viral species:
 *   1. Use the EsViritu metadata TSV to find ALL accessions for that species
 *      in the database (not just the detected ones — gives fuller coverage).
 *      Uses grep -F to match the exact species string, cut -f1 for accession.
 *   2. Extract those reference sequences from the .fna with seqkit grep.
 *   3. Build a small per-species BLAST nucleotide DB with makeblastdb.
 *      (Small input avoids the segfault seen with the full ~20k-entry .fna.)
 *   4. BLAST the assembled contigs (or raw reads) against that per-species DB.
 *
 * Output format columns: qseqid sseqid stitle pident length qcovs evalue bitscore
 */

process BLASTN_VALIDATE {
    tag "${meta.id}:${meta.safe_species}"
    label 'process_medium'

    publishDir "${params.outdir}/validation/${meta.id}/blast", mode: 'copy'

    input:
    tuple val(meta), path(query)
    path(db_dir)    // staged EsViritu DB directory (contains .fna and metadata TSV)

    output:
    tuple val(meta), path("${meta.id}_${meta.safe_species}_blastn.tsv"),    emit: blast_results
    path "${meta.id}_${meta.safe_species}_ref_lengths.tsv",                 emit: ref_lengths
    path "${meta.id}_${meta.safe_species}_has_contigs.txt",                 emit: has_contigs

    script:
    def species_val = meta.species
    """
    fna=\$(find -L ${db_dir} -name "*.fna"      | head -1)
    meta_tsv=\$(find -L ${db_dir} -name "*.tsv" | head -1)

    # Step 1: collect all accessions for this species from the metadata TSV.
    # Use grep -F to match the exact species string, then cut the accession (col 1).
    grep -F "${species_val}" "\${meta_tsv}" | cut -f1 > species_accessions.txt

    # Guard: if no accessions found in DB, write empty result and exit cleanly
    if [ ! -s species_accessions.txt ]; then
        echo -e "Sample\tSpecies\tScaffold_ID\tMatched_Reference\tIdentity_%\tAlign_Len\tQuery_Len\tMismatches\tGap_Opens\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score\tCov_%" \\
            > ${meta.id}_${meta.safe_species}_blastn.tsv
        echo -e "Accession\tLength" > ${meta.id}_${meta.safe_species}_ref_lengths.tsv
        printf '%s\t%s\ttrue\n' "${meta.id}" "${meta.safe_species}" \\
            > ${meta.id}_${meta.safe_species}_has_contigs.txt
        echo "WARNING: no DB entries found for species '${species_val}' -- skipping BLAST"
        exit 0
    fi

    # Step 2: extract all matching reference sequences from the .fna by accession ID.
    seqkit grep -f species_accessions.txt "\${fna}" > species_refs.fasta

    # Guard: if extraction yielded nothing, write empty result and exit cleanly
    if [ ! -s species_refs.fasta ]; then
        echo -e "Sample\tSpecies\tScaffold_ID\tMatched_Reference\tIdentity_%\tAlign_Len\tQuery_Len\tMismatches\tGap_Opens\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score\tCov_%" \\
            > ${meta.id}_${meta.safe_species}_blastn.tsv
        echo -e "Accession\tLength" > ${meta.id}_${meta.safe_species}_ref_lengths.tsv
        printf '%s\t%s\ttrue\n' "${meta.id}" "${meta.safe_species}" \\
            > ${meta.id}_${meta.safe_species}_has_contigs.txt
        echo "WARNING: no sequences extracted for species '${species_val}' -- skipping BLAST"
        exit 0
    fi

    # Reference sequence lengths (used to set accurate x-axis limits in contig plots
    # and as denominators for blast_genome_cov_pct).  Written here so the file is
    # guaranteed to exist before any further exit.
    echo -e "Accession\tLength" > ${meta.id}_${meta.safe_species}_ref_lengths.tsv
    awk '/^>/{if(len>0)print name"\\t"len; name=substr(\$0,2); gsub(/ .*/,"",name); len=0} !/^>/{len+=length(\$0)} END{if(len>0)print name"\\t"len}' \\
        species_refs.fasta >> ${meta.id}_${meta.safe_species}_ref_lengths.tsv

    # Step 3: build a small per-species BLAST DB
    makeblastdb \\
        -in species_refs.fasta \\
        -dbtype nucl \\
        -out species_blastdb \\
        -title "${meta.safe_species}"

    # Guard: if the query FASTA has no sequences (SPAdes produced no contigs or
    # read count was below threshold), write header-only result and exit cleanly.
    # Mark has_contigs=false so the overview table can distinguish this from
    # "contigs assembled but no BLAST hits".
    if [ \$(grep -c "^>" ${query} 2>/dev/null || echo 0) -eq 0 ]; then
        echo -e "Sample\tSpecies\tScaffold_ID\tMatched_Reference\tIdentity_%\tAlign_Len\tQuery_Len\tMismatches\tGap_Opens\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score\tCov_%" \\
            > ${meta.id}_${meta.safe_species}_blastn.tsv
        printf '%s\t%s\tfalse\n' "${meta.id}" "${meta.safe_species}" \\
            > ${meta.id}_${meta.safe_species}_has_contigs.txt
        echo "WARNING: empty query FASTA for species '${species_val}' -- skipping BLAST"
        exit 0
    fi

    # Contigs present — mark before BLAST so the file exists even if blastn fails
    printf '%s\t%s\ttrue\n' "${meta.id}" "${meta.safe_species}" \\
        > ${meta.id}_${meta.safe_species}_has_contigs.txt

    # Step 4: BLAST contigs against the per-species DB
    echo -e "Sample\tSpecies\tScaffold_ID\tMatched_Reference\tIdentity_%\tAlign_Len\tQuery_Len\tMismatches\tGap_Opens\tQ_Start\tQ_End\tS_Start\tS_End\tE-value\tBit_Score\tCov_%" \\
        > ${meta.id}_${meta.safe_species}_blastn.tsv

    blastn \\
        -query ${query} \\
        -db species_blastdb \\
        -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \\
        -evalue 1e-5 \\
        -num_threads ${task.cpus} \\
        -out blast_raw.tsv

    # Prepend Sample and Species columns to every hit row
    awk -v sample="${meta.id}" -v species="${meta.safe_species}" \\
        '{print sample "\t" species "\t" \$0}' blast_raw.tsv \\
        >> ${meta.id}_${meta.safe_species}_blastn.tsv
    """
}
