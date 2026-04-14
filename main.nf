#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HOST_FILTER          } from './modules/local/host_filter'
include { FASTP_TRIM           } from './modules/local/fastp_trim'
include { FASTP_DEDUP          } from './modules/local/fastp_dedup'
include { ESVIRITU             } from './modules/local/esviritu'
include { SUMMARIZE_ESV        } from './modules/local/summarize_esv'
include { COLLECT_READ_STATS   } from './modules/local/collect_read_stats'
include { SUMMARIZE_READ_STATS } from './modules/local/summarize_read_stats'
include { MAKE_OVERVIEW_TABLE  } from './modules/local/make_overview_table'

// Validation sub-workflow modules (only loaded when --validate is enabled)
include { SPLIT_VIRAL_READS      } from './modules/local/split_viral_reads'
include { SPADES_ASSEMBLY        } from './modules/local/spades_assembly'
include { BLASTN_VALIDATE        } from './modules/local/blastn_validate'
include { SUMMARIZE_VALIDATION   } from './modules/local/summarize_validation'
include { VISUALIZE_VALIDATION   } from './modules/local/visualize_validation'

// ── Logging ──────────────────────────────────────────────────────────────────
log.info """
    ╔═══════════════════════════════════════════════╗
    ║           EsViritu Nextflow Pipeline          ║
    ╠═══════════════════════════════════════════════╣
    ║  samplesheet : ${params.samplesheet}
    ║  host_index  : ${params.host_index}
    ║  esviritu_db : ${params.esviritu_db}
    ║  outdir      : ${params.outdir}
    ║  validate    : ${params.validate}
    ║  spades_mode : ${params.spades_mode ?: (params.esviritu_db ==~ /(?i).*HEV.*/ ? 'rnaviral (auto)' : 'meta (auto)')}
    ╚═══════════════════════════════════════════════╝
    """.stripIndent()

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // ── Parse samplesheet ────────────────────────────────────────────────────
    // Expected format: two columns (sample;fastq_dir)
    // sep is configurable via params.samplesheet_sep (default: ';')
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: params.samplesheet_sep, strip: true)
        .map { row ->
            // Validate required columns
            if (!row.containsKey('sample') || !row.containsKey('fastq_dir')) {
                error "Samplesheet must contain 'sample' and 'fastq_dir' columns. Found: ${row.keySet()}"
            }

            def meta     = [id: row.sample.replaceAll(/\s/, '_')]
            def fastq_dir = file(row.fastq_dir, checkIfExists: true)

            // Glob for R1 and R2
            def r1_files = fastq_dir.listFiles().findAll { it.name =~ /(?i).*R1.*\.fastq\.gz$/ }.sort()
            def r2_files = fastq_dir.listFiles().findAll { it.name =~ /(?i).*R2.*\.fastq\.gz$/ }.sort()

            if (r1_files.size() != 1 || r2_files.size() != 1) {
                log.warn "[SKIP] ${meta.id}: expected exactly 1 R1 and 1 R2 in ${fastq_dir} " +
                         "(found ${r1_files.size()} R1, ${r2_files.size()} R2)"
                return null
            }

            return [meta, r1_files[0], r2_files[0]]
        }
        .filter { it != null }
        .set { ch_reads }

    // ── Stage bowtie2 index files ────────────────────────────────────────────
    Channel
        .fromPath("${params.host_index}*.bt2", checkIfExists: true)
        .collect()
        .set { ch_host_index }

    // ── Stage EsViritu database directory ───────────────────────────────────
    // Using Channel.value so the path is broadcast and reusable across multiple
    // processes (ESVIRITU, BLASTN_VALIDATE) without being consumed.
    Channel
        .fromPath(params.esviritu_db, checkIfExists: true)
        .first()
        .set { ch_esviritu_db }

    // ── Derive SPAdes assembly mode ──────────────────────────────────────────
    // Explicit --spades_mode always wins. Otherwise auto-detect from DB name:
    //   DB path/name contains 'HEV' (case-insensitive) → rnaviral
    //   anything else                                   → meta
    // To add future RNA-virus DBs, extend the pattern with | e.g. HEV|RSV
    def effective_spades_mode = params.spades_mode ?:
        (params.esviritu_db ==~ /(?i).*HEV.*/ ? 'rnaviral' : 'meta')
    Channel.value(effective_spades_mode).set { ch_spades_mode }

    // ── Pipeline steps ───────────────────────────────────────────────────────────────────────────────
    HOST_FILTER(ch_reads, ch_host_index)
    FASTP_TRIM(HOST_FILTER.out.reads)
    FASTP_DEDUP(FASTP_TRIM.out.reads)
    ESVIRITU(FASTP_DEDUP.out.reads, ch_esviritu_db)

    // ── Read-count funnel (per sample) ─────────────────────────────
    HOST_FILTER.out.raw_read_count
        .join(FASTP_TRIM.out.json)
        .join(FASTP_DEDUP.out.json)
        .set { ch_read_stats_input }

    COLLECT_READ_STATS(ch_read_stats_input)

    // ── Batch summary ────────────────────────────────────
    // Collect all per-sample TSV outputs then run a single summarise step
    ESVIRITU.out.tsv_files
        .collect()
        .set { ch_esv_all }

    SUMMARIZE_ESV(ch_esv_all)

    // ── Batch read-stats + enriched detection table ───────────────────
    SUMMARIZE_READ_STATS(
        COLLECT_READ_STATS.out.tsv.map { meta, tsv -> tsv }.collect(),
        SUMMARIZE_ESV.out.info_tsv
    )

    // ── Validation sub-workflow ──────────────────────────────────────────────
    // Enabled with --validate. For each detected viral species in every sample:
    //   1. Extract reads from the third-pass BAM that map to any accession in
    //      that species group (multiple strains pooled per species).
    //   2. De novo assemble with SPAdes; filter contigs < validate_min_contig_len.
    //      If read count < validate_min_reads, convert R1 to FASTA instead.
    //   3. BLAST assembled contigs against the full EsViritu nucleotide DB.
    //   4. Summarise per-species BLAST hits into a per-sample TSV.
    if (params.validate) {

        // Fan out ESVIRITU results to one channel item per (sample × species):
        //   parse detected_virus.info.tsv, group all accessions by species.
        //   Accessions are stored in meta so they are available in BLASTN_VALIDATE.
        ESVIRITU.out.tsv_info_meta
            .splitCsv(header: true, sep: '\t', elem: 1)
            .map  { meta, row ->
                // Resolve grouping taxon: subspecies when available and requested,
                // otherwise fall back to species.
                def taxon = (params.assembly_taxon_level == 'subspecies' &&
                             row.subspecies && row.subspecies.trim() != '')
                             ? row.subspecies
                             : row.species
                [meta, taxon, row.Accession]
            }
            .groupTuple(by: [0, 1])
            // [meta, taxon_str, [acc1, acc2, ...]]
            .map  { meta, taxon, accs ->
                // Strip rank prefix (s__ / t__) and replace all non-alphanumeric
                // characters (incl. brackets used in Rotavirus serotype names)
                // with underscores for safe use in file names.
                // IMPORTANT: this regex MUST stay in sync with to_safe() in
                // bin/make_overview_table.R — any divergence will silently break
                // safe_taxon joins in the overview table.
                def safe = taxon.replaceAll(/^[st]__/, '').replaceAll(/[^A-Za-z0-9._-]/, '_')
                [meta + [species: taxon, safe_species: safe, accessions: accs], accs]
            }
            .set { ch_species_accs }

        // Key the species channel and BAM channel both on meta.id so we can
        // join them despite meta having different extra fields.
        ch_species_accs
            .map { meta, accs -> [meta.id, meta, accs] }
            .combine(
                ESVIRITU.out.third_bam
                    .join(ESVIRITU.out.third_bai)
                    .map { meta, bam, bai -> [meta.id, bam, bai] },
                by: 0
            )
            .map { id, meta, accs, bam, bai -> [meta, accs, bam, bai] }
            .set { ch_split_input }

        SPLIT_VIRAL_READS(ch_split_input)

        SPADES_ASSEMBLY(SPLIT_VIRAL_READS.out.reads, ch_spades_mode)

        BLASTN_VALIDATE(SPADES_ASSEMBLY.out.query, ch_esviritu_db)

        // Collect per-species sidecar files produced by BLASTN_VALIDATE
        BLASTN_VALIDATE.out.has_contigs
            .collect()
            .set { ch_has_contigs_all }

        BLASTN_VALIDATE.out.ref_lengths
            .collect()
            .set { ch_ref_lengths_all }

        // Group all per-species BLAST TSVs for each sample, then summarise
        BLASTN_VALIDATE.out.blast_results
            .map  { meta, tsv -> [meta.id, meta, tsv] }
            .groupTuple(by: 0)
            .map  { id, metas, tsvs -> [metas[0].subMap('id'), tsvs] }
            .set  { ch_summary_input }

        SUMMARIZE_VALIDATION(ch_summary_input)

        VISUALIZE_VALIDATION(SUMMARIZE_VALIDATION.out.summary, ch_ref_lengths_all)

        // ── Comprehensive overview table (with Part 3 BLAST data) ──────────
        MAKE_OVERVIEW_TABLE(
            SUMMARIZE_READ_STATS.out.read_stats,
            SUMMARIZE_ESV.out.assembly_summary_tsv,
            SUMMARIZE_READ_STATS.out.info_enriched,
            SUMMARIZE_VALIDATION.out.summary.map { meta, tsv -> tsv }.collect(),
            ch_has_contigs_all,
            ch_ref_lengths_all,
            true,
            params.validate_min_reads
        )
    } else {
        // ── Overview table without validation (Part 3 columns = NA) ────────
        Channel
            .of("# no validation run")
            .collectFile(name: "no_validation.txt")
            .set { ch_no_val }

        Channel
            .of("# no validation run")
            .collectFile(name: "no_has_contigs.txt")
            .set { ch_no_has_contigs }

        Channel
            .of("# no validation run")
            .collectFile(name: "no_ref_lengths.txt")
            .set { ch_no_ref_lengths }

        MAKE_OVERVIEW_TABLE(
            SUMMARIZE_READ_STATS.out.read_stats,
            SUMMARIZE_ESV.out.assembly_summary_tsv,
            SUMMARIZE_READ_STATS.out.info_enriched,
            ch_no_val,
            ch_no_has_contigs,
            ch_no_ref_lengths,
            false,
            params.validate_min_reads
        )
    }
}
