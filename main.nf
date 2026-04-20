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
include { SPADES_ASSEMBLY        } from './modules/local/spades_assembly'
include { BLASTN_VALIDATE        } from './modules/local/blastn_validate'
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
    ║  spades_cap  : ${params.validate_spades_max_pairs ? "${params.validate_spades_max_pairs} pairs" : 'disabled'}
    ║  spades_mem  : 220.GB
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

    Channel
        .fromPath("${params.esviritu_db}/*.tsv", checkIfExists: true)
        .first()
        .set { ch_esviritu_db_meta }

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
    // Enabled with --validate. For each sample:
    //   1. De novo assemble all quality-filtered reads with SPAdes.
    //      Contigs < validate_min_contig_len are filtered; if total read count
    //      < validate_min_reads the assembly is skipped (empty query FASTA).
    //   2. BLAST assembled contigs against the full EsViritu .fna.
    //      Each contig is assigned exclusively to the DB taxon with the
    //      highest total BLAST bitscore.
    //   3. Taxonomy grouping (subspecies / species) and safe_taxon derivation
    //      are performed inside BLASTN_VALIDATE using params.assembly_taxon_level.
    if (params.validate) {

        // Assemble all quality-filtered reads per sample in one SPAdes job.
        SPADES_ASSEMBLY(FASTP_DEDUP.out.reads, ch_spades_mode)

        SPADES_ASSEMBLY.out.query
            .set { ch_blast_input }

        BLASTN_VALIDATE(ch_blast_input, ch_esviritu_db)

        // Collect per-sample sidecar files produced by BLASTN_VALIDATE
        BLASTN_VALIDATE.out.has_contigs
            .collect()
            .set { ch_has_contigs_all }

        BLASTN_VALIDATE.out.ref_lengths
            .collect()
            .set { ch_ref_lengths_all }

        VISUALIZE_VALIDATION(BLASTN_VALIDATE.out.blast_results, ch_ref_lengths_all)

        // ── Comprehensive overview table (with Part 3 BLAST data) ──────────
        MAKE_OVERVIEW_TABLE(
            SUMMARIZE_READ_STATS.out.read_stats,
            SUMMARIZE_ESV.out.assembly_summary_tsv,
            SUMMARIZE_READ_STATS.out.info_enriched,
            BLASTN_VALIDATE.out.blast_results.map { meta, tsv -> tsv }.collect(),
            ch_has_contigs_all,
            ch_ref_lengths_all,
            ch_esviritu_db_meta,
            true,
            params.validate_min_reads,
            params.assembly_taxon_level
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
            ch_esviritu_db_meta,
            false,
            params.validate_min_reads,
            params.assembly_taxon_level
        )
    }
}
