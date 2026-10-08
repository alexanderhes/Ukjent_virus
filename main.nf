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
include { FILTER_HOST_CONTIGS    } from './modules/local/filter_host_contigs'
include { EXTRACT_HOST_FASTA; MAKE_HOST_BLASTDB } from './modules/local/host_blastdb'
include { BLASTN_VALIDATE        } from './modules/local/blastn_validate'
include { VISUALIZE_VALIDATION   } from './modules/local/visualize_validation'
include { SCAFFOLD_GROUPS        } from './modules/local/scaffold_groups'
include { SCAFFOLD_MINIMAP2      } from './modules/local/scaffold_minimap2'
include { POLISH_SCAFFOLD        } from './modules/local/polish_scaffold'

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
    ║  scaffold    : ${params.validate && params.scaffold ? "minimap2 (polish: ${params.polish_rounds} rounds, min depth ${params.polish_min_depth})" : 'disabled'}
    ╚═══════════════════════════════════════════════╝
    """.stripIndent()

// ── Provisional verdict thresholds ───────────────────────────────────────────
// Shown at start and end of every run until the verdict_* thresholds are
// validated and params.verdict_thresholds_provisional is set to false.
if (params.verdict_thresholds_provisional) {
    def verdictWarning = """
        WORK IN PROGRESS: overview verdict/flag thresholds are provisional and not validated.
        Treat verdict, esv_verdict, blast_verdict, flags, esv_flags and blast_flags as indicative only.
          esv      : min_reads=${params.verdict_esv_min_reads} min_breadth_pct=${params.verdict_esv_min_breadth_pct} divergent_identity_pct=${params.verdict_esv_divergent_identity_pct}
          blast    : contig_min_aln_bp=${params.verdict_contig_min_aln_bp} contig_min_aln_pct=${params.verdict_contig_min_aln_pct} artefact_max_aln_bp=${params.verdict_artefact_max_aln_bp} artefact_max_aln_pct=${params.verdict_artefact_max_aln_pct} divergent_identity_pct=${params.verdict_blast_divergent_identity_pct}
          recurrent: min_samples=${params.verdict_recurrent_min_samples}
        Set verdict_thresholds_provisional = false in conf/params.config once validated.
        """.stripIndent()
    log.warn verdictWarning
    workflow.onComplete { log.warn verdictWarning }
}

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
    // Order: trim -> dedup -> host filter.
    // Trimming runs before host filtering: with adapters still on, short-insert
    // pairs read through into adapter, their soft-clipped mates "dovetail" and
    // bowtie2 does not count them as concordant, so host pairs leak through.
    // Deduplication runs before host filtering so the bowtie2 step (the most
    // compute-intensive pre-processing step) gets fewer reads.
    FASTP_TRIM(ch_reads)
    FASTP_DEDUP(FASTP_TRIM.out.reads)
    HOST_FILTER(FASTP_DEDUP.out.reads, ch_host_index)
    ESVIRITU(HOST_FILTER.out.reads, ch_esviritu_db)

    // ── Read-count funnel (per sample) ─────────────────────────────
    // raw / trimmed from FASTP_TRIM, dedup from FASTP_DEDUP, host-filtered
    // (= analysed) from the HOST_FILTER output count
    FASTP_TRIM.out.json
        .join(FASTP_DEDUP.out.json)
        .join(HOST_FILTER.out.read_count)
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
    //   1b. (--host_contig_filter, default on) Remove host contigs -- mainly
    //      human rRNA that assembles into contigs -- by megablast against the
    //      host BLAST database (--host_blastdb) in FILTER_HOST_CONTIGS.
    //   2. BLAST assembled contigs against the full EsViritu .fna.
    //      Each contig is assigned exclusively to the DB taxon with the
    //      highest total BLAST bitscore.
    //   3. Taxonomy grouping (subspecies / species) and safe_taxon derivation
    //      are performed inside BLASTN_VALIDATE using params.assembly_taxon_level.
    //   4. (--scaffold, default on) Contigs of the same taxon (or segment) are
    //      combined against their best reference with minimap2 + samtools
    //      consensus, then polished with the sample's reads into one genome per
    //      sample x taxon (x segment), published in validation/<sample>/genomes/.
    if (params.validate) {

        // Assemble all quality-filtered reads per sample in one SPAdes job.
        SPADES_ASSEMBLY(HOST_FILTER.out.reads, ch_spades_mode)

        // Remove human rRNA that assembles into contigs (and any other host
        // contigs) before BLAST and scaffolding; see FILTER_HOST_CONTIGS.
        if (params.host_contig_filter) {
            if (!params.host_blastdb) {
                error "--host_contig_filter needs --host_blastdb (host BLAST database prefix, set in the host_<alias> profile)"
            }
            // Build the host BLAST database from the bowtie2 index on first use;
            // later runs find it at params.host_blastdb and skip the build.
            if (file("${params.host_blastdb}.n*")) {
                Channel
                    .fromPath("${params.host_blastdb}.*")
                    .collect()
                    .set { ch_host_blastdb }
            } else {
                log.info "Host BLAST database ${params.host_blastdb} not found -- building it from the bowtie2 host index"
                EXTRACT_HOST_FASTA(ch_host_index)
                MAKE_HOST_BLASTDB(EXTRACT_HOST_FASTA.out.fasta)
                MAKE_HOST_BLASTDB.out.db
                    .collect()
                    .set { ch_host_blastdb }
            }
            FILTER_HOST_CONTIGS(SPADES_ASSEMBLY.out.query, ch_host_blastdb)
            FILTER_HOST_CONTIGS.out.query.set { ch_blast_input }
        } else {
            SPADES_ASSEMBLY.out.query.set { ch_blast_input }
        }

        BLASTN_VALIDATE(ch_blast_input, ch_esviritu_db)

        // Collect per-sample sidecar files produced by BLASTN_VALIDATE
        BLASTN_VALIDATE.out.has_contigs
            .collect()
            .set { ch_has_contigs_all }

        BLASTN_VALIDATE.out.ref_lengths
            .map { meta, tsv -> tsv }
            .collect()
            .set { ch_ref_lengths_all }

        VISUALIZE_VALIDATION(BLASTN_VALIDATE.out.blast_results, ch_ref_lengths_all)

        // ── Combine contigs per taxon into polished genomes ────────────────
        if (params.scaffold) {
            SCAFFOLD_GROUPS(
                BLASTN_VALIDATE.out.blast_results.join(BLASTN_VALIDATE.out.ref_lengths)
            )

            ch_blast_input
                .join(SCAFFOLD_GROUPS.out.groups)
                .set { ch_scaffold_input }

            SCAFFOLD_MINIMAP2(ch_scaffold_input, ch_esviritu_db)

            // One polishing task per sample, on the host-filtered reads
            SCAFFOLD_MINIMAP2.out.templates
                .join(SCAFFOLD_GROUPS.out.contigs)
                .join(HOST_FILTER.out.reads)
                .set { ch_polish_input }

            POLISH_SCAFFOLD(ch_polish_input)

            POLISH_SCAFFOLD.out.stats
                .map { meta, tsv -> tsv }
                .set { ch_scaffold_stats }

            // Batch table of all genomes of the run
            ch_scaffold_stats
                .collectFile(
                    name: 'genome_summary.tsv',
                    keepHeader: true,
                    sort: { it.name },
                    storeDir: "${params.outdir}/validation"
                )

            ch_scaffold_stats
                .collect()
                .ifEmpty([])
                .set { ch_scaffold_stats_all }
        } else {
            Channel.value([]).set { ch_scaffold_stats_all }
        }

        // ── Comprehensive overview table (with Part 3 BLAST data) ──────────
        MAKE_OVERVIEW_TABLE(
            SUMMARIZE_READ_STATS.out.read_stats,
            SUMMARIZE_ESV.out.assembly_summary_tsv,
            SUMMARIZE_READ_STATS.out.info_enriched,
            BLASTN_VALIDATE.out.blast_results.map { meta, tsv -> tsv }.collect(),
            ch_has_contigs_all,
            ch_ref_lengths_all,
            ch_scaffold_stats_all,
            BLASTN_VALIDATE.out.tblastx.collect(),
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
            [],
            [],
            ch_esviritu_db_meta,
            false,
            params.validate_min_reads,
            params.assembly_taxon_level
        )
    }
}
