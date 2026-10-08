# EsViritu Nextflow Pipeline — User Manual

## Table of Contents

1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Installation](#3-installation)
4. [Quick Start](#4-quick-start)
5. [Samplesheet Format](#5-samplesheet-format)
6. [Parameters Reference](#6-parameters-reference)
7. [Pipeline Steps](#7-pipeline-steps)
8. [Output Structure](#8-output-structure)
9. [Output Files Reference](#9-output-files-reference)
10. [Validation Sub-workflow](#10-validation-sub-workflow)
11. [Overview Table Columns](#11-overview-table-columns)
12. [Automated Production Wrapper](#12-automated-production-wrapper)
13. [Resource Configuration](#13-resource-configuration)
14. [Troubleshooting](#14-troubleshooting)

---

## 1. Overview

This pipeline wraps the [EsViritu](https://github.com/dsamoht/esviritu) virus detection tool in a fully automated, reproducible Nextflow workflow. It handles:

- Host-read removal (human T2T + PhiX)
- Quality trimming and de-duplication
- Virus detection using EsViritu
- Read-count funnel tracking across pre-processing steps
- A comprehensive per-sample overview table (read counts, coverage metrics, normalised abundances)
- An optional **validation sub-workflow**: de novo assembly of virus-specific reads and BLAST confirmation of detected species

All steps run inside Docker containers — no conda environments or manual tool installation are required beyond Nextflow and Docker.

---

## 2. Requirements

| Dependency | Minimum version | Notes |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | 23.04 | Run inside the `NEXTFLOW` conda environment on the production server |
| Docker | 20.10 | Must be accessible to the user running Nextflow |
| Java | 11 | Required by Nextflow |

Disk space: ~10–20 GB per sample (intermediates are kept in `work/`; run `nextflow clean -f` after a successful run to free space).

---

## 3. Installation

### From GitHub (production)

The pipeline is pulled automatically by `NGS_wrapper.sh` from GitHub:

```bash
nextflow pull alexanderhes/Ukjent_virus -r main
```

A custom Docker image for the EsViritu process must be built once (or after a `docker/Dockerfile` change):

```bash
docker build -t esviritu_pipeline:latest docker/
```

The wrapper script handles this automatically. A VS Code task is also provided:
**Tasks → Build EsViritu Docker image**

### Local development

Clone the repository and work directly from the local directory:

```bash
git clone https://github.com/alexanderhes/Ukjent_virus.git
cd Ukjent_virus
```

---

## 4. Quick Start

### Minimal run (detection only)

```bash
nextflow run main.nf \
    --samplesheet assets/test_samplesheet.csv \
    --host_index  assets/host_ref/host_index \
    --esviritu_db assets/db/esviritu_DB/v3.2.4 \
    --outdir      results/my_run
```

### With assembly validation

```bash
nextflow run main.nf \
    --samplesheet assets/test_samplesheet.csv \
    --host_index  assets/host_ref/host_index \
    --host_blastdb assets/host_ref/blastdb/host_ref \
    --esviritu_db assets/db/esviritu_DB/v3.2.4 \
    --outdir      results/my_run \
    --validate
```

### Resume a failed/interrupted run

Add `-resume` to reuse cached results for completed steps:

```bash
nextflow run main.nf ... -resume
```

### Run the test dataset

```bash
bash EsViritu_test.sh
```

---

## 5. Samplesheet Format

The samplesheet is a delimited text file (default delimiter: `;`) with a header row containing at minimum `sample` and `fastq_dir` columns.

```
sample;fastq_dir
sample1-UV;/path/to/raw_data/sample1-UV
sample2-UV;/path/to/raw_data/sample2-UV
```

**Rules:**
- `sample` — unique identifier; must not contain spaces. Whitespace is automatically replaced with underscores.
- `fastq_dir` — path to a directory containing **exactly one** R1 and one R2 FASTQ file matching the pattern `*R1*.fastq.gz` / `*R2*.fastq.gz` (case-insensitive). Samples with zero or multiple matches are skipped with a warning.
- The column delimiter can be changed with `--samplesheet_sep` (`;`, `,` or `\t`).

---

## 6. Parameters Reference

### Required

| Parameter | Description |
|---|---|
| `--samplesheet` | Path to the samplesheet file |
| `--host_index` | Path to the bowtie2 index prefix (without `.bt2` extension) |
| `--esviritu_db` | Path to the EsViritu virus database directory |

### Host filtering

| Parameter | Default | Description |
|---|---|---|
| `--sensitive_host_filter` | `false` | Use high-sensitivity unpaired mode (`--very-sensitive-local`, R1+R2 separately, seqkit pair re-sync). Default `false` uses faster paired-end mode (`--sensitive-local`). |

### Output

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `./results` | Root output directory |

### Samplesheet

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet_sep` | `;` | Column delimiter |

### Read pre-processing

| Parameter | Default | Description |
|---|---|---|
| `--min_length` | `50` | Minimum read length after trimming (fastp) |
| `--quality_cutoff` | `20` | Phred quality threshold for trimming (fastp) |
| `--complexity_threshold` | `30` | Low-complexity filter threshold, 0–100 (fastp) |

### Validation sub-workflow

| Parameter | Default | Description |
|---|---|---|
| `--validate` | `false` | Enable the validation sub-workflow |
| `--validate_min_reads` | `50` | Minimum virus-specific reads required to attempt SPAdes assembly |
| `--validate_min_contig_len` | `800` | Minimum contig length (bp) to pass to BLAST |
| `--validate_spades_max_pairs` | `500000` | If greater than 0, sample-level deduplicated read pairs above this cap are deterministically subsampled before SPAdes |
| `--validate_spades_sample_seed` | `11` | Random seed used for deterministic SPAdes input subsampling |
| `--assembly_taxon_level` | `subspecies` | Taxonomic grouping level for read extraction before assembly. `subspecies` uses the finest available rank (recommended for diverse groups such as Enteroviruses and Rotaviruses); `species` always groups at species level. |
| `--tblastx_rescue` | `true` | Experimental: tblastx of contigs without nucleotide support; passing hits are added to the overview |
| `--host_contig_filter` | `true` | Remove host contigs (mainly human rRNA that assembles into contigs) before BLAST |
| `--host_blastdb` | set by `host_<alias>` profile | Host BLAST database prefix; built automatically from the bowtie2 host index on first use if missing |
| `--host_contig_max_pct` | `50` | Remove a contig when at least this % of its length is covered by host BLAST hits |
| `--host_contig_min_identity` | `90` | Minimum identity (%) of host BLAST hits counted towards that coverage |
| `--scaffold` | `true` | Combine the contigs of each taxon (or segment) into one read-polished genome per sample (runs only with `--validate`) |
| `--genome_complete_pct` | `90` | `genome_assessment` = `complete` when at least this % of the reference is called |
| `--genome_partial_pct` | `50` | `genome_assessment` = `partial` when at least this %; below = `fragmented` |
| `--polish_rounds` | `2` | Read-mapping polish iterations |
| `--polish_min_depth` | `5` | Minimum read depth to call a base; positions below become `N` |

### Overview verdict thresholds

Thresholds for the `verdict`, `esv_verdict` and `blast_verdict` columns and the flag columns (see [Section 11](#verdict)).

| Parameter | Default | Description |
|---|---|---|
| `--verdict_esv_min_reads` | `10` | Mapping evidence: minimum EsViritu read count for `esv_verdict = supported` |
| `--verdict_esv_min_breadth_pct` | `5` | Mapping evidence: minimum EsViritu breadth of coverage (%) for `esv_verdict = supported` |
| `--verdict_contig_min_aln_bp` | `300` | De novo evidence: minimum aligned bases on the best contig for `blast_verdict = confirmed` |
| `--verdict_contig_min_aln_pct` | `50` | De novo evidence: minimum % of the best contig aligned for `blast_verdict = confirmed` |
| `--verdict_artefact_max_aln_bp` | `200` | De novo evidence: best contig with fewer aligned bases than this is `artefact_suspected` |
| `--verdict_artefact_max_aln_pct` | `20` | De novo evidence: best contig with a smaller aligned fraction (%) than this is `artefact_suspected` |
| `--verdict_esv_divergent_identity_pct` | `90` | Flag only: `esv_flags` gets `divergent` when EsViritu read identity (%) is below this |
| `--verdict_blast_divergent_identity_pct` | `90` | Flag only: `blast_flags` gets `divergent` when BLAST identity (%) is below this |
| `--verdict_recurrent_min_samples` | `2` | Flag only: `flags` gets `in_<k>/<N>_samples` when the taxon is reported in at least this many samples |

The thresholds above also drive the flags: `esv_flags` `low_reads` / `low_breadth` use the two `verdict_esv_*` minimums, and `blast_flags` `short_alignment` / `low_contig_cov` use the two `verdict_artefact_*` maximums.

### Resources

| Parameter | Default | Description |
|---|---|---|
| `--max_cpus` | `8` | Maximum CPUs per process |
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_time` | `72.h` | Maximum wall time per process |

---

## 7. Pipeline Steps

```
Raw reads (R1 + R2)
       │
       ▼
 FASTP_TRIM           adapter auto-detection, quality/length/complexity filter
       │               poly-G/X trimming
       ▼
 FASTP_DEDUP          deduplication on trimmed reads
       │               (~40–50% of trimmed reads are duplicates)
       ▼
 HOST_FILTER          bowtie2 --sensitive-local (paired-end mode, default)
       │               removes human T2T + PhiX reads from trimmed, deduplicated reads
       ▼
 ESVIRITU             virus detection 
       │               produces per-sample HTML report, BAM files, detection TSVs
       │
       ├──────────────────────────────────────────────────────────────┐
       │  (always)                                                    │  (--validate)
       ▼                                                              ▼
 COLLECT_READ_STATS   parses fastp JSONs & raw counts    SPLIT_VIRAL_READS     extract per-species reads
 SUMMARIZE_READ_STATS combine across all samples         SPADES_ASSEMBLY       de novo assembly (metaSPAdes)
 SUMMARIZE_ESV        batch detection summary            BLASTN_VALIDATE       BLAST contigs vs species DB
 MAKE_OVERVIEW_TABLE  37-column overview TSV             SUMMARIZE_VALIDATION  per-sample BLAST summary
                                                         VISUALIZE_VALIDATION  per-sample contig PDF
```

### Host filtering strategy

Each host has two references, both set in the `host_<alias>` profile in `nextflow.config`: a bowtie2 index (`--host_index`, read-level `HOST_FILTER`) and a BLAST database (`--host_blastdb`, contig-level `FILTER_HOST_CONTIGS` in the validation sub-workflow).

The BLAST database does not have to be built by hand. If no database exists at `--host_blastdb` when a `--validate` run starts, the pipeline reconstructs the host sequences from the bowtie2 index (`EXTRACT_HOST_FASTA`, `bowtie2-inspect`) and builds the database there (`MAKE_HOST_BLASTDB`, `makeblastdb`, written via `storeDir`). Later runs find it and skip the build. The user running the pipeline needs write access to that directory. To build it manually instead:

```bash
makeblastdb -in host_ref.fasta -dbtype nucl -out <host_dir>/blastdb/host_ref
```

By default, bowtie2 is run in **paired-end mode** with `--sensitive-local --no-discordant --no-mixed`. Unmapped pairs are written directly via `--un-conc-gz`, producing synchronised R1/R2 output in a single pass.

Read order is `FASTP_TRIM` → `FASTP_DEDUP` → `HOST_FILTER`. Deduplicating before host filtering means bowtie2, the most compute-intensive pre-processing step, aligns fewer reads. Host filtering runs **after** trimming for a correctness reason: on untrimmed short-insert libraries (e.g. 2×301 bp reads from ~140 bp fragments) both mates read into adapter; in local mode their soft-clipped alignments extend past each other ("dovetail"), which bowtie2 does not count as concordant, so such human pairs were kept as non-host. In one test sample, re-filtering the old host-filtered reads with dovetailing allowed removed a further 78% as human. Trimming first removes the adapter read-through.

To enable the original high-sensitivity mode (independently filtering R1 and R2 with `--very-sensitive-local`, then re-synchronising with `seqkit pair`), pass `--sensitive_host_filter true` to Nextflow or use `--sensitive-filter` in the wrapper script. This is slower but retains reads where only one mate maps to host.

### Trimming and deduplication (two-pass fastp)

Deduplication is intentionally separated from trimming. Running dedup **after** trimming exposes the full duplicate population (typically 40–50%) that is hidden at the raw-read stage by adapter and quality differences between otherwise identical sequences.

---

## 8. Output Structure

```
results/
└── <analysis_name>/
    ├── pipeline_info/
    │   ├── report.html          # Nextflow execution report
    │   ├── timeline.html        # Task timeline
    │   └── trace.txt            # Resource usage per task
    ├── host_filtered/
    │   └── <sample>/            # Host-depleted R1/R2 FASTQs + bowtie2 log
    ├── fastp_trim/
    │   └── <sample>/            # Trimmed R1/R2 + fastp JSON/HTML/log
    ├── fastp_dedup/
    │   └── <sample>/            # Deduplicated R1/R2 + fastp JSON/HTML/log
    ├── esviritu/
    │   └── <sample>/            # Full EsViritu output (HTML report, BAMs, TSVs)
    ├── esviritu_batch/
    │   └── esv_summary/         # Batch summary TSVs across all samples
    ├── overview/
    │   ├── <sample>_overview.tsv          # Per-sample 37-column summary
    │   └── esv_staged.overview.tsv        # All samples combined
    └── validation/              # (--validate only)
        ├── <sample>_validation_contigs.pdf    # Contig alignment plot (multi-page PDF)
        └── <sample>/
            ├── <sample>_validation_summary.tsv  # All BLAST hits for the sample
            ├── assembly/        # SPAdes contigs, host-filtered contigs (*_query.nohost.fasta),
            │                    # host contig report (*_host_contigs.tsv) and raw host BLAST hits
            ├── blast/           # Per-species BLAST TSVs + ref_lengths + has_contigs
            ├── reads/           # Per-species R1/R2 FASTQs extracted for assembly
            └── genomes/         # (--scaffold) combined, read-polished genomes
                ├── <sample>_<taxon>[__seg<label>].fasta   # one per virus (segment)
                ├── <sample>_genomes.fasta                  # all genomes of the sample
                ├── <sample>_genome_stats.tsv               # statistics per genome
                └── <sample>_scaffold_contigs.tsv           # which contigs were used
        genome_summary.tsv       # (--scaffold) statistics of all genomes in the run
```

---

## 9. Output Files Reference

### `overview/<sample>_overview.tsv`

The primary result file. See [Section 11](#11-overview-table-columns) for a full column description.

### `validation/<sample>_validation_contigs.pdf`

A multi-page PDF with one page per detected viral species. Each contig assembled by SPAdes is drawn as a horizontal bar spanning its aligned region on the best reference genome. Bars are coloured by matched reference accession and ordered by alignment length (longest first). The x-axis spans 0 to the true reference genome length.

### `validation/<sample>/<sample>_validation_summary.tsv`

All BLAST hits across all species for the sample. Columns:

| Column | Description |
|---|---|
| `Sample` | Sample identifier |
| `Species` | Safe-name species string (underscores) |
| `Scaffold_ID` | SPAdes contig name |
| `Matched_Reference` | BLAST subject accession |
| `Identity_%` | Nucleotide identity percentage |
| `Align_Len` | Alignment length (bp) |
| `Query_Len` | Contig length (bp) |
| `Mismatches` | Number of mismatches |
| `Gap_Opens` | Number of gap openings |
| `Q_Start` / `Q_End` | Contig alignment coordinates |
| `S_Start` / `S_End` | Reference alignment coordinates |
| `E-value` | BLAST E-value |
| `Bit_Score` | BLAST bit score |
| `Cov_%` | Query coverage percentage |

### `validation/<sample>/genomes/` and `validation/genome_summary.tsv`

Polished genomes are in `validation/<sample>/genomes/<sample>_<group_id>.fasta` (one per virus, or per segment), with headers `<sample>|<taxon>|<segment>|<reference>`; `<sample>_genomes.fasta` holds all of them. `<sample>_genome_stats.tsv` (all samples combined in `validation/genome_summary.tsv`) has one row per genome:

| Column | Description |
|---|---|
| `sample`, `safe_taxon`, `segment`, `group_id` | Group identity; `segment` is empty for non-segmented viruses |
| `ref_accession`, `ref_length` | Reference the contigs were ordered against |
| `n_contigs` / `n_contigs_used` | Contigs in the group / contigs minimap2 aligned to the reference |
| `contig_bp`, `contig_cov_pct` | Reference bases covered by contigs, and as % of the reference |
| `genome_len`, `n_count`, `N_pct` | Final polished genome length, `N` count and `N` % |
| `ref_cov_pct` | Called (non-`N`) bases as % of the reference length, capped at 100 |
| `mean_depth` | Mean read depth over the genome's called positions (depth ≥ `--polish_min_depth`) in the final polish round |

---

## 10. Validation Sub-workflow

Enable with `--validate`. The current sub-workflow runs once per sample on the trimmed, deduplicated, host-filtered read pair produced by `HOST_FILTER`.

### Step-by-step

1. **De novo assembly** (`SPADES_ASSEMBLY`): SPAdes is run on the full sample-level deduplicated read pair from `FASTP_DEDUP`, using the configured `spades_mode` (`meta` by default, or `rnaviral` for HEV databases when auto-detected). Assembly is skipped (empty query FASTA produced) if the sample has fewer than `--validate_min_reads` read pairs. When the sample exceeds `--validate_spades_max_pairs`, the FASTQ pair is deterministically subsampled with `seqkit sample` before assembly. `SPADES_ASSEMBLY` has an explicit `220.GB` memory allocation so the SPAdes `-m` setting resolves to `220`. Contigs shorter than `--validate_min_contig_len` bp are filtered out with `seqkit seq --min-len`. If SPAdes produces no contigs meeting the length threshold the output FASTA is empty and the downstream BLAST step is skipped gracefully.

   > **Note on assembly failures**: metaSPAdes requires successful insert-size estimation, which depends on FR-oriented read pairs with insert sizes larger than the reads themselves. For virus groups where the template fragments are very short (insert size ≈ read length), R1/R2 pairs may overlap completely, preventing insert-size estimation and resulting in 0 assembled contigs. This is a library preparation characteristic, not a pipeline bug. The `assembly_status` column will report `no_contigs_assembled` in this case.

   **Host contig filter** (`FILTER_HOST_CONTIGS`, `--host_contig_filter`, on by default): removes human rRNA that assembles into contigs. Reads from the multi-copy rRNA arrays can pass `HOST_FILTER` because their mates align to different repeat copies, so the pair is never concordant. Their contigs BLAST-hit viral references whose GenBank records carry human rRNA at their ends (e.g. OL738674.1, KJ716849.1), giving recurring `artefact_suspected` rows. The contigs are searched with megablast against the host BLAST database (`--host_blastdb`); a contig is removed when host hits with ≥ `--host_contig_min_identity` % identity cover ≥ `--host_contig_max_pct` % of its length. The kept contigs (`<sample>_query.nohost.fasta`) go to BLAST and scaffolding; `<sample>_host_contigs.tsv` lists every contig with its host coverage, best host hit and whether it was removed.

2. **BLAST validation** (`BLASTN_VALIDATE`): the assembled contigs are BLASTed against the full EsViritu `.fna` file with `blastn -task blastn` (11 nt seeds; the default megablast needs 28 nt exact matches and misses divergent viruses). Each contig is assigned exclusively to the database taxon with the highest total BLAST bitscore. Results use E-value ≤ 1×10⁻⁵.

   **Protein-level rescue** (`--tblastx_rescue`, on by default, experimental): contigs with no blastn hit, or whose hits cover less than `--verdict_artefact_max_aln_pct` % / `--verdict_artefact_max_aln_bp` bp of the contig, are searched with `tblastx` (both sides translated in six frames) against the same database, to catch viruses too divergent for a nucleotide search. A candidate passes when its hits to the best accession (highest total bitscore) cover at least the same thresholds. All candidates with a protein hit are listed in `blast/<sample>_tblastx.tsv` (raw HSPs in `<sample>_tblastx_raw.tsv`); passing hits appear in the overview (see section 11).

3. **Visualise** (`VISUALIZE_VALIDATION`): one PDF page per matched taxon showing contig coverage of the reference genome.

4. **Combine contigs per taxon** (`--scaffold`, on by default). The contigs of one virus in a sample are joined into one genome against a reference, then corrected with the sample's reads:
   - `SCAFFOLD_GROUPS` first leaves out artefact contigs: a contig is only scaffolded when its BLAST hits cover at least `--verdict_artefact_max_aln_pct` % (default 20) of its length and at least `--verdict_artefact_max_aln_bp` (default 200) bp — the same thresholds the overview uses for `artefact_suspected`. This stops long non-viral contigs (e.g. bacterial rRNA) whose only match is a short stretch at the end of a viral reference from being placed on that virus. The decision per contig is in `validation/<sample>/genomes/<sample>_scaffold_contigs.tsv`. It then groups the remaining contigs per taxon, or per segment for segmented viruses (segment labels harmonised as in the overview table; contigs on an accession without a segment label form their own `?` group). Each group's reference is the accession with the highest total BLAST bitscore in the group, so every segment of a reassortant keeps its own reference.
   - `SCAFFOLD_MINIMAP2` aligns the contigs to the reference (`minimap2 -x map-ont -k 11 -w 5`, which aligns contigs only ~80–86% identical to their reference end to end; the assembly preset `asm20` dropped much of such contigs) and builds a contig consensus in reference coordinates with `samtools consensus`; overlapping contigs merge. Gaps and genome ends are filled with the reference (lowercase) to make a polishing template.
   - `POLISH_SCAFFOLD` maps the sample's host-filtered reads to all templates of the sample at once (`bowtie2 --local`) and calls a new consensus (`samtools consensus`), for `--polish_rounds` rounds; rounds after the first use only the read pairs that mapped in round 1. Positions with fewer than `--polish_min_depth` reads become `N`, so every base in the final genome is supported by reads, never copied from the reference. End `N`s are trimmed.

   Because the template covers the whole reference (gaps and both ends filled with reference sequence), reads can extend the genome beyond the assembled contigs wherever the sample's virus is similar enough to the reference for them to align; elsewhere the genome stays `N`. RagTag and ABACAS were tested as alternatives (October 2026): with their templates ending at the first/last contig the genomes were shorter, and RagTag placed long non-viral contigs on a virus from a short match.

### `assembly_status` values

| Value | Meaning |
|---|---|
| `too_few_reads` | Fewer than `validate_min_reads` deduplicated read pairs available — assembly not attempted |
| `no_contigs_assembled` | SPAdes ran but produced no contigs ≥ `validate_min_contig_len` bp |
| `no_blast_hits` | Contigs were assembled but none had significant BLAST hits against the species reference database |
| `assembled` | At least one contig had a significant BLAST hit |

### Taxonomic level for assembly (`--assembly_taxon_level`)

By default (`subspecies`), reads are grouped and assembled at the finest available taxonomic resolution — using the subspecies rank when EsViritu assigns one (e.g. `hepatitis C virus genotype 1a`), falling back to species otherwise. The subspecies label is prefixed with the species name (e.g. `Rotavirus alphagastroenteritidis 1`), because some subspecies labels (e.g. serotype `1`) are shared by unrelated species; labels that already contain the species name are used as they are.

Use `--assembly_taxon_level species` to always group at species level.

---

## 11. Overview Table Columns

The overview table (`overview/<sample>_overview.tsv`) contains 37 columns, in this order:

1. **Findings:** `sample_ID`, `virus_name`, `verdict`, `esv_verdict`, `blast_verdict`, `flags`
2. **Mapping evidence (EsViritu):** `esv_read_count`, `esv_breadth_pct`, `RPM`, `esv_flags`
3. **De novo evidence (SPAdes + BLAST):** `blast_coverage`, `n_contigs`, `blast_flags`
4. **Technical detail, EsViritu:** `esv_accession`, `genome_length_bp`, `esv_covered_bases`, `esv_ani`, `pi`, `RPKMF`, `RPKMR`
5. **Technical detail, SPAdes + BLAST:** `assembly_status`, `best_blast_reference`, `longest_contig_bp`, `contig_aln_bp`, `contig_aln_pct`, `blast_identity_pct`
6. **Taxonomy:** `family`, `genus`, `species`, `subspecies`
7. **Read funnel** (sample-level, same on every row of a sample)

The sections below describe the columns grouped by source.

### Identity

| Column | Description |
|---|---|
| `sample_ID` | Sample identifier |
| `virus_name` | Display name: species + subspecies at `subspecies` level (e.g. `Alphainfluenzavirus influenzae H3N2`), otherwise species |

### Verdict

Mapping evidence (EsViritu) and de novo evidence (SPAdes + BLAST) are judged separately, then combined. Thresholds are the `--verdict_*` parameters.

| Column | Description |
|---|---|
| `verdict` | Combined call: `confirmed` (`blast_verdict` confirmed) > `probable` (`esv_verdict` supported) > `artefact_suspected` (BLAST-only hit judged an artefact) > `weak` (everything else) |
| `esv_verdict` | Mapping evidence: `supported` if `esv_read_count` ≥ `verdict_esv_min_reads` and `esv_breadth_pct` ≥ `verdict_esv_min_breadth_pct`, otherwise `weak`. `NA` for BLAST-only rows |
| `blast_verdict` | De novo evidence, judged on the best contig (`contig_aln_bp`, `contig_aln_pct`): `confirmed` (≥ `verdict_contig_min_aln_bp` bp and ≥ `verdict_contig_min_aln_pct` % of the contig aligned), `artefact_suspected` (< `verdict_artefact_max_aln_bp` bp or < `verdict_artefact_max_aln_pct` % aligned: a short viral match inside a longer, likely non-viral contig), `inconclusive` (in between), `no_contig` (no contig assigned to this taxon). `NA` without `--validate` |
| `flags` | General flags: `blast_only` (no EsViritu hit) and `in_<k>/<N>_samples` when the taxon is reported in at least `--verdict_recurrent_min_samples` samples of the batch (informational: a real outbreak and a shared contaminant both recur). Empty when no flag applies |
| `esv_flags` | Reasons behind `esv_verdict`, each with the value and the threshold it failed: `low_reads(5<10)`, `low_breadth(2.35%<5%)`, `divergent(86.57%<90%)` (EsViritu read identity). Empty when none apply |
| `blast_flags` | Reasons behind `blast_verdict`, each with the value and the threshold it failed: `short_alignment(81bp<200bp)`, `low_contig_cov(1.9%<20%)` (best contig), `divergent(82.18%<90%)` (BLAST identity). Empty when none apply |

### Read funnel

| Column | Description |
|---|---|
| `raw_reads` | Total reads (R1 + R2) before any filtering |
| `trimmed_reads` | Reads after quality trimming (first step) |
| `trim_removed_pct` | Percentage of raw reads removed by trimming |
| `dedup_reads` | Reads after deduplication |
| `dup_rate_pct` | Percentage of trimmed reads identified as duplicates |
| `host_filtered_reads` | Reads remaining after host removal = reads analysed by EsViritu and SPAdes |
| `host_removal_pct` | Percentage of deduplicated reads removed as host |

### EsViritu detection

| Column | Description |
|---|---|
| `family` | Viral family |
| `genus` | Viral genus |
| `species` | Viral species (ICTV taxonomy, prefix stripped) |
| `subspecies` | Subspecies / strain-level classification when available; `NA` if absent from the EsViritu database entry |
| `esv_accession` | EsViritu reference accession(s) matched for this row; comma-separated for segmented viruses. Multiple rows for the same subspecies indicate distinct reference strains all detected in the sample. |
| `genome_length_bp` | Reference assembly length (sum of all segments for multi-segment viruses) |
| `esv_read_count` | Reads assigned to this virus by EsViritu |
| `esv_covered_bases` | Number of reference bases covered by at least one read |
| `esv_breadth_pct` | Breadth of coverage: `esv_covered_bases / genome_length_bp × 100` |
| `esv_ani` | Average nucleotide identity of aligned reads to the reference |
| `pi` | π (nucleotide diversity): mean pairwise nucleotide differences per site |
| `RPKMF` | Reads Per Kilobase per Million filtered reads (denominator = `host_filtered_reads`) |
| `RPM` | Reads Per Million filtered reads |
| `RPKMR` | Reads Per Kilobase per Million raw reads (denominator = `raw_reads`) |

### Assembly & BLAST validation (`--validate` only)

| Column | Description |
|---|---|
| `assembly_status` | See [Section 10](#assembly_status-values) |
| `n_contigs` | Number of assembled contigs with BLAST hits |
| `longest_contig_bp` | Length of the longest assembled contig |
| `contig_aln_bp` | Query bases of the best contig (most aligned bases) covered by BLAST hits to its assigned reference |
| `contig_aln_pct` | `contig_aln_bp` as % of that contig's length |
| `best_blast_reference` | Non-segmented: accession of the reference with the highest total BLAST bit score. Segmented: best accession per segment, `segment:accession;...` (e.g. `L:KF974361.1;M:KF974359.1`). Segments may come from different assemblies (e.g. reassortants). |
| `blast_coverage` | Reference coverage by assembled contigs. Non-segmented viruses: % of the best reference genome covered (e.g. `92.5%`). Segmented viruses: per-segment coverage, `seg<label>:<pct>%;...` (e.g. `segL:85%;segM:85%;segS:no_hit`); expected segments are those of the best-supported reference assembly, segment labels are harmonised across assemblies (`L RNA` → `L`, `RNA 2` → `2`), and `seg?` = hit to a reference without a segment annotation. `NA` when no contig was assigned |
| `blast_identity_pct` | Nucleotide identity of the best BLAST hit |

### Protein-level rescue (`--validate` with `--tblastx_rescue`, experimental)

| Column | Description |
|---|---|
| `tblastx_n_contigs` | Contigs without nucleotide support that passed the tblastx coverage thresholds for this virus |
| `tblastx_aa_identity_pct` | Amino-acid identity of their tblastx hits (weighted by aligned length) |
| `tblastx_contig_aln_pct` | % of the best such contig covered by tblastx hits |
| `tblastx_reference` | Accession with the highest total tblastx bitscore |

A virus found **only** by tblastx gets its own row with `verdict = protein_hit_only` and `assembly_status = tblastx_only_detection` (a possible divergent/novel virus to review by hand); a virus already in the table gets the columns above and the flag `tblastx_hit`. `protein_hit_only` rows are not shown in the HTML report yet.

### Combined genome (`--validate` with `--scaffold`)

`NA` when no genome was attempted for the virus (no contig passed the scaffolding filter).

| Column | Description |
|---|---|
| `genome_assessment` | Placed right after the verdicts. Overall called (non-`N`) bases over all segments as % of the reference length: `complete` (≥ `--genome_complete_pct`, default 90), `partial` (≥ `--genome_partial_pct`, default 50), `fragmented` (below), `no_genome` (scaffolding ran but no base was called). A segmented virus with an expected segment without any contig (`no_hit` in `blast_coverage`) is at most `partial` |
| `genome_coverage` | Called (non-`N`) bases of the polished genome as % of the reference. Segmented viruses: per segment, `seg<label>:<pct>%;...`, `no_genome` when a segment had contigs but no genome was built |
| `genome_N_pct` | % `N` in the polished genome (all segments together) |
| `genome_mean_depth` | Length-weighted mean read depth |
| `genome_len` | Total polished genome length (bp, all segments) |
| `genome_n_contigs` | Contigs placed on the reference (all segments) |

---

## 12. Automated Production Wrapper

`NGS_wrapper.sh` is a fully automated wrapper for the FHI production environment. It handles the complete upstream/downstream chain:

```
N-drive (SMB)      →  local temp storage  →  Nextflow pipeline  →  N-drive (results)
```

### Usage

```bash
bash NGS_wrapper.sh -r <analysis_name> -a <AGENS> -y <YEAR>
```

**Example:**

```bash
bash NGS_wrapper.sh -r NGS_SEQ-20260210-01 -a UkjentVirus -y 2026
```

**Arguments:**

| Flag | Description |
|---|---|
| `-r` | Analysis run name — used to name the output folder on the N-drive and local status files |
| `-a` | Agens subfolder on the N-drive results tree (e.g. `UkjentVirus`) |
| `-y` | Year subfolder on the N-drive results tree (e.g. `2026`) |

### What the wrapper does

1. **Downloads** the samplesheet (`<RUN>_samplesheet.csv`) from the N-drive.
2. **Downloads** per-sample FASTQ directories: for each row in the samplesheet, the `fastq_dir` path is resolved against the N-drive SMB share and the directory is downloaded into local temp storage.
3. **Builds** the Nextflow-compatible samplesheet (rewrites `fastq_dir` to local paths).
4. **Pulls** the latest pipeline version from GitHub (`alexanderhes/Ukjent_virus`).
5. **Builds** the custom Docker image (`esviritu_pipeline:latest`).
6. **Runs** the Nextflow pipeline.
7. **Uploads** results back to the N-drive under `…/2-Resultater/<AGENS>/<YEAR>/<RUN>/`.
8. **Cleans up** local temp files and Nextflow work directories.

### Status file

The wrapper maintains a per-run status file at `~/esv_<RUN>_status.txt`. This is updated at each major step and on any error. A main log is appended to `~/esv_wrapper.log`.

### N-drive samplesheet format

The samplesheet file (`<RUN>_samplesheet.csv`) must be placed in:
```
Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/<AGENS>/samplesheets/
```

Required format — semicolon-delimited, two columns:

```
sample;fastq_dir
sample1-UV;/mnt/N/Virologi/NGS/0-Sekvenseringsbiblioteker/Illumina_Run/NGS_SEQ-20260210-01/sample1-UV
sample2-UV;/mnt/N/Virologi/NGS/0-Sekvenseringsbiblioteker/Illumina_Run/NGS_SEQ-20260210-01/sample2-UV
```

- `sample` — unique sample identifier (no spaces)
- `fastq_dir` — absolute path to the per-sample FASTQ directory using the `/mnt/N/` mount prefix

Samples can come from different sequencing runs — each row is downloaded independently. The `fastq_dir` directory must contain exactly one file matching `*R1*.fastq.gz` and one matching `*R2*.fastq.gz`. Any BOM (byte order mark) from Windows-created files is stripped automatically.

---

## 13. Resource Configuration

Process resources are assigned by label in `conf/base.config`. Failed processes exit with memory-related codes (104, 134, 137, 139, 143, 247) are automatically retried up to 2 times with increased memory.

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_low` | 2 | 4 GB (×attempt) | 4 h (×attempt) |
| `process_medium` | 6 | 16 GB (×attempt) | 8 h (×attempt) |
| `process_high` | 12 | 32 GB (×attempt) | 24 h (×attempt) |

Label-based values are capped at `--max_cpus`, `--max_memory`, and `--max_time`. `SPADES_ASSEMBLY` additionally has an explicit `220.GB` memory override so the SPAdes `-m` setting matches the host allocation.

### Docker containers

| Process(es) | Container |
|---|---|
| HOST_FILTER, FASTP_TRIM, FASTP_DEDUP, COLLECT_READ_STATS, SUMMARIZE_ESV, VISUALIZE_VALIDATION, SUMMARIZE_VALIDATION | `community.wave.seqera.io/library/bowtie2_esviritu_samtools_seqkit_r-tidyverse:3ee4a52f7d6ae7d9` |
| ESVIRITU | `esviritu_pipeline:latest` (built locally from `docker/Dockerfile`) |
| SPLIT_VIRAL_READS, SPADES_ASSEMBLY, FILTER_HOST_CONTIGS, BLASTN_VALIDATE | `community.wave.seqera.io/library/blast_samtools_seqkit_spades:fc92dccb1ec56163` |
| SUMMARIZE_READ_STATS, MAKE_OVERVIEW_TABLE | `community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186` |
| SCAFFOLD_GROUPS, SCAFFOLD_MINIMAP2, POLISH_SCAFFOLD | default image (`bowtie2_esviritu_samtools_seqkit_r-tidyverse`, includes minimap2 and samtools ≥ 1.16) |

---

## 14. Troubleshooting

### "Expected exactly 1 R1 and 1 R2"

The pipeline skips samples where `fastq_dir` does not contain exactly one file matching `*R1*.fastq.gz` and one matching `*R2*.fastq.gz`. Check that:
- The path in the samplesheet is correct
- Files are named consistently (e.g. `SAMPLEID_R1_001.fastq.gz`)
- No stray FASTQ files exist in the directory

### ESVIRITU fails to find the database

Confirm `--esviritu_db` points to a directory that contains a `.fna` file, a `.mmi` indexed file and a `.tsv` metadata file. The pipeline resolves these by glob (`find -L ... -name "*.fna"`).

### SPAdes assembles 0 contigs for a virus (`no_contigs_assembled`)

metaSPAdes requires paired reads with an insert size larger than the read length (FR orientation, non-overlapping). When viral fragments are very short (insert size ≈ read length), reads overlap entirely and SPAdes cannot estimate insert size — the assembly graph remains unresolved. Possible approaches:
- Confirm the library was prepared with an appropriate insert size for the expected viral genome
- Review `spades.log` in the relevant `work/` directory for the specific SPAdes warning

### `no_blast_hits` despite assembled contigs

BLAST is run against a per-species reference database extracted from the EsViritu `.fna`. If the assembled contigs are divergent enough that no alignment passes E-value ≤ 1×10⁻⁵, no hits are recorded. This may indicate:
- A highly divergent strain not well represented in the reference database
- Mis-assembly artefacts producing non-viral sequence

### `blast_coverage` is lower than `esv_breadth_pct`

`esv_breadth_pct` is based on all reads mapping to the reference (EsViritu metric). `blast_coverage` is based only on assembled contigs ≥ 800 bp with significant BLAST hits. Read-level coverage will always be more complete than contig-level coverage, especially for low-coverage samples where reads are too sparse to assemble long contigs.

### Cleaning up work files

After a successful run, remove intermediate files:

```bash
nextflow clean -f
```

This deletes the `work/` directory contents. The `results/` directory is unaffected.
