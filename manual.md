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
| `--assembly_taxon_level` | `subspecies` | Taxonomic grouping level for read extraction before assembly. `subspecies` uses the finest available rank (recommended for diverse groups such as Enteroviruses and Rotaviruses); `species` always groups at species level. |

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
 HOST_FILTER          bowtie2 --very-sensitive-local (unpaired mode)
       │               removes human T2T + PhiX reads; seqkit pair re-syncs
       ▼
 FASTP_TRIM           adapter auto-detection, quality/length/complexity filter
       │               poly-G/X trimming
       ▼
 FASTP_DEDUP          deduplication on already-trimmed reads
       │               (~40–50% of trimmed reads are duplicates)
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
 MAKE_OVERVIEW_TABLE  30-column overview TSV             SUMMARIZE_VALIDATION  per-sample BLAST summary
                                                         VISUALIZE_VALIDATION  per-sample contig PDF
```

### Host filtering strategy

bowtie2 is run independently on R1 and R2 (unpaired mode) with `--very-sensitive-local`. This judges each read on its own merits, delivering higher sensitivity than paired-mode filtering. After filtering, `seqkit pair` re-synchronises the two files, discarding any read whose mate was removed.

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
    │   ├── <sample>_overview.tsv          # Per-sample 30-column summary
    │   └── esv_staged.overview.tsv        # All samples combined
    └── validation/              # (--validate only)
        ├── <sample>_validation_contigs.pdf    # Contig alignment plot (multi-page PDF)
        └── <sample>/
            ├── <sample>_validation_summary.tsv  # All BLAST hits for the sample
            ├── assembly/        # SPAdes query FASTA per species
            ├── blast/           # Per-species BLAST TSVs + ref_lengths + has_contigs
            └── reads/           # Per-species R1/R2 FASTQs extracted for assembly
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

---

## 10. Validation Sub-workflow

Enable with `--validate`. The sub-workflow runs for every (sample × detected species) pair where EsViritu reports a detection.

### Step-by-step

1. **Read extraction** (`SPLIT_VIRAL_READS`): reads mapping to any reference accession belonging to the target species are extracted from the EsViritu third-pass BAM using `samtools view` + `samtools fastq`. Secondary and supplementary alignments are excluded. Mate synchronisation is enforced.

2. **De novo assembly** (`SPADES_ASSEMBLY`): SPAdes is run in `--meta` mode. Assembly is skipped (empty query FASTA produced) if the species has fewer than `--validate_min_reads` reads. Contigs shorter than `--validate_min_contig_len` bp are filtered out with `seqkit seq --min-len`. If SPAdes produces no contigs meeting the length threshold the output FASTA is empty and the downstream BLAST step is skipped gracefully.

   > **Note on assembly failures**: metaSPAdes requires successful insert-size estimation, which depends on FR-oriented read pairs with insert sizes larger than the reads themselves. For virus groups where the template fragments are very short (insert size ≈ read length), R1/R2 pairs may overlap completely, preventing insert-size estimation and resulting in 0 assembled contigs. This is a library preparation characteristic, not a pipeline bug. The `assembly_status` column will report `no_contigs_assembled` in this case.

3. **BLAST validation** (`BLASTN_VALIDATE`): the assembled contigs are BLASTed against a small per-species reference database extracted from the EsViritu `.fna` file. Only accessions belonging to the detected species are included. Results use E-value ≤ 1×10⁻⁵.

4. **Summarise** (`SUMMARIZE_VALIDATION`): per-species BLAST TSVs are concatenated into a single per-sample summary.

5. **Visualise** (`VISUALIZE_VALIDATION`): one PDF page per species showing contig coverage of the reference genome.

### `assembly_status` values

| Value | Meaning |
|---|---|
| `too_few_reads` | Fewer than `validate_min_reads` reads mapped — assembly not attempted |
| `no_contigs_assembled` | SPAdes ran but produced no contigs ≥ `validate_min_contig_len` bp |
| `no_blast_hits` | Contigs were assembled but none had significant BLAST hits against the species reference database |
| `assembled` | At least one contig had a significant BLAST hit |

### Taxonomic level for assembly (`--assembly_taxon_level`)

By default (`subspecies`), reads are grouped and assembled at the finest available taxonomic resolution — using the subspecies rank when EsViritu assigns one (e.g. `hepatitis C virus genotype 1a`), falling back to species otherwise. 

Use `--assembly_taxon_level species` to always group at species level.

---

## 11. Overview Table Columns

The overview table (`overview/<sample>_overview.tsv`) contains 30 columns:

### Identity

| Column | Description |
|---|---|
| `sample_ID` | Sample identifier |

### Read funnel

| Column | Description |
|---|---|
| `raw_reads` | Total read pairs before any filtering |
| `host_filtered_reads` | Read pairs remaining after host removal |
| `host_removal_pct` | Percentage of raw reads removed as host |
| `trimmed_reads` | Read pairs after quality trimming |
| `trim_removed_pct` | Percentage of host-filtered reads removed by trimming |
| `dedup_reads` | Read pairs after deduplication |
| `dup_rate_pct` | Percentage of trimmed reads identified as duplicates |

### EsViritu detection

| Column | Description |
|---|---|
| `virus_name` | Most informative display name (subspecies when available, otherwise species) |
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
| `RPKMF` | Reads Per Kilobase per Million filtered reads (denominator = `dedup_reads`) |
| `RPM` | Reads Per Million filtered reads |
| `RPKMR` | Reads Per Kilobase per Million raw reads (denominator = `raw_reads`) |

### Assembly & BLAST validation (`--validate` only)

| Column | Description |
|---|---|
| `assembly_status` | See [Section 10](#assembly_status-values) |
| `n_contigs` | Number of assembled contigs with BLAST hits |
| `longest_contig_bp` | Length of the longest assembled contig |
| `best_blast_reference` | Accession of the reference with the highest total BLAST bit score |
| `blast_genome_cov_pct` | Fraction of the best reference genome covered by assembled contigs |
| `blast_segment_coverage` | For segmented viruses: comma-separated list of `segment:contig_count` (e.g. `VP1:2,VP2:1`); `NA` for non-segmented viruses |
| `blast_identity_pct` | Nucleotide identity of the best BLAST hit |

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

All values are capped at `--max_cpus`, `--max_memory`, and `--max_time`.

### Docker containers

| Process(es) | Container |
|---|---|
| HOST_FILTER, FASTP_TRIM, FASTP_DEDUP, COLLECT_READ_STATS, SUMMARIZE_ESV, VISUALIZE_VALIDATION, SUMMARIZE_VALIDATION | `community.wave.seqera.io/library/bowtie2_esviritu_samtools_seqkit_r-tidyverse:3ee4a52f7d6ae7d9` |
| ESVIRITU | `esviritu_pipeline:latest` (built locally from `docker/Dockerfile`) |
| SPLIT_VIRAL_READS, SPADES_ASSEMBLY, BLASTN_VALIDATE | `community.wave.seqera.io/library/blast_samtools_seqkit_spades:fc92dccb1ec56163` |
| SUMMARIZE_READ_STATS, MAKE_OVERVIEW_TABLE | `community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186` |

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

### `blast_genome_cov_pct` is lower than `esv_breadth_pct`

`esv_breadth_pct` is based on all reads mapping to the reference (EsViritu metric). `blast_genome_cov_pct` is based only on assembled contigs ≥ 800 bp with significant BLAST hits. Read-level coverage will always be more complete than contig-level coverage, especially for low-coverage samples where reads are too sparse to assemble long contigs.

### Cleaning up work files

After a successful run, remove intermediate files:

```bash
nextflow clean -f
```

This deletes the `work/` directory contents. The `results/` directory is unaffected.
