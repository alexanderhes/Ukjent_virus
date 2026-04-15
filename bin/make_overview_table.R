#!/usr/bin/env Rscript
# make_overview_table.R
#
# Produces a comprehensive per-sample and batch-level overview TSV combining:
#
#   Part 1 — Read funnel:
#     raw_reads → host_filtered_reads → trimmed_reads → dedup_reads,
#     with percentage removal at each step.
#
#   Part 2 — EsViritu results (one row per assembly-set, one row per accession/accession-set):
#     Breadth of coverage (%), RPM, RPKMR, RPKMF, ANI, Pi, read count.
#     Segmented viruses are handled correctly via assembly_summary.tsv
#     (one row per assembly-set, combined genome length across all segments).
#
#   Part 3 — Assembly + BLAST validation (NA when --validate not used):
#     assembly_status, n_contigs, longest_contig_bp,
#     blast_genome_cov_pct (non-segmented) or blast_segment_coverage
#     (per-segment inline string "seg1:X%;seg2:no_hit;..."),
#     blast_identity_pct.
#
# Params are passed as command-line arguments by the Nextflow process:
#   make_overview_table.R <run_validate> <validate_min_reads>
# Input TSVs are discovered by filename pattern in the working directory.

suppressPackageStartupMessages(library(tidyverse))

# ── Parameters ────────────────────────────────────────────────────────────────
args               <- commandArgs(trailingOnly = TRUE)
run_validate       <- as.logical(args[1])
validate_min_reads <- as.integer(args[2])

read_stats_file  <- list.files(".", pattern = "read_stats\\.tsv$",        full.names = TRUE)[1]
asm_sum_file     <- list.files(".", pattern = "assembly_summary\\.tsv$",   full.names = TRUE)[1]
info_enrich_file <- list.files(".", pattern = "info\\.enriched\\.tsv$",    full.names = TRUE)[1]

# ── Helpers ───────────────────────────────────────────────────────────────────
strip_prefix <- function(x) sub("^[a-z]__", "", as.character(x))

# Replace all non-alphanumeric (except . _ -) characters → underscore.
# IMPORTANT: must match the Nextflow safe-name logic in main.nf exactly:
#   taxon.replaceAll(/^[st]__/, '').replaceAll(/[^A-Za-z0-9._-]/, '_')
# Any divergence between this function and the Groovy code above will
# silently break safe_taxon joins in the overview table.
to_safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Compute total covered bases from a set of (possibly overlapping) intervals.
union_covered <- function(starts, ends) {
  s <- pmin(as.integer(starts), as.integer(ends))
  e <- pmax(as.integer(starts), as.integer(ends))
  ok <- !is.na(s) & !is.na(e)
  s <- s[ok]; e <- e[ok]
  if (length(s) == 0L) return(0L)
  ord <- order(s); s <- s[ord]; e <- e[ord]
  cur_s <- s[1]; cur_e <- e[1]; cov <- 0L
  for (i in seq_along(s)) {
    if (s[i] > cur_e) {
      cov   <- cov + (cur_e - cur_s + 1L)
      cur_s <- s[i]; cur_e <- e[i]
    } else {
      cur_e <- max(cur_e, e[i])
    }
  }
  cov + (cur_e - cur_s + 1L)
}

# ── Part 1: Read funnel ───────────────────────────────────────────────────────
read_stats <- read_tsv(read_stats_file, col_types = cols(.default = "c"),
                       show_col_types = FALSE) %>%
  mutate(across(-sample_ID, as.numeric))

# ── Part 2: EsViritu results ─────────────────────────────────────────────────
asm_raw <- read_tsv(asm_sum_file, col_types = cols(.default = "c"),
                    show_col_types = FALSE) %>%
  mutate(across(c(read_count, covered_bases, Asm_length, RPKMF,
                  avg_read_identity), as.numeric))

info_raw <- read_tsv(info_enrich_file, col_types = cols(.default = "c"),
                     show_col_types = FALSE) %>%
  mutate(across(c(read_count, Length, Pi, raw_reads), as.numeric))

# Per-assembly mean Pi (from per-accession info table)
pi_tbl <- info_raw %>%
  select(sample_ID, Assembly, Pi) %>%
  group_by(sample_ID, Assembly) %>%
  summarise(pi = round(mean(Pi, na.rm = TRUE), 4), .groups = "drop")

asm <- asm_raw %>%
  mutate(
    # Choose most informative display name: subspecies > species.
    use_ssp        = !is.na(subspecies) & nzchar(subspecies),
    virus_name     = ifelse(use_ssp, strip_prefix(subspecies), strip_prefix(species)),
    # Subspecies display name: populated when a distinct subspecies exists;
    # NA when the detection is only resolved to species level.
    subspecies_col = ifelse(use_ssp, strip_prefix(subspecies), NA_character_),
    safe_taxon     = to_safe(virus_name),
    # Detect segmented viruses: Segment column is non-empty
    is_segmented = !is.na(Segment) & nzchar(gsub('"', "", trimws(Segment))),
    # Breadth of coverage = covered_bases / total assembly length (all segments)
    esv_breadth_pct = round(covered_bases / Asm_length * 100, 2),
    family  = strip_prefix(family),
    genus   = strip_prefix(genus),
    species = strip_prefix(species)
  ) %>%
  left_join(pi_tbl, by = c("sample_ID", "Assembly"))

# ── Part 3: BLAST validation ──────────────────────────────────────────────────
# Empty template used when no blast data is available
blast_empty <- tibble(
  sample_ID            = character(),
  safe_taxon           = character(),
  n_contigs            = integer(),
  longest_contig_bp    = integer(),
  best_blast_reference = character(),
  blast_genome_cov_pct = numeric(),
  blast_segment_coverage = character(),
  blast_identity_pct   = numeric()
)

if (run_validate) {
  val_files <- list.files(".", pattern = "_blastn\\.tsv$",
                          full.names = TRUE)

  blast_raw <- if (length(val_files) > 0) {
    map_dfr(val_files, function(f) {
      df <- read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
      # Normalise the special-character column name Identity_% → identity_pct
      names(df)[names(df) == "Identity_%"] <- "identity_pct"
      df
    }) %>%
      rename(sample_ID = Sample, safe_taxon = Species) %>%
      filter(!is.na(Matched_Reference)) %>%
      mutate(
        Bit_Score    = as.numeric(Bit_Score),
        Align_Len    = as.integer(Align_Len),
        Query_Len    = as.integer(Query_Len),
        S_Start      = as.integer(S_Start),
        S_End        = as.integer(S_End),
        identity_pct = as.numeric(identity_pct)
      )
  } else {
    tibble()
  }

  # ── has_contigs table: distinguishes "SPAdes found no contigs" from
  # "contigs assembled but BLAST found no hits" ────────────────────────────
  # Files contain three tab-separated fields: sample_ID, safe_taxon, has_contigs
  hc_files <- list.files(".", pattern = "_has_contigs\\.txt$", full.names = TRUE)
  has_contigs_tbl <- if (length(hc_files) > 0) {
    map_dfr(hc_files, function(f) {
      tryCatch(
        read_tsv(f, col_names = c("sample_ID", "safe_taxon", "has_contigs"),
                 col_types = "ccc", show_col_types = FALSE) %>%
          mutate(has_contigs = tolower(trimws(has_contigs)) == "true"),
        error = function(e) tibble(sample_ID = character(),
                                   safe_taxon = character(),
                                   has_contigs = logical())
      )
    })
  } else {
    tibble(sample_ID = character(), safe_taxon = character(),
           has_contigs = logical())
  }

  if (nrow(blast_raw) > 0) {

    # Accession → Length mapping: start from info table (detected accessions),
    # then supplement with per-species reference lengths from BLASTN_VALIDATE.
    # This ensures non-detected accessions (a closely-related strain chosen as
    # best BLAST reference) receive their true genome length as denominator.
    ref_len_map <- info_raw %>% select(Accession, Length) %>% distinct()
    rl_files <- list.files(".", pattern = "_ref_lengths\\.tsv$", full.names = TRUE)
    if (length(rl_files) > 0) {
      extra_lens <- map_dfr(rl_files, function(f) {
        tryCatch(
          read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
            mutate(Length = as.numeric(Length)) %>%
            filter(!is.na(Length)),
          error = function(e) tibble(Accession = character(), Length = numeric())
        )
      }) %>%
        select(Accession, Length) %>% distinct()
      # Union: prefer known DB length over seqkit measurement when both exist
      ref_len_map <- bind_rows(
        extra_lens,
        ref_len_map
      ) %>%
        group_by(Accession) %>%
        slice(1) %>%
        ungroup()
    }

    # ── Non-segmented viruses ──────────────────────────────────────────────
    nonseg_taxa <- asm %>% filter(!is_segmented) %>%
      select(sample_ID, safe_taxon) %>% distinct()

    blast_ns <- semi_join(blast_raw, nonseg_taxa,
                          by = c("sample_ID", "safe_taxon"))

    stats_ns <- if (nrow(blast_ns) > 0) {
      # Best reference = accession with highest total bit-score across all contigs
      best_ref_ns <- blast_ns %>%
        group_by(sample_ID, safe_taxon, Matched_Reference) %>%
        summarise(total_bs = sum(Bit_Score, na.rm = TRUE), .groups = "drop") %>%
        group_by(sample_ID, safe_taxon) %>%
        slice_max(total_bs, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(sample_ID, safe_taxon, best_ref = Matched_Reference)

      blast_ns %>%
        inner_join(best_ref_ns %>% rename(Matched_Reference = best_ref),
                   by = c("sample_ID", "safe_taxon", "Matched_Reference")) %>%
        group_by(sample_ID, safe_taxon) %>%
        summarise(
          n_contigs            = n_distinct(Scaffold_ID),
          longest_contig_bp    = max(Query_Len,    na.rm = TRUE),
          blast_identity_pct   = round(
            sum(identity_pct * Align_Len, na.rm = TRUE) /
            sum(Align_Len,               na.rm = TRUE), 2),
          cov_bases            = union_covered(S_Start, S_End),
          max_s_end            = max(S_End, na.rm = TRUE),
          best_blast_reference = first(Matched_Reference),
          .groups = "drop"
        ) %>%
        left_join(ref_len_map, by = c("best_blast_reference" = "Accession")) %>%
        mutate(
          # Use known DB length when available; fall back to max(S_End) as proxy
          ref_length             = coalesce(Length, max_s_end),
          blast_genome_cov_pct   = round(cov_bases / ref_length * 100, 2),
          blast_segment_coverage = NA_character_
        ) %>%
        select(-cov_bases, -max_s_end, -Length, -ref_length)
    } else {
      blast_empty
    }

    # ── Segmented viruses ──────────────────────────────────────────────────
    seg_taxa <- asm %>% filter(is_segmented) %>%
      select(sample_ID, safe_taxon) %>% distinct()

    # Expand comma-separated Accession/Segment lists → one row per accession
    seg_map <- asm %>%
      filter(is_segmented) %>%
      select(sample_ID, safe_taxon, accessions = Accession, segments = Segment) %>%
      mutate(
        acc_list = strsplit(accessions, ","),
        seg_list = strsplit(segments,   ",")
      ) %>%
      select(-accessions, -segments) %>%
      unnest(c(acc_list, seg_list)) %>%
      transmute(sample_ID, safe_taxon,
                acc = trimws(acc_list),
                seg = trimws(seg_list))

    blast_seg <- blast_raw %>%
      semi_join(seg_taxa, by = c("sample_ID", "safe_taxon")) %>%
      inner_join(seg_map, by = c("sample_ID", "safe_taxon",
                                 "Matched_Reference" = "acc"))

    stats_seg <- if (nrow(blast_seg) > 0) {
      # Per segment: best accession (max total bit-score) then union coverage
      best_acc_seg <- blast_seg %>%
        group_by(sample_ID, safe_taxon, seg, Matched_Reference) %>%
        summarise(total_bs = sum(Bit_Score, na.rm = TRUE), .groups = "drop") %>%
        group_by(sample_ID, safe_taxon, seg) %>%
        slice_max(total_bs, n = 1, with_ties = FALSE) %>%
        ungroup()

      cov_per_seg <- blast_seg %>%
        inner_join(best_acc_seg %>%
                     select(sample_ID, safe_taxon, seg, Matched_Reference),
                   by = c("sample_ID", "safe_taxon", "seg", "Matched_Reference")) %>%
        left_join(ref_len_map, by = c("Matched_Reference" = "Accession")) %>%
        group_by(sample_ID, safe_taxon, seg, Length) %>%
        summarise(cov_bases = union_covered(S_Start, S_End), .groups = "drop") %>%
        mutate(cov_pct = round(cov_bases / Length * 100, 1))

      # Include all segments even those with no blast hit ("no_hit")
      seg_cov_str <- seg_map %>%
        select(sample_ID, safe_taxon, seg) %>%
        distinct() %>%
        left_join(cov_per_seg %>% select(sample_ID, safe_taxon, seg, cov_pct),
                  by = c("sample_ID", "safe_taxon", "seg")) %>%
        mutate(
          seg_num = suppressWarnings(as.integer(seg)),
          cov_str = ifelse(is.na(cov_pct), "no_hit",
                           paste0(round(cov_pct, 0), "%")),
          seg_str = paste0("seg", seg, ":", cov_str)
        ) %>%
        arrange(sample_ID, safe_taxon, seg_num) %>%
        group_by(sample_ID, safe_taxon) %>%
        summarise(blast_segment_coverage = paste(seg_str, collapse = ";"),
                  .groups = "drop")

      blast_seg %>%
        group_by(sample_ID, safe_taxon) %>%
        summarise(
          n_contigs          = n_distinct(Scaffold_ID),
          longest_contig_bp  = max(Query_Len, na.rm = TRUE),
          blast_identity_pct = round(
            sum(identity_pct * Align_Len, na.rm = TRUE) /
            sum(Align_Len,               na.rm = TRUE), 2),
          .groups = "drop"
        ) %>%
        left_join(seg_cov_str, by = c("sample_ID", "safe_taxon")) %>%
        mutate(
          best_blast_reference = NA_character_,
          blast_genome_cov_pct = NA_real_
        )
    } else {
      blast_empty
    }

    blast_stats <- bind_rows(stats_ns, stats_seg) %>%
      select(sample_ID, safe_taxon, n_contigs, longest_contig_bp,
             best_blast_reference, blast_genome_cov_pct,
             blast_segment_coverage, blast_identity_pct)

  } else {
    blast_stats <- blast_empty
  }

  # Determine assembly_status for every ESV taxon (including too-few-reads cases)
  # Collapse to one row per (sample_ID, safe_taxon) using max read_count across
  # assemblies to avoid many-to-many explosion when multiple assemblies share a taxon.
  blast_part <- asm %>%
    group_by(sample_ID, safe_taxon) %>%
    summarise(esv_read_count = max(read_count), .groups = "drop") %>%
    left_join(blast_stats, by = c("sample_ID", "safe_taxon")) %>%
    left_join(has_contigs_tbl, by = c("sample_ID", "safe_taxon")) %>%
    mutate(
      assembly_status = case_when(
        esv_read_count < validate_min_reads  ~ "too_few_reads",
        !is.na(has_contigs) & !has_contigs   ~ "no_contigs_assembled",
        is.na(n_contigs) | n_contigs == 0    ~ "no_blast_hits",
        TRUE                                 ~ "assembled"
      )
    ) %>%
    select(-esv_read_count, -has_contigs)

} else {
  # --validate not used: all Part 3 columns are NA
  # Use distinct to avoid many-to-many when multiple assemblies share a safe_taxon
  blast_part <- asm %>%
    distinct(sample_ID, safe_taxon) %>%
    mutate(
      assembly_status        = NA_character_,
      n_contigs              = NA_integer_,
      longest_contig_bp      = NA_integer_,
      best_blast_reference   = NA_character_,
      blast_genome_cov_pct   = NA_real_,
      blast_segment_coverage = NA_character_,
      blast_identity_pct     = NA_real_
    )
}

# ── Combine all parts and write outputs ───────────────────────────────────────
overview <- asm %>%
  left_join(read_stats, by = "sample_ID") %>%
  mutate(
    RPM   = round(read_count / (raw_reads / 1e6), 2),
    RPKMR = round((read_count / (Asm_length / 1000)) / (raw_reads / 1e6), 2)
  ) %>%
  left_join(blast_part, by = c("sample_ID", "safe_taxon")) %>%
  select(
    # Identity
    sample_ID,
    # Part 1: Read funnel
    raw_reads, host_filtered_reads, host_removal_pct,
    trimmed_reads, trim_removed_pct, dedup_reads, dup_rate_pct,
    # Part 2: EsViritu
    virus_name, family, genus, species, subspecies = subspecies_col,
    esv_accession    = Accession,
    genome_length_bp = Asm_length,
    esv_read_count   = read_count,
    esv_covered_bases = covered_bases,
    esv_breadth_pct,
    esv_ani = avg_read_identity,
    pi, RPKMF, RPM, RPKMR,
    # Part 3: BLAST / assembly
    assembly_status, n_contigs, longest_contig_bp,
    best_blast_reference, blast_genome_cov_pct,
    blast_segment_coverage, blast_identity_pct
  )

write_tsv(overview, "esv_staged.overview.tsv")

for (sid in unique(overview$sample_ID)) {
  write_tsv(filter(overview, sample_ID == sid), paste0(sid, "_overview.tsv"))
}
