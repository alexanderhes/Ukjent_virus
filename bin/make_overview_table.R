#!/usr/bin/env Rscript
# make_overview_table.R
#
# Produces a comprehensive per-sample and batch-level overview TSV combining:
#
#   Part 1 - Read funnel:
#     raw_reads -> host_filtered_reads -> trimmed_reads -> dedup_reads,
#     with percentage removal at each step.
#
#   Part 2 - EsViritu results (one row per assembly-set, one row per accession/accession-set):
#     Breadth of coverage (%), RPM, RPKMR, RPKMF, ANI, Pi, read count.
#     Segmented viruses are handled correctly via assembly_summary.tsv
#     (one row per assembly-set, combined genome length across all segments).
#
#   Part 3 - Assembly + BLAST validation (NA when --validate not used):
#     assembly_status, n_contigs, longest_contig_bp,
#     blast_genome_cov_pct (non-segmented) or blast_segment_coverage
#     (per-segment inline string "seg1:X%;seg2:no_hit;..."),
#     blast_identity_pct.
#
# Params are passed as command-line arguments by the Nextflow process:
#   make_overview_table.R <run_validate> <validate_min_reads> <assembly_taxon_level>
# Input TSVs are discovered by filename pattern in the working directory.

suppressPackageStartupMessages(library(tidyverse))

# ── Parameters ────────────────────────────────────────────────────────────────
args                 <- commandArgs(trailingOnly = TRUE)
run_validate         <- as.logical(args[1])
validate_min_reads   <- as.integer(args[2])
assembly_taxon_level <- args[3]

if (!assembly_taxon_level %in% c("species", "subspecies")) {
  stop("assembly_taxon_level must be 'species' or 'subspecies'")
}

read_stats_file <- list.files(".", pattern = "read_stats\\.tsv$", full.names = TRUE)[1]
asm_sum_file <- list.files(".", pattern = "assembly_summary\\.tsv$", full.names = TRUE)[1]
info_enrich_file <- list.files(".", pattern = "info\\.enriched\\.tsv$", full.names = TRUE)[1]
db_meta_file <- list.files(".", pattern = "all_metadata\\.tsv$", full.names = TRUE)[1]

if (is.na(db_meta_file) || !nzchar(db_meta_file)) {
  stop("Could not find staged EsViritu DB metadata TSV")
}

# ── Helpers ───────────────────────────────────────────────────────────────────
strip_prefix <- function(x) sub("^[a-z]__", "", as.character(x))

to_safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

pick_taxon_label <- function(species, subspecies, taxon_level) {
  use_ssp <- taxon_level == "subspecies" &
    !is.na(subspecies) & nzchar(subspecies) & subspecies != "NA"
  ifelse(use_ssp, subspecies, species)
}

first_non_missing <- function(x) {
  values <- x[!is.na(x) & nzchar(x)]
  if (length(values) == 0L) {
    return(NA_character_)
  }
  values[1]
}

unique_or_na <- function(x) {
  values <- unique(x[!is.na(x) & nzchar(x)])
  if (length(values) == 1L) {
    return(values[1])
  }
  NA_character_
}

first_numeric <- function(x) {
  values <- x[!is.na(x)]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  values[1]
}

weighted_identity <- function(identity_pct, align_len) {
  denom <- sum(align_len, na.rm = TRUE)
  if (is.na(denom) || denom == 0) {
    return(NA_real_)
  }
  round(sum(identity_pct * align_len, na.rm = TRUE) / denom, 2)
}

union_covered <- function(starts, ends) {
  s <- pmin(as.integer(starts), as.integer(ends))
  e <- pmax(as.integer(starts), as.integer(ends))
  ok <- !is.na(s) & !is.na(e)
  s <- s[ok]
  e <- e[ok]
  if (length(s) == 0L) return(0L)
  ord <- order(s)
  s <- s[ord]
  e <- e[ord]
  cur_s <- s[1]
  cur_e <- e[1]
  cov <- 0L
  for (i in seq_along(s)) {
    if (s[i] > cur_e) {
      cov <- cov + (cur_e - cur_s + 1L)
      cur_s <- s[i]
      cur_e <- e[i]
    } else {
      cur_e <- max(cur_e, e[i])
    }
  }
  cov + (cur_e - cur_s + 1L)
}

# ── Shared metadata lookup from the full EsViritu DB ──────────────────────────
db_meta_acc <- read_tsv(db_meta_file, col_types = cols(.default = "c"),
                        show_col_types = FALSE) %>%
  mutate(
    Length = as.numeric(Length),
    Asm_length = as.numeric(Asm_length),
    family = strip_prefix(family),
    genus = strip_prefix(genus),
    species = strip_prefix(species),
    subspecies_clean = ifelse(
      !is.na(subspecies) & nzchar(subspecies) & subspecies != "NA",
      strip_prefix(subspecies),
      NA_character_
    ),
    taxon_label = pick_taxon_label(species, subspecies_clean, assembly_taxon_level),
    virus_name = taxon_label,
    safe_taxon = to_safe(taxon_label),
    segment_clean = gsub('"', "", trimws(Segment)),
    is_segmented = !is.na(segment_clean) & nzchar(segment_clean)
  )

taxon_meta <- db_meta_acc %>%
  group_by(safe_taxon) %>%
  summarise(
    virus_name = first_non_missing(virus_name),
    family = first_non_missing(family),
    genus = first_non_missing(genus),
    species = first_non_missing(species),
    subspecies_meta = unique_or_na(subspecies_clean),
    genome_length_bp = first_numeric(Asm_length),
    is_segmented = any(is_segmented),
    .groups = "drop"
  )

seg_map <- db_meta_acc %>%
  filter(is_segmented) %>%
  transmute(
    safe_taxon,
    acc = trimws(Accession),
    seg = trimws(segment_clean)
  ) %>%
  distinct()

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

pi_tbl <- info_raw %>%
  select(sample_ID, Assembly, Pi) %>%
  group_by(sample_ID, Assembly) %>%
  summarise(pi = round(mean(Pi, na.rm = TRUE), 4), .groups = "drop")

asm <- asm_raw %>%
  mutate(
    family = strip_prefix(family),
    genus = strip_prefix(genus),
    species = strip_prefix(species),
    subspecies_col = ifelse(
      !is.na(subspecies) & nzchar(subspecies) & subspecies != "NA",
      strip_prefix(subspecies),
      NA_character_
    ),
    taxon_label = pick_taxon_label(species, subspecies_col, assembly_taxon_level),
    virus_name = taxon_label,
    safe_taxon = to_safe(virus_name),
    segment_clean = gsub('"', "", trimws(Segment)),
    is_segmented = !is.na(segment_clean) & nzchar(segment_clean),
    esv_breadth_pct = round(covered_bases / Asm_length * 100, 2)
  ) %>%
  left_join(pi_tbl, by = c("sample_ID", "Assembly"))

# ── Part 3: BLAST validation ──────────────────────────────────────────────────
blast_empty <- tibble(
  sample_ID = character(),
  safe_taxon = character(),
  n_contigs = integer(),
  longest_contig_bp = integer(),
  best_blast_reference = character(),
  blast_genome_cov_pct = numeric(),
  blast_segment_coverage = character(),
  blast_identity_pct = numeric()
)

if (run_validate) {
  val_files <- list.files(".", pattern = "_blastn\\.tsv$", full.names = TRUE)

  blast_raw <- if (length(val_files) > 0) {
    map_dfr(val_files, function(f) {
      df <- read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
      names(df)[names(df) == "Identity_%"] <- "identity_pct"
      df
    }) %>%
      rename(sample_ID = Sample, safe_taxon = Species) %>%
      filter(!is.na(Matched_Reference)) %>%
      mutate(
        Bit_Score = as.numeric(Bit_Score),
        Align_Len = as.integer(Align_Len),
        Query_Len = as.integer(Query_Len),
        S_Start = as.integer(S_Start),
        S_End = as.integer(S_End),
        identity_pct = as.numeric(identity_pct)
      )
  } else {
    tibble()
  }

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
    tibble(sample_ID = character(), safe_taxon = character(), has_contigs = logical())
  }

  sample_has_contigs_tbl <- has_contigs_tbl %>%
    filter(safe_taxon == "__sample_assembly__") %>%
    transmute(sample_ID, sample_has_contigs = has_contigs) %>%
    distinct()

  if (nrow(blast_raw) > 0) {
    rl_files <- list.files(".", pattern = "_ref_lengths\\.tsv$", full.names = TRUE)
    extra_lens <- if (length(rl_files) > 0) {
      map_dfr(rl_files, function(f) {
        tryCatch(
          read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
            mutate(Length = as.numeric(Length)) %>%
            filter(!is.na(Length), !is.na(Accession)),
          error = function(e) tibble(Accession = character(), Length = numeric())
        )
      })
    } else {
      tibble(Accession = character(), Length = numeric())
    }

    ref_len_map <- bind_rows(
      db_meta_acc %>% select(Accession, Length),
      extra_lens %>% select(Accession, Length),
      info_raw %>% select(Accession, Length)
    ) %>%
      filter(!is.na(Accession)) %>%
      group_by(Accession) %>%
      summarise(Length = first_numeric(Length), .groups = "drop")

    nonseg_taxa <- taxon_meta %>%
      filter(!is_segmented) %>%
      select(safe_taxon) %>%
      distinct()

    blast_ns <- blast_raw %>%
      semi_join(nonseg_taxa, by = "safe_taxon")

    stats_ns <- if (nrow(blast_ns) > 0) {
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
          n_contigs = n_distinct(Scaffold_ID),
          longest_contig_bp = max(Query_Len, na.rm = TRUE),
          blast_identity_pct = weighted_identity(identity_pct, Align_Len),
          cov_bases = union_covered(S_Start, S_End),
          max_s_end = max(S_End, na.rm = TRUE),
          best_blast_reference = first(Matched_Reference),
          .groups = "drop"
        ) %>%
        left_join(ref_len_map, by = c("best_blast_reference" = "Accession")) %>%
        mutate(
          ref_length = coalesce(Length, as.numeric(max_s_end)),
          blast_genome_cov_pct = round(cov_bases / ref_length * 100, 2),
          blast_segment_coverage = NA_character_
        ) %>%
        select(-cov_bases, -max_s_end, -Length, -ref_length)
    } else {
      blast_empty
    }

    seg_taxa <- seg_map %>%
      select(safe_taxon) %>%
      distinct()

    blast_seg <- blast_raw %>%
      semi_join(seg_taxa, by = "safe_taxon") %>%
      inner_join(seg_map, by = c("safe_taxon", "Matched_Reference" = "acc"))

    stats_seg <- if (nrow(blast_seg) > 0) {
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
        mutate(
          cov_pct = ifelse(!is.na(Length) & Length > 0,
                           round(cov_bases / Length * 100, 1),
                           NA_real_)
        )

      seg_cov_str <- seg_map %>%
        distinct(safe_taxon, seg) %>%
        tidyr::crossing(sample_ID = unique(blast_seg$sample_ID)) %>%
        left_join(cov_per_seg %>% select(sample_ID, safe_taxon, seg, cov_pct),
                  by = c("sample_ID", "safe_taxon", "seg")) %>%
        mutate(
          seg_num = suppressWarnings(as.integer(seg)),
          cov_str = ifelse(is.na(cov_pct), "no_hit", paste0(round(cov_pct, 0), "%")),
          seg_str = paste0("seg", seg, ":", cov_str)
        ) %>%
        arrange(sample_ID, safe_taxon, is.na(seg_num), seg_num, seg) %>%
        group_by(sample_ID, safe_taxon) %>%
        summarise(blast_segment_coverage = paste(seg_str, collapse = ";"),
                  .groups = "drop")

      blast_seg %>%
        group_by(sample_ID, safe_taxon) %>%
        summarise(
          n_contigs = n_distinct(Scaffold_ID),
          longest_contig_bp = max(Query_Len, na.rm = TRUE),
          blast_identity_pct = weighted_identity(identity_pct, Align_Len),
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
      distinct() %>%
      select(sample_ID, safe_taxon, n_contigs, longest_contig_bp,
             best_blast_reference, blast_genome_cov_pct,
             blast_segment_coverage, blast_identity_pct)
  } else {
    ref_len_map <- bind_rows(
      db_meta_acc %>% select(Accession, Length),
      info_raw %>% select(Accession, Length)
    ) %>%
      filter(!is.na(Accession)) %>%
      group_by(Accession) %>%
      summarise(Length = first_numeric(Length), .groups = "drop")

    blast_stats <- blast_empty
  }

  blast_part <- asm %>%
    group_by(sample_ID, safe_taxon) %>%
    summarise(esv_read_count = max(read_count), .groups = "drop") %>%
    left_join(blast_stats, by = c("sample_ID", "safe_taxon")) %>%
    left_join(sample_has_contigs_tbl, by = "sample_ID") %>%
    mutate(
      assembly_status = case_when(
        esv_read_count < validate_min_reads    ~ "too_few_reads",
        !is.na(sample_has_contigs) & !sample_has_contigs ~ "no_contigs_assembled",
        is.na(n_contigs) | n_contigs == 0      ~ "no_blast_hits",
        TRUE                                   ~ "assembled"
      )
    ) %>%
    select(-esv_read_count, -sample_has_contigs)
} else {
  ref_len_map <- bind_rows(
    db_meta_acc %>% select(Accession, Length),
    info_raw %>% select(Accession, Length)
  ) %>%
    filter(!is.na(Accession)) %>%
    group_by(Accession) %>%
    summarise(Length = first_numeric(Length), .groups = "drop")

  blast_stats <- blast_empty

  blast_part <- asm %>%
    distinct(sample_ID, safe_taxon) %>%
    mutate(
      assembly_status = NA_character_,
      n_contigs = NA_integer_,
      longest_contig_bp = NA_integer_,
      best_blast_reference = NA_character_,
      blast_genome_cov_pct = NA_real_,
      blast_segment_coverage = NA_character_,
      blast_identity_pct = NA_real_
    )
}

# ── Combine all parts and write outputs ───────────────────────────────────────
overview <- asm %>%
  left_join(read_stats, by = "sample_ID") %>%
  mutate(
    RPM = round(read_count / (raw_reads / 1e6), 2),
    RPKMR = round((read_count / (Asm_length / 1000)) / (raw_reads / 1e6), 2)
  ) %>%
  left_join(blast_part, by = c("sample_ID", "safe_taxon")) %>%
  select(
    sample_ID,
    raw_reads, host_filtered_reads, host_removal_pct,
    trimmed_reads, trim_removed_pct, dedup_reads, dup_rate_pct,
    virus_name, family, genus, species, subspecies = subspecies_col,
    esv_accession = Accession,
    genome_length_bp = Asm_length,
    esv_read_count = read_count,
    esv_covered_bases = covered_bases,
    esv_breadth_pct,
    esv_ani = avg_read_identity,
    pi, RPKMF, RPM, RPKMR,
    assembly_status, n_contigs, longest_contig_bp,
    best_blast_reference, blast_genome_cov_pct,
    blast_segment_coverage, blast_identity_pct
  )

blast_only_rows <- overview[0, ]

if (run_validate && nrow(blast_stats) > 0) {
  esv_taxa <- asm %>%
    distinct(sample_ID, safe_taxon)

  blast_only_rows <- blast_stats %>%
    anti_join(esv_taxa, by = c("sample_ID", "safe_taxon")) %>%
    left_join(taxon_meta, by = "safe_taxon") %>%
    left_join(read_stats, by = "sample_ID") %>%
    left_join(ref_len_map %>%
                rename(best_blast_reference = Accession,
                       best_ref_length = Length),
              by = "best_blast_reference") %>%
    transmute(
      sample_ID,
      raw_reads, host_filtered_reads, host_removal_pct,
      trimmed_reads, trim_removed_pct, dedup_reads, dup_rate_pct,
      virus_name = coalesce(virus_name, gsub("_", " ", safe_taxon)),
      family, genus, species, subspecies = subspecies_meta,
      esv_accession = NA_character_,
      genome_length_bp = as.numeric(coalesce(genome_length_bp, best_ref_length)),
      esv_read_count = NA_real_,
      esv_covered_bases = NA_real_,
      esv_breadth_pct = NA_real_,
      esv_ani = NA_real_,
      pi = NA_real_,
      RPKMF = NA_real_,
      RPM = NA_real_,
      RPKMR = NA_real_,
      assembly_status = "blast_only_detection",
      n_contigs, longest_contig_bp,
      best_blast_reference, blast_genome_cov_pct,
      blast_segment_coverage, blast_identity_pct
    )
}

overview <- bind_rows(overview, blast_only_rows)

write_tsv(overview, "esv_staged.overview.tsv")

for (sid in unique(overview$sample_ID)) {
  write_tsv(filter(overview, sample_ID == sid), paste0(sid, "_overview.tsv"))
}
