#!/usr/bin/env Rscript
# summarize_read_stats.R
#
# Produces two outputs in the working directory (published by SUMMARIZE_READ_STATS):
#
#   esv_staged.read_stats.tsv
#     Batch-level read-count funnel: one row per sample with
#     raw_reads, host_filtered_reads, trimmed_reads, dedup_reads,
#     and the three derived removal percentages.
#
#   esv_staged.detected_virus.info.enriched.tsv
#     The EsViritu per-accession info table joined with raw_reads from the
#     read_stats table. RPM/RPKMR are NOT computed here — make_overview_table.R
#     calculates them at assembly-set level, correctly summing segment lengths
#     for segmented viruses.

suppressPackageStartupMessages(library(tidyverse))

# ── Read funnel table ──────────────────────────────────────────────────────────
read_stats_files <- list.files(".", pattern = "_read_stats\\.tsv$", full.names = TRUE)

read_stats <- map_dfr(read_stats_files, read_tsv,
                      col_types = cols(.default = "c")) %>%
  mutate(across(c(raw_reads, host_filtered_reads, trimmed_reads, dedup_reads,
                  host_removal_pct, trim_removed_pct, dup_rate_pct),
                as.numeric))

write_tsv(read_stats, "esv_staged.read_stats.tsv")

# ── Enrich EsViritu info table with raw_reads ──────────────────────────────────
info <- read_tsv("esv_staged.detected_virus.info.tsv",
                 col_types = cols(.default = "c")) %>%
  mutate(read_count = as.numeric(read_count),
         Length     = as.numeric(Length))

info_enriched <- info %>%
  left_join(read_stats %>% select(sample_ID, raw_reads), by = "sample_ID")

write_tsv(info_enriched, "esv_staged.detected_virus.info.enriched.tsv")
