#!/usr/bin/env Rscript
# prepare_scaffold_groups.R
#
# Groups the BLAST-assigned contigs of one sample into scaffolding groups: one
# group per taxon for non-segmented viruses, one per segment for segmented ones.
# Each group gets the reference accession its contigs will be ordered against:
# the accession with the highest total BLAST bitscore within the group, so for
# segmented viruses (and reassortants) every segment keeps its own reference.
#
# A taxon counts as segmented when any of its matched accessions carries a
# segment label in ref_lengths.tsv. Contigs of a segmented taxon whose accession
# has no label are grouped per accession with segment "?".
#
# Artefact contigs are left out before grouping: a contig is only scaffolded
# when its BLAST hits cover >= <min_aln_pct> % of its length AND >= <min_aln_bp>
# bp (union of hit intervals on the contig). Typical artefacts are long contigs
# (human/bacterial rRNA) whose only match is a short stretch at the end of a
# viral reference; a scaffolder would otherwise place the whole contig on the
# virus. The thresholds are the overview's artefact thresholds
# (verdict_artefact_max_aln_pct / _bp), so a contig is scaffolded exactly when
# it would not be called artefact_suspected.
#
# Usage:
#   prepare_scaffold_groups.R <sample_id> <blastn.tsv> <ref_lengths.tsv> <out.tsv>
#                             <contigs_out.tsv> <min_aln_pct> <min_aln_bp>
#
# Output columns:
#   sample  safe_taxon  segment  group_id  ref_accession  ref_length  n_contigs  contig_ids
#   segment is empty for non-segmented taxa; contig_ids is comma-separated;
#   group_id is a filesystem-safe name used for per-group output files.
#   Taxa whose contigs are all excluded get no group.
# contigs_out.tsv (one row per BLAST-assigned contig):
#   sample  safe_taxon  contig  contig_len  aln_bp  aln_pct  scaffolded

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: prepare_scaffold_groups.R <sample_id> <blastn.tsv> <ref_lengths.tsv> <out.tsv> <contigs_out.tsv> <min_aln_pct> <min_aln_bp>")
}
sample_id    <- args[1]
blast_file   <- args[2]
ref_file     <- args[3]
out_file     <- args[4]
contigs_file <- args[5]
min_aln_pct  <- as.numeric(args[6])
min_aln_bp   <- as.numeric(args[7])
if (anyNA(c(min_aln_pct, min_aln_bp))) stop("min_aln_pct and min_aln_bp must be numeric")

# Copy of union_covered() in bin/make_overview_table.R -- keep in sync.
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

# Copy of normalize_segment() in bin/make_overview_table.R -- keep in sync.
normalize_segment <- function(x) {
  x <- trimws(gsub('"', "", as.character(x)))
  lx <- tolower(x)
  size <- str_match(
    lx,
    "^(?:segment[ -]?)?([lms]|large|medium|middle|small)(?:[ ;-]+(?:rna|segment|large|medium|middle|small))?$"
  )[, 2]
  unknown_num <- str_match(lx, "^unknown ?([0-9]+)$")[, 2]
  num <- str_match(lx, "\\b(?:rna|seg|segment)[ _-]?([0-9]+)\\b")[, 2]
  case_when(
    is.na(x) | !nzchar(x) | lx == "na"        ~ NA_character_,
    size %in% c("l", "large")                 ~ "L",
    size %in% c("m", "medium", "middle")      ~ "M",
    size %in% c("s", "small")                 ~ "S",
    !is.na(unknown_num)                       ~ paste0("unknown", unknown_num),
    !is.na(num)                               ~ num,
    TRUE                                      ~ x
  )
}

to_safe <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

empty_groups <- tibble(
  sample = character(), safe_taxon = character(), segment = character(),
  group_id = character(), ref_accession = character(), ref_length = numeric(),
  n_contigs = integer(), contig_ids = character()
)

blast <- read_tsv(blast_file, col_types = cols(.default = "c"), show_col_types = FALSE)

empty_contigs <- tibble(
  sample = character(), safe_taxon = character(), contig = character(),
  contig_len = integer(), aln_bp = integer(), aln_pct = numeric(),
  scaffolded = logical()
)

if (nrow(blast) == 0) {
  write_tsv(empty_groups, out_file, na = "")
  write_tsv(empty_contigs, contigs_file, na = "")
  quit(save = "no")
}

# Aligned share of each contig (hits to its assigned reference only, as in
# blastn.tsv); contigs below either threshold are not scaffolded.
contig_support <- blast %>%
  group_by(safe_taxon = Species, contig = Scaffold_ID) %>%
  summarise(
    contig_len = as.integer(first(Query_Len)),
    aln_bp = union_covered(Q_Start, Q_End),
    .groups = "drop"
  ) %>%
  mutate(
    sample = sample_id,
    aln_pct = round(aln_bp / contig_len * 100, 1),
    scaffolded = aln_pct >= min_aln_pct & aln_bp >= min_aln_bp
  ) %>%
  select(all_of(names(empty_contigs)))

write_tsv(contig_support, contigs_file, na = "")

blast <- blast %>%
  semi_join(contig_support %>% filter(scaffolded),
            by = c("Species" = "safe_taxon", "Scaffold_ID" = "contig"))

if (nrow(blast) == 0) {
  write_tsv(empty_groups, out_file, na = "")
  quit(save = "no")
}

ref_lengths <- read_tsv(ref_file, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
  transmute(Matched_Reference = Accession,
            ref_length = as.numeric(Length),
            seg = normalize_segment(Segment)) %>%
  distinct(Matched_Reference, .keep_all = TRUE)

hits <- blast %>%
  transmute(safe_taxon = Species, contig = Scaffold_ID, Matched_Reference,
            Bit_Score = as.numeric(Bit_Score)) %>%
  left_join(ref_lengths, by = "Matched_Reference") %>%
  group_by(safe_taxon) %>%
  mutate(segmented = any(!is.na(seg))) %>%
  ungroup() %>%
  mutate(
    segment = case_when(
      !segmented   ~ "",
      is.na(seg)   ~ "?",
      TRUE         ~ seg
    ),
    # Unlabelled accessions of a segmented taxon each form their own group.
    group_key = if_else(segment == "?", paste0("?", Matched_Reference), segment)
  )

best_ref <- hits %>%
  group_by(safe_taxon, group_key, Matched_Reference, ref_length) %>%
  summarise(total_bs = sum(Bit_Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(safe_taxon, group_key) %>%
  slice_max(total_bs, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(safe_taxon, group_key, ref_accession = Matched_Reference, ref_length)

groups <- hits %>%
  group_by(safe_taxon, group_key, segment) %>%
  summarise(
    n_contigs = n_distinct(contig),
    contig_ids = paste(sort(unique(contig)), collapse = ","),
    .groups = "drop"
  ) %>%
  inner_join(best_ref, by = c("safe_taxon", "group_key")) %>%
  mutate(
    sample = sample_id,
    group_id = case_when(
      segment == ""  ~ safe_taxon,
      segment == "?" ~ to_safe(paste0(safe_taxon, "__seg_unlabelled_", ref_accession)),
      TRUE           ~ to_safe(paste0(safe_taxon, "__seg", segment))
    )
  ) %>%
  arrange(safe_taxon, segment) %>%
  select(all_of(names(empty_groups)))

write_tsv(groups, out_file, na = "")
