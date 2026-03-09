#!/usr/bin/env Rscript
# visualize_validation.R
#
# Reads the per-sample validation BLAST summary TSV produced by the EsViritu
# pipeline and generates a multi-page PDF: one page per detected viral species.
# Each contig is drawn as a horizontal bar spanning its matched region on the
# reference genome.  Bars are coloured by the matched reference accession and
# ordered longest-alignment-first (top of page = best assembled contig).
#
# Usage: visualize_validation.R <input_validation_summary.tsv> <output.pdf>

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: visualize_validation.R <input.tsv> <output.pdf> [ref_lengths_dir]")
}

input_tsv       <- args[1]
output_pdf      <- args[2]
ref_lengths_dir <- if (length(args) >= 3) args[3] else NULL

# ── Load reference sequence lengths ───────────────────────────────────────────
# Produced by BLASTN_VALIDATE (seqkit fx2tab on the per-species reference FASTA).
# Used to set x-axis limits so the full genome is always shown (0 → genome_len).
ref_lengths_tbl <- if (!is.null(ref_lengths_dir)) {
  rl_files <- list.files(ref_lengths_dir, pattern = "_ref_lengths\\.tsv$",
                         full.names = TRUE)
  if (length(rl_files) > 0) {
    suppressWarnings(
      map_dfr(rl_files, function(f) {
        tryCatch(
          read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
            mutate(Length = as.numeric(Length)) %>%
            filter(!is.na(Length), !is.na(Accession)),
          error = function(e) tibble(Accession = character(), Length = numeric())
        )
      }) %>% distinct()
    )
  } else {
    tibble(Accession = character(), Length = numeric())
  }
} else {
  tibble(Accession = character(), Length = numeric())
}

# ── Load data ──────────────────────────────────────────────────────────────────
blast_data <- read_tsv(input_tsv, show_col_types = FALSE) %>%
  rename(
    Identity_pct = `Identity_%`,
    Evalue       = `E-value`,
    Cov_pct      = `Cov_%`
  )

# Empty input guard
if (nrow(blast_data) == 0) {
  pdf(output_pdf, width = 10, height = 5)
  plot.new()
  text(0.5, 0.5, "No BLAST hits found", cex = 1.5, col = "grey50")
  dev.off()
  message("No BLAST hits — blank PDF written.")
  quit(status = 0)
}

# ── Pre-process ────────────────────────────────────────────────────────────────
# Normalise reference coordinates so Plot_Start < Plot_End (handles minus-strand hits)
blast_data <- blast_data %>%
  mutate(
    Plot_Start      = pmin(S_Start, S_End),
    Plot_End        = pmax(S_Start, S_End),
    # Human-readable species name: replace underscores with spaces
    Species_display = str_replace_all(Species, "_", " ")
  )

# Keep best hit per contig per species (highest Bit_Score)
best_hits <- blast_data %>%
  group_by(Sample, Species, Scaffold_ID) %>%
  slice_max(Bit_Score, n = 1, with_ties = FALSE) %>%
  ungroup()

species_list <- unique(best_hits$Species)  # safe-taxon form used for grouping

message(sprintf(
  "Plotting %d species | %d total contigs -> %s",
  length(species_list), nrow(best_hits), output_pdf
))

# ── Plot: per-species variable-height PDFs ────────────────────────────────────
# Each species is written to its own temporary PDF sized to fit its contigs
# (0.40 inch per contig, minimum 3 inches) then merged into one output file.
tmp_pdfs <- character(0)

for (sp in species_list) {

  sp_data <- best_hits %>%
    filter(Species == sp) %>%
    # Longest alignment at the top
    arrange(desc(Align_Len)) %>%
    mutate(Scaffold_Factor = factor(Scaffold_ID, levels = unique(Scaffold_ID)))

  sp_label <- unique(sp_data$Species_display)[1]  # spaces, for titles / messages

  # Per-species page height: 0.40 inch per contig, minimum 3 inches
  n_sp       <- nrow(sp_data)
  page_height <- max(3, n_sp * 0.40 + 2)
  tmp_pdf    <- tempfile(fileext = ".pdf")
  tmp_pdfs   <- c(tmp_pdfs, tmp_pdf)
  pdf(tmp_pdf, width = 13, height = page_height)

  # Determine x-axis upper limit: use true reference genome length when available,
  # fall back to 5% beyond the rightmost observed S_End coordinate.
  best_ref_sp <- sp_data %>%
    mutate(Bit_Score = suppressWarnings(as.numeric(Bit_Score))) %>%
    group_by(Matched_Reference) %>%
    summarise(total_bs = sum(Bit_Score, na.rm = TRUE), .groups = "drop") %>%
    slice_max(total_bs, n = 1, with_ties = FALSE) %>%
    pull(Matched_Reference)

  genome_len <- if (nrow(ref_lengths_tbl) > 0 &&
                    length(best_ref_sp) > 0 &&
                    best_ref_sp %in% ref_lengths_tbl$Accession) {
    ref_lengths_tbl$Length[ref_lengths_tbl$Accession == best_ref_sp][1]
  } else {
    max(sp_data$Plot_End, na.rm = TRUE) * 1.05
  }

  p <- ggplot(sp_data) +
    geom_rect(aes(
      xmin = Plot_Start,
      xmax = Plot_End,
      ymin = as.numeric(Scaffold_Factor) - 0.4,
      ymax = as.numeric(Scaffold_Factor) + 0.4,
      fill = Matched_Reference
    ), alpha = 0.85, colour = "black", linewidth = 0.08) +
    scale_y_continuous(
      breaks = seq_len(nrow(sp_data)),
      labels = sp_data$Scaffold_ID,
      trans  = "reverse",
      expand = c(0.05, 0.05)
    ) +
    scale_x_continuous(
      limits = c(0, genome_len),
      labels = scales::comma,
      expand = expansion(0)
    ) +
    labs(
      title    = paste(unique(sp_data$Sample), "--", sp_label),
      subtitle = sprintf(
        "%d contigs ordered by alignment length (longest first)",
        nrow(sp_data)
      ),
      x    = "Reference position (bp)",
      y    = "Contig ID",
      fill = "Matched reference"
    ) +
    theme_bw() +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.4, "cm"),
      axis.text.y      = element_text(size = 7),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 3, title.position = "top"))

  print(p)
  dev.off()
  message(sprintf("  Plotted: %s (%d contigs, %.1f in tall)", sp_label, n_sp, page_height))
}

# ── Merge per-species PDFs into a single output ───────────────────────────────
merged_ok <- FALSE
if (requireNamespace("pdftools", quietly = TRUE)) {
  tryCatch({
    pdftools::pdf_combine(tmp_pdfs, output = output_pdf)
    merged_ok <- TRUE
  }, error = function(e) message("pdftools merge failed: ", e$message))
}
if (!merged_ok) {
  gs_bin <- Sys.which("gs")
  if (nchar(gs_bin) > 0) {
    ret <- system2(gs_bin,
                   c("-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite",
                     paste0("-sOutputFile=", shQuote(output_pdf)),
                     shQuote(tmp_pdfs)))
    merged_ok <- (ret == 0)
  }
}
if (!merged_ok) {
  stop("Could not merge per-species PDFs: install R package 'pdftools' or ghostscript (gs)")
}

file.remove(tmp_pdfs[file.exists(tmp_pdfs)])
message("Done.")
