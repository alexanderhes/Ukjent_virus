#!/usr/bin/env Rscript
# visualize_validation.R
#
# Reads the per-sample validation BLAST summary TSV produced by the EsViritu
# pipeline and generates a multi-page PDF: one page per detected viral species.
# Each contig is drawn as a horizontal bar spanning its matched region on the
# reference genome.  Bars are coloured by the matched reference accession and
# ordered longest-alignment-first (top of page = best assembled contig).
# Segmented viruses get one panel per matched segment reference, each spanning
# that segment's full length (segment labels come from the ref_lengths files).
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
# Segment column (raw DB label) is absent in ref_lengths files from older runs.
if (!"Segment" %in% names(ref_lengths_tbl)) ref_lengths_tbl$Segment <- NA_character_
ref_lengths_tbl <- ref_lengths_tbl %>%
  select(Accession, Length, Segment) %>%
  distinct(Accession, .keep_all = TRUE)

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
    left_join(ref_lengths_tbl, by = c("Matched_Reference" = "Accession")) %>%
    # Longest alignment at the top
    arrange(desc(Align_Len)) %>%
    mutate(Scaffold_Factor = factor(Scaffold_ID, levels = unique(Scaffold_ID)))

  sp_label <- unique(sp_data$Species_display)[1]  # spaces, for titles / messages
  n_sp     <- nrow(sp_data)
  tmp_pdf  <- tempfile(fileext = ".pdf")
  tmp_pdfs <- c(tmp_pdfs, tmp_pdf)

  is_segmented <- any(!is.na(sp_data$Segment))

  if (is_segmented) {
    # Segmented virus: one panel per matched reference accession (= segment),
    # each with its own x-axis spanning that segment's full length.
    sp_data <- sp_data %>%
      mutate(
        seg_label = coalesce(Segment, "?"),
        Panel = sprintf("Segment %s (%s)", seg_label, Matched_Reference),
        seg_num = suppressWarnings(as.integer(seg_label))
      ) %>%
      group_by(Matched_Reference) %>%
      mutate(Ref_Len = coalesce(first(Length), max(Plot_End, na.rm = TRUE))) %>%
      ungroup()

    panel_levels <- sp_data %>%
      distinct(Panel, seg_label, seg_num) %>%
      arrange(seg_label == "?", is.na(seg_num), seg_num, seg_label) %>%
      pull(Panel)

    # Discrete y: first level is drawn at the bottom, so reverse to put the
    # longest alignment at the top of each panel.
    sp_data <- sp_data %>%
      mutate(
        Panel = factor(Panel, levels = panel_levels),
        Scaffold_Factor = factor(Scaffold_ID, levels = rev(unique(Scaffold_ID)))
      )

    # Invisible points at 0 and the segment length fix each panel's x range.
    panel_extent <- bind_rows(
      distinct(sp_data, Panel, x_end = Ref_Len),
      distinct(sp_data, Panel) %>% mutate(x_end = 0)
    )

    n_panels    <- length(panel_levels)
    max_per_pan <- max(table(sp_data$Panel))
    page_height <- max(3, n_panels * (max_per_pan * 0.30 + 0.9) + 2)
    pdf(tmp_pdf, width = 13, height = page_height)

    p <- ggplot(sp_data) +
      geom_tile(aes(
        x      = (Plot_Start + Plot_End) / 2,
        width  = Plot_End - Plot_Start,
        y      = Scaffold_Factor,
        height = 0.8,
        fill   = Matched_Reference
      ), alpha = 0.85, colour = "black", linewidth = 0.08) +
      geom_blank(data = panel_extent, aes(x = x_end)) +
      facet_wrap(~Panel, ncol = 1, scales = "free") +
      scale_x_continuous(labels = scales::comma, expand = expansion(0)) +
      labs(
        title    = paste(unique(sp_data$Sample), "--", sp_label),
        subtitle = sprintf(
          "%d contigs on %d segment references; each panel spans the full segment",
          n_sp, n_panels
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
        strip.text       = element_text(size = 8),
        panel.grid.minor = element_blank()
      ) +
      guides(fill = guide_legend(ncol = 3, title.position = "top"))

    print(p)
    dev.off()
    message(sprintf("  Plotted: %s (%d contigs, %d segments, %.1f in tall)",
                    sp_label, n_sp, n_panels, page_height))
    next
  }

  # Per-species page height: 0.40 inch per contig, minimum 3 inches
  page_height <- max(3, n_sp * 0.40 + 2)
  pdf(tmp_pdf, width = 13, height = page_height)

  # x-axis upper limit: longest reference among all matched accessions (contigs
  # may hit different strains of different lengths), and never shorter than
  # the rightmost aligned position so no contig is clipped.
  genome_len <- max(c(sp_data$Length, sp_data$Plot_End), na.rm = TRUE)
  if (all(is.na(sp_data$Length))) genome_len <- genome_len * 1.05

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
