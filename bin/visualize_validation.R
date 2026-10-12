#!/usr/bin/env Rscript
# visualize_validation.R
#
# Draws what was assembled for each virus on its reference genome, as
# multi-page PDFs with one page per detected viral species. Segmented viruses
# get one panel per matched segment reference, each spanning that segment's
# full length (segment labels come from the ref_lengths files).
#
# contigs mode (VISUALIZE_VALIDATION): each contig is a horizontal bar spanning
#   its best BLAST hit on the reference, coloured by the matched reference
#   accession and ordered longest-alignment-first. Contigs are split by the
#   artefact rule of SCAFFOLD_GROUPS (prepare_scaffold_groups.R): a contig whose
#   BLAST hits cover < <min_aln_pct> % of its length or < <min_aln_bp> bp is a
#   likely artefact. Passing contigs go to <contigs.pdf>, likely artefacts to
#   <artefacts.pdf>, so the contig plot shows the contigs the genomes are built
#   from.
#
# genomes mode (VISUALIZE_GENOMES), two PDFs:
#   <genomes.pdf>: one page per species; each full-length polished genome
#     (POLISH_SCAFFOLD) is drawn on the reference it was built against, every
#     base at its reference position (positions file). Bases are coloured by
#     source: contig base kept (too few reads), read-called inside a contig,
#     read-called outside the contigs (what the reads added) or N (case and
#     in_contig in the positions file). End Ns are part of the full-length genome; contig
#     overhangs extend beyond the reference ends (dotted lines).
#   <coverage.pdf>: one page per genome (per segment for segmented viruses),
#     three tracks on the reference axis: read depth on the genome (dashed line
#     at <min_depth>), the consensus by source (N in black) and the contigs the
#     genome was built from (their minimap2 placement in SCAFFOLD_MINIMAP2).
#
# Usage:
#   visualize_validation.R contigs <blastn.tsv> <ref_lengths_dir> <contigs.pdf>
#                                  <artefacts.pdf> <min_aln_pct> <min_aln_bp>
#   visualize_validation.R genomes <genome_stats.tsv> <genome_positions.tsv>
#                                  <contig_ref.tsv> <ref_lengths_dir> <genomes.pdf>
#                                  <coverage.pdf> <min_depth>

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

usage <- paste(
  "Usage:",
  "  visualize_validation.R contigs <blastn.tsv> <ref_lengths_dir> <contigs.pdf> <artefacts.pdf> <min_aln_pct> <min_aln_bp>",
  "  visualize_validation.R genomes <genome_stats.tsv> <genome_positions.tsv> <contig_ref.tsv> <ref_lengths_dir> <genomes.pdf> <coverage.pdf> <min_depth>",
  sep = "\n"
)
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else ""
if (!(mode == "contigs" && length(args) == 7) && !(mode == "genomes" && length(args) == 8)) {
  stop(usage)
}

# ── Reference sequence lengths ────────────────────────────────────────────────
# Produced by BLASTN_VALIDATE (seqkit fx2tab on the per-species reference FASTA).
# Used to set x-axis limits so the full genome is always shown (0 → genome_len).
load_ref_lengths <- function(ref_lengths_dir) {
  rl_files <- list.files(ref_lengths_dir, pattern = "_ref_lengths\\.tsv$",
                         full.names = TRUE)
  tbl <- if (length(rl_files) > 0) {
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
  # Segment column (raw DB label) is absent in ref_lengths files from older runs.
  if (!"Segment" %in% names(tbl)) tbl$Segment <- NA_character_
  tbl %>%
    select(Accession, Length, Segment) %>%
    distinct(Accession, .keep_all = TRUE)
}

blank_pdf <- function(output_pdf, text) {
  pdf(output_pdf, width = 10, height = 5)
  plot.new()
  text(0.5, 0.5, text, cex = 1.5, col = "grey50")
  dev.off()
  message(text, " -- blank PDF written: ", output_pdf)
}

# ── One species page ──────────────────────────────────────────────────────────
# d: one species, one or more bars per row, with columns Sample,
# Species_display, Row_ID (factor, levels in top-to-bottom order), Plot_Start,
# Plot_End, Fill, Matched_Reference, Length (reference length) and Segment (raw
# DB label, NA for non-segmented references). Plot_Start may be negative (contig
# overhangs before the reference start); ref_ends draws dotted lines at 0 and
# the reference length.
draw_species_page <- function(d, tmp_pdf, title, subtitle_flat, subtitle_seg,
                              y_lab, fill_lab, fill_scale = NULL, outline = "black",
                              ref_ends = FALSE) {
  n_rows <- nlevels(droplevels(d$Row_ID))
  theme_page <- theme_bw() +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.4, "cm"),
      axis.text.y      = element_text(size = 7),
      strip.text       = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  if (any(!is.na(d$Segment))) {
    # Segmented virus: one panel per matched reference accession (= segment),
    # each with its own x-axis spanning that segment's full length.
    d <- d %>%
      mutate(
        seg_label = coalesce(Segment, "?"),
        Panel = sprintf("Segment %s (%s)", seg_label, Matched_Reference),
        seg_num = suppressWarnings(as.integer(seg_label))
      ) %>%
      group_by(Matched_Reference) %>%
      mutate(Ref_Len = coalesce(first(Length), max(Plot_End, na.rm = TRUE))) %>%
      ungroup()

    panel_levels <- d %>%
      distinct(Panel, seg_label, seg_num) %>%
      arrange(seg_label == "?", is.na(seg_num), seg_num, seg_label) %>%
      pull(Panel)

    # Discrete y: first level is drawn at the bottom, so reverse to put the
    # first row at the top of each panel.
    d <- d %>%
      mutate(
        Panel = factor(Panel, levels = panel_levels),
        Row_ID = factor(Row_ID, levels = rev(levels(Row_ID)))
      )

    # Invisible points at 0 (or the first overhang base) and the segment length
    # fix each panel's x range.
    panel_extent <- bind_rows(
      distinct(d, Panel, x_end = Ref_Len),
      d %>% group_by(Panel) %>% summarise(x_end = min(0, Plot_Start), .groups = "drop")
    )
    end_lines <- if (ref_ends) {
      geom_vline(data = bind_rows(distinct(d, Panel, x = Ref_Len), distinct(d, Panel) %>% mutate(x = 0)),
                 aes(xintercept = x), linetype = "dotted", colour = "grey30", linewidth = 0.4)
    } else NULL

    n_panels    <- length(panel_levels)
    max_per_pan <- max(tapply(as.character(d$Row_ID), d$Panel, n_distinct), na.rm = TRUE)
    page_height <- max(3, n_panels * (max_per_pan * 0.30 + 0.9) + 2)
    pdf(tmp_pdf, width = 13, height = page_height)

    p <- ggplot(d) +
      geom_tile(aes(
        x      = (Plot_Start + Plot_End) / 2,
        width  = Plot_End - Plot_Start,
        y      = Row_ID,
        height = 0.8,
        fill   = Fill
      ), alpha = 0.85, colour = outline, linewidth = 0.08) +
      geom_blank(data = panel_extent, aes(x = x_end)) +
      end_lines +
      facet_wrap(~Panel, ncol = 1, scales = "free") +
      scale_x_continuous(labels = scales::comma, expand = expansion(0)) +
      labs(title = title, subtitle = subtitle_seg(n_rows, n_panels),
           x = "Reference position (bp)", y = y_lab, fill = fill_lab) +
      fill_scale +
      theme_page +
      guides(fill = guide_legend(ncol = 3, title.position = "top"))

    print(p)
    dev.off()
    return(sprintf("%d rows, %d segments, %.1f in tall", n_rows, n_panels, page_height))
  }

  # Page height: 0.40 inch per row, minimum 3 inches
  page_height <- max(3, n_rows * 0.40 + 2)
  pdf(tmp_pdf, width = 13, height = page_height)

  # x-axis upper limit: longest reference among all matched accessions (contigs
  # may hit different strains of different lengths), and never shorter than
  # the rightmost aligned position so no bar is clipped.
  genome_len <- max(c(d$Length, d$Plot_End), na.rm = TRUE)
  if (all(is.na(d$Length))) genome_len <- genome_len * 1.05
  x_start <- min(0, d$Plot_Start, na.rm = TRUE)
  end_lines <- if (ref_ends) {
    geom_vline(xintercept = c(0, max(d$Length, na.rm = TRUE)), linetype = "dotted",
               colour = "grey30", linewidth = 0.4)
  } else NULL

  p <- ggplot(d) +
    geom_rect(aes(
      xmin = Plot_Start,
      xmax = Plot_End,
      ymin = as.numeric(Row_ID) - 0.4,
      ymax = as.numeric(Row_ID) + 0.4,
      fill = Fill
    ), alpha = 0.85, colour = outline, linewidth = 0.08) +
    end_lines +
    scale_y_continuous(
      breaks = seq_len(nlevels(d$Row_ID)),
      labels = levels(d$Row_ID),
      trans  = "reverse",
      expand = c(0.05, 0.05)
    ) +
    scale_x_continuous(
      limits = c(x_start, genome_len),
      labels = scales::comma,
      expand = expansion(0)
    ) +
    labs(title = title, subtitle = subtitle_flat(n_rows),
         x = "Reference position (bp)", y = y_lab, fill = fill_lab) +
    fill_scale +
    theme_page +
    guides(fill = guide_legend(ncol = 3, title.position = "top"))

  print(p)
  dev.off()
  sprintf("%d rows, %.1f in tall", n_rows, page_height)
}

# ── Merge per-species PDFs into a single output ───────────────────────────────
merge_pdfs <- function(tmp_pdfs, output_pdf) {
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
  invisible(file.remove(tmp_pdfs[file.exists(tmp_pdfs)]))
}

# Each species is written to its own temporary PDF sized to fit its rows, then
# merged into one output file. title_suffix is appended to "<sample> -- <species>".
write_species_pdf <- function(bars, output_pdf, empty_text, title_suffix = "", ...) {
  if (nrow(bars) == 0) {
    blank_pdf(output_pdf, empty_text)
    return(invisible())
  }
  species_list <- unique(bars$Species)
  message(sprintf("Plotting %d species -> %s", length(species_list), output_pdf))
  tmp_pdfs <- character(0)
  for (sp in species_list) {
    d <- bars %>%
      filter(Species == sp) %>%
      mutate(Row_ID = droplevels(Row_ID))
    sp_label <- d$Species_display[1]
    tmp_pdf  <- tempfile(fileext = ".pdf")
    tmp_pdfs <- c(tmp_pdfs, tmp_pdf)
    info <- draw_species_page(
      d, tmp_pdf,
      title = paste0(paste(d$Sample[1], "--", sp_label), title_suffix),
      ...
    )
    message(sprintf("  Plotted: %s (%s)", sp_label, info))
  }
  merge_pdfs(tmp_pdfs, output_pdf)
}

# ── Contigs mode ──────────────────────────────────────────────────────────────
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

run_contigs <- function(input_tsv, ref_lengths_dir, contigs_pdf, artefacts_pdf,
                        min_aln_pct, min_aln_bp) {
  if (anyNA(c(min_aln_pct, min_aln_bp))) stop("min_aln_pct and min_aln_bp must be numeric")
  ref_lengths_tbl <- load_ref_lengths(ref_lengths_dir)

  blast_data <- read_tsv(input_tsv, show_col_types = FALSE) %>%
    rename(
      Identity_pct = `Identity_%`,
      Evalue       = `E-value`,
      Cov_pct      = `Cov_%`
    )

  if (nrow(blast_data) == 0) {
    blank_pdf(contigs_pdf, "No BLAST hits found")
    blank_pdf(artefacts_pdf, "No BLAST hits found")
    return(invisible())
  }

  # Aligned share of each contig, computed as in prepare_scaffold_groups.R
  # (union of all hit intervals on the contig, % rounded to one decimal), so
  # a contig is plotted as passing exactly when it is scaffolded.
  contig_support <- blast_data %>%
    group_by(Species, Scaffold_ID) %>%
    summarise(
      contig_len = as.integer(first(Query_Len)),
      aln_bp = union_covered(Q_Start, Q_End),
      .groups = "drop"
    ) %>%
    mutate(
      aln_pct = round(aln_bp / contig_len * 100, 1),
      passes = aln_pct >= min_aln_pct & aln_bp >= min_aln_bp
    )

  # Normalise reference coordinates so Plot_Start < Plot_End (handles
  # minus-strand hits), keep the best hit per contig per species (highest
  # Bit_Score) and order longest alignment first.
  best_hits <- blast_data %>%
    mutate(
      Plot_Start      = pmin(S_Start, S_End),
      Plot_End        = pmax(S_Start, S_End),
      # Human-readable species name: replace underscores with spaces
      Species_display = str_replace_all(Species, "_", " ")
    ) %>%
    group_by(Sample, Species, Scaffold_ID) %>%
    slice_max(Bit_Score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(contig_support, by = c("Species", "Scaffold_ID")) %>%
    left_join(ref_lengths_tbl, by = c("Matched_Reference" = "Accession")) %>%
    arrange(Species, desc(Align_Len)) %>%
    mutate(Fill = Matched_Reference)

  n_art <- sum(!best_hits$passes)
  message(sprintf("%d contigs: %d pass, %d likely artefacts (BLAST hits < %s %% or < %s bp of the contig)",
                  nrow(best_hits), nrow(best_hits) - n_art, n_art, min_aln_pct, min_aln_bp))

  passing <- best_hits %>%
    filter(passes) %>%
    mutate(Row_ID = factor(Scaffold_ID, levels = unique(Scaffold_ID)))
  write_species_pdf(
    passing, contigs_pdf,
    empty_text    = "No contigs passed the artefact thresholds",
    subtitle_flat = function(n) sprintf("%d contigs ordered by alignment length (longest first)", n),
    subtitle_seg  = function(n, p) sprintf(
      "%d contigs on %d segment references; each panel spans the full segment", n, p),
    y_lab = "Contig ID", fill_lab = "Matched reference"
  )

  # Artefact rows show why the contig failed: the aligned share of the contig.
  artefacts <- best_hits %>%
    filter(!passes) %>%
    mutate(
      Row_ID = sprintf("%s  (%s of %s bp aligned, %s%%)", Scaffold_ID,
                       scales::comma(aln_bp), scales::comma(contig_len), aln_pct),
      Row_ID = factor(Row_ID, levels = unique(Row_ID))
    )
  rule <- sprintf("BLAST hits cover < %s%% or < %s bp of the contig; not used for the genomes",
                  min_aln_pct, min_aln_bp)
  write_species_pdf(
    artefacts, artefacts_pdf,
    empty_text    = "No likely artefact contigs",
    title_suffix  = " -- likely artefact contigs",
    subtitle_flat = function(n) sprintf("%d contigs: %s", n, rule),
    subtitle_seg  = function(n, p) sprintf("%d contigs on %d segment references: %s", n, p, rule),
    y_lab = "Contig ID", fill_lab = "Matched reference"
  )
}


# ── Genomes mode ──────────────────────────────────────────────────────────────
genome_levels <- c("Contig base (too few reads)", "Read-called, in contig",
                   "Read-called, outside contigs", "N")
genome_colours <- c("Contig base (too few reads)" = "#08519c",
                    "Read-called, in contig" = "#6baed6",
                    "Read-called, outside contigs" = "#fd8d3c", "N" = "grey80")

# Category of each genome base (index into genome_levels) from its case
# (uppercase = contig base kept because too few reads cover it, lowercase =
# read-called) and whether it lies on a contig position of the unpolished
# template (in_contig).
base_category <- function(base, in_contig) {
  case_when(
    grepl("^[ACGT]$", base)                 ~ 1L,
    grepl("^[acgt]$", base) & in_contig == 1 ~ 2L,
    grepl("^[acgt]$", base)                 ~ 3L,
    TRUE                                    ~ 4L
  )
}

# Places one genome on its reference from the reference position of each base
# (POLISH_SCAFFOLD positions file). Returns one row per reference position from
# min(1, first overhang position) to max(ref length, last overhang position):
# category (0 = no genome base there) and depth. Inserted bases ("<p>i") have
# no reference position and are left out; a reference position deleted in the
# genome takes the category of the base before it.
genome_on_reference <- function(g, ref_len) {
  g <- g %>% filter(!grepl("i$", refpos)) %>% mutate(ref_pos = as.integer(refpos))
  x <- seq(min(1L, g$ref_pos), max(ref_len, g$ref_pos))
  out <- tibble(ref_pos = x) %>%
    left_join(g %>% select(ref_pos, base, depth, in_contig), by = "ref_pos") %>%
    mutate(category = base_category(base, in_contig))
  out$category[is.na(out$base)] <- NA_integer_
  # Deletions inside the genome span carry the category before them; positions
  # outside it (only possible when the genome is empty at an end) stay 0.
  inside <- out$ref_pos >= min(g$ref_pos) & out$ref_pos <= max(g$ref_pos)
  ok <- !is.na(out$category)
  filled <- out$category[ok][cumsum(ok)]
  out$category <- if_else(inside, filled, 0L)
  out
}

# Run-length intervals of a category vector along reference positions,
# 0-based start / end for geom_rect (position p spans p - 1 to p).
category_runs <- function(ref_pos, category) {
  r <- rle(category)
  ends <- cumsum(r$lengths)
  tibble(Plot_Start = ref_pos[ends - r$lengths + 1L] - 1, Plot_End = ref_pos[ends],
         value = r$values) %>%
    filter(value > 0) %>%
    mutate(Fill = factor(genome_levels[value], levels = genome_levels))
}

genome_fill <- scale_fill_manual(values = genome_colours, drop = TRUE)

# One coverage page: depth, consensus and contig tracks stacked with gtable so
# their x-axes line up whatever the width of the y-axis labels. Dotted lines
# mark the reference ends; bases beyond them are contig overhangs.
draw_coverage_page <- function(g, pos, contigs, min_depth, tmp_pdf) {
  ref_len <- g$Length
  x_lim   <- c(min(0, min(pos$ref_pos) - 1), max(ref_len, max(pos$ref_pos)))
  x_scale <- scale_x_continuous(limits = x_lim, labels = scales::comma,
                                expand = expansion(0))
  ref_ends <- geom_vline(xintercept = c(0, ref_len), linetype = "dotted",
                         colour = "grey30", linewidth = 0.4)
  theme_track <- theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7),
          legend.position = "bottom", legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))
  no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                axis.title.x = element_blank())

  # Depth on the genome at the reference position of each genome base; breaks
  # where no genome base sits (deletions) so gaps are not bridged.
  depth <- pos %>%
    filter(!is.na(depth)) %>%
    mutate(run = cumsum(c(1L, diff(ref_pos) != 1L)))
  overhang <- if (g$overhang_5p + g$overhang_3p > 0)
    sprintf("Contig overhang %s bp 5', %s bp 3'. ", g$overhang_5p, g$overhang_3p) else ""
  p_depth <- ggplot(depth) +
    geom_ribbon(aes(x = ref_pos, ymin = 0, ymax = depth, group = run),
                fill = "#2c7fb8", alpha = 0.6) +
    geom_hline(yintercept = min_depth, linetype = "dashed", colour = "#d7301f",
               linewidth = 0.4) +
    ref_ends + x_scale +
    scale_y_continuous(labels = scales::comma, expand = expansion(c(0, 0.05))) +
    labs(
      title = sprintf("%s -- %s%s -- %s", g$sample, g$Species_display,
                      if (!is.na(g$Segment)) paste0(" -- segment ", g$Segment) else "",
                      g$ref_accession),
      subtitle = sprintf(
        paste0("Genome %s bp (%s bp trimmed): %s%% of the reference called; N %s inside, %s at the ends; ",
               "mean depth %sx; %s of %s contigs used.\n",
               "%sDashed line: min depth %s. Dotted lines: reference ends."),
        scales::comma(g$genome_len), scales::comma(g$trimmed_len), g$ref_cov_pct,
        scales::comma(g$n_internal), scales::comma(g$n_end_5p + g$n_end_3p),
        round(g$mean_depth), g$n_contigs_used, g$n_contigs, overhang, min_depth),
      y = "Read depth"
    ) +
    theme_track + no_x

  # Consensus by source; N in black so isolated Ns show as lines.
  cons <- category_runs(pos$ref_pos, pos$category)
  p_cons <- ggplot(cons) +
    geom_rect(aes(xmin = Plot_Start, xmax = Plot_End, ymin = 0, ymax = 1, fill = Fill)) +
    ref_ends + x_scale +
    scale_y_continuous(breaks = 0.5, labels = "Consensus", expand = expansion(0)) +
    scale_fill_manual(values = replace(genome_colours, "N", "black"), drop = TRUE) +
    labs(y = NULL, fill = NULL) +
    theme_track + no_x

  # Contigs (aligned part), one row each, in reference order.
  contigs <- contigs %>%
    arrange(ref_start, contig) %>%
    mutate(row = factor(contig, levels = rev(unique(contig))))
  p_contig <- ggplot(contigs) +
    geom_rect(aes(xmin = ref_start - 1, xmax = ref_end,
                  ymin = as.numeric(row) - 0.4, ymax = as.numeric(row) + 0.4),
              fill = "#41ab5d", colour = "black", linewidth = 0.08) +
    ref_ends + x_scale +
    scale_y_continuous(breaks = seq_len(nlevels(contigs$row)), labels = levels(contigs$row),
                       limits = c(0.5, max(1, nlevels(contigs$row)) + 0.5),
                       expand = expansion(0)) +
    labs(x = "Reference position (bp)", y = "Contigs (aligned part)") +
    theme_track

  n_contigs <- max(1, nlevels(contigs$row))
  grobs <- lapply(list(p_depth, p_cons, p_contig), ggplotGrob)
  combined <- do.call(rbind, c(grobs, size = "max"))
  panels <- combined$layout$t[grepl("^panel", combined$layout$name)]
  combined$heights[panels] <- unit(c(3, 0.5, 0.35 * n_contigs), "null")

  page_height <- 4.7 + 0.35 * n_contigs
  pdf(tmp_pdf, width = 13, height = page_height)
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
}

run_genomes <- function(stats_tsv, positions_tsv, contig_tsv, ref_lengths_dir,
                        genomes_pdf, coverage_pdf, min_depth) {
  ref_lengths_tbl <- load_ref_lengths(ref_lengths_dir)
  stats <- read_tsv(stats_tsv, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    mutate(across(c(ref_length, n_contigs, n_contigs_used, genome_len, n_count,
                    ref_cov_pct, mean_depth, trimmed_len, n_end_5p, n_end_3p,
                    n_internal, overhang_5p, overhang_3p), as.numeric)) %>%
    filter(genome_len > 0) %>%
    left_join(ref_lengths_tbl, by = c("ref_accession" = "Accession")) %>%
    mutate(
      Species         = safe_taxon,
      Species_display = str_replace_all(safe_taxon, "_", " "),
      Length          = coalesce(Length, ref_length),
      # Segment from the DB label as in the contig plot; the normalised label
      # from the genome stats when the reference has none in ref_lengths.
      Segment         = coalesce(Segment, na_if(segment, ""))
    ) %>%
    arrange(Species, Segment)

  # Depth, base (case = source) and reference position per genome base
  positions <- read_tsv(positions_tsv, col_types = "ciicci", show_col_types = FALSE)
  contig_aln <- read_tsv(contig_tsv, col_types = "ccciic", show_col_types = FALSE) %>%
    select(group_id, contig, ref_start, ref_end)

  stats <- stats %>% filter(group_id %in% positions$group_id)
  if (nrow(stats) == 0) {
    blank_pdf(genomes_pdf, "No genomes assembled")
    blank_pdf(coverage_pdf, "No genomes assembled")
    return(invisible())
  }

  pos <- stats %>%
    select(group_id, Length) %>%
    mutate(p = map2(group_id, Length, ~ genome_on_reference(
      positions %>% filter(group_id == .x), .y))) %>%
    select(group_id, p) %>%
    unnest(p)

  # ── Genome overview: one page per species ──
  bars <- pos %>%
    group_by(group_id) %>%
    summarise(runs = list(category_runs(ref_pos, category)), .groups = "drop") %>%
    unnest(runs) %>%
    inner_join(stats, by = "group_id") %>%
    mutate(
      Sample            = sample,
      Matched_Reference = ref_accession,
      Row_ID = sprintf(
        "%s bp (%s trimmed), %s%% of ref called, N %s inside / %s at ends, %sx%s",
        scales::comma(genome_len), scales::comma(trimmed_len), ref_cov_pct,
        scales::comma(n_internal), scales::comma(n_end_5p + n_end_3p),
        round(mean_depth),
        if_else(overhang_5p + overhang_3p > 0,
                sprintf("\ncontig overhang %s bp 5', %s bp 3'", overhang_5p, overhang_3p), "")
      )
    ) %>%
    arrange(Species, Segment)

  # Every row is one genome, so its stats label is unique within a species page.
  bars$Row_ID <- factor(bars$Row_ID, levels = unique(bars$Row_ID))

  write_species_pdf(
    bars, genomes_pdf,
    empty_text    = "No genomes assembled",
    title_suffix  = " -- polished genome",
    subtitle_flat = function(n) "Full-length genome on its reference; dotted lines: reference ends (bases beyond them are contig overhangs)",
    subtitle_seg  = function(n, p) sprintf(
      "Full-length genomes of %d segments; each panel spans the segment; dotted lines: reference ends", p),
    y_lab = "Polished genome", fill_lab = "Genome base",
    fill_scale = genome_fill,
    # No bar outlines: N and read-supported stretches can alternate every few
    # bases, and their outlines would merge into black blocks.
    outline = NA,
    ref_ends = TRUE
  )

  # ── Coverage: one page per genome ──
  message(sprintf("Plotting %d genome coverage pages -> %s", nrow(stats), coverage_pdf))
  tmp_pdfs <- character(0)
  for (i in seq_len(nrow(stats))) {
    g <- stats[i, ]
    tmp_pdf  <- tempfile(fileext = ".pdf")
    tmp_pdfs <- c(tmp_pdfs, tmp_pdf)
    draw_coverage_page(g, pos %>% filter(group_id == g$group_id),
                       contig_aln %>% filter(group_id == g$group_id),
                       min_depth, tmp_pdf)
    message(sprintf("  Plotted: %s", g$group_id))
  }
  merge_pdfs(tmp_pdfs, coverage_pdf)
}

if (mode == "contigs") {
  run_contigs(args[2], args[3], args[4], args[5],
              as.numeric(args[6]), as.numeric(args[7]))
} else {
  run_genomes(args[2], args[3], args[4], args[5], args[6], args[7],
              as.numeric(args[8]))
}
message("Done.")
