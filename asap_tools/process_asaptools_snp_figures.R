#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  library(patchwork)
  library(genbankr)
  library(Biostrings)
})

# Resolve path to local function files relative to this script
.script_path   <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
.functions_dir <- file.path(dirname(.script_path), "asap_tools_functions")
source(file.path(.functions_dir, "_extract_gene_table.R"))
source(file.path(.functions_dir, "_split_genbank_records.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: process_asaptools_snp_figures.R <rdata> <prefix> <poi_csv> [<snp_threshold>] [<snp_depth>] [<breadth_threshold>] [<aa_rdata>] [<gb1> ...]")
}

RDATA_INPUT       <- args[1]
PREFIX            <- args[2]
POI_CSV           <- args[3]
SNP_THRESHOLD     <- if (length(args) >= 4 && !args[4] %in% c("NULL", "NA", "")) as.numeric(args[4]) else 0.03
MIN_DEPTH         <- if (length(args) >= 5 && !args[5] %in% c("NULL", "NA", "")) as.numeric(args[5]) else 100
BREADTH_THRESHOLD <- if (length(args) >= 6 && !args[6] %in% c("NULL", "NA", "")) as.numeric(args[6]) else 0.8
# Trailing args (7+) are, in order and all optional: an interactive-export toggle
# (TRUE/FALSE), the amino-acid table (SNP_Amino_Acid_Table.Rdata, used to color
# SNP tiles by predicted effect), then GenBank references. Each is auto-detected
# by shape (TRUE/FALSE, *.Rdata, otherwise a GenBank path) so the toggle and AA
# table can be omitted and legacy `... breadth <gb...>` invocations still work.
.tail <- if (length(args) >= 7) args[7:length(args)] else character(0)

# Interactive HTML widgets are expensive (selfcontained ggplotly). Off by default;
# only exported when the leading tail arg is TRUE.
if (length(.tail) >= 1 && toupper(.tail[1]) %in% c("TRUE", "FALSE")) {
  EXPORT_INTERACTIVE <- toupper(.tail[1]) == "TRUE"
  .tail <- .tail[-1]
} else {
  EXPORT_INTERACTIVE <- FALSE
}

if (length(.tail) == 0) {
  AA_RDATA <- NA_character_
  GB_FILES <- character(0)
} else if (.tail[1] %in% c("NULL", "NA", "")) {          # explicit "no AA table"
  AA_RDATA <- NA_character_
  GB_FILES <- if (length(.tail) >= 2) .tail[-1] else character(0)
} else if (grepl("\\.Rdata$", .tail[1], ignore.case = TRUE)) {
  AA_RDATA <- .tail[1]
  GB_FILES <- if (length(.tail) >= 2) .tail[-1] else character(0)
} else {                                                  # legacy: GenBank files here on
  AA_RDATA <- NA_character_
  GB_FILES <- .tail
}

BREADTH_THRESHOLD_PCT <- BREADTH_THRESHOLD * 100

safe_plot <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[WARN] %s failed: %s", label, conditionMessage(e)))
  })
}

wrap_cap <- function(txt, w = 110) {
  paste(sapply(strsplit(txt, "\n")[[1]], stringr::str_wrap, width = w), collapse = "\n")
}

# --- Automatic figure sizing -------------------------------------------------
# Every JPG dimension below is derived from the number of facets (assays) and
# samples actually present, so figures stay legible whether a run has 1 amplicon
# or 50 and a handful of samples or hundreds.

# Clamp a numeric to the range [lo, hi]
clamp <- function(x, lo, hi) max(lo, min(hi, x))

# Rows/cols ggplot2::facet_wrap() will use for n panels. Mirrors ggplot's
# default heuristic (~square grid) unless ncol/nrow is pinned.
facet_grid_dims <- function(n, ncol = NULL, nrow = NULL) {
  n <- max(1, n)
  if (!is.null(ncol))      { nc <- ncol;              nr <- ceiling(n / nc) }
  else if (!is.null(nrow)) { nr <- nrow;              nc <- ceiling(n / nr) }
  else                     { nc <- ceiling(sqrt(n));  nr <- ceiling(n / nc) }
  list(nrow = max(1, nr), ncol = max(1, nc))
}

# Width/height for a facet_wrap grid: scale with the number of facet columns/rows
size_facet_grid <- function(n_facets, w_per = 5.5, h_per = 3.8,
                            w_base = 2, h_base = 2,
                            min_w = 8, max_w = 40, min_h = 6, max_h = 49,
                            ncol = NULL, nrow = NULL) {
  d <- facet_grid_dims(n_facets, ncol = ncol, nrow = nrow)
  list(width  = clamp(w_base + w_per * d$ncol, min_w, max_w),
       height = clamp(h_base + h_per * d$nrow, min_h, max_h))
}

# --- SNP effect classification (for genome-track tile coloring) --------------
# Categories, ordered by increasing impact so a position with several variants
# takes its highest-impact class (see aa_category_map build below).
AA_CATEGORY_LEVELS <- c("Non-coding", "Synonymous", "AA Change", "Frameshift")
AA_CATEGORY_COLORS <- c("Non-coding"   = "#95a5a6",  # grey
                        "Synonymous"   = "#3498db",  # blue
                        "AA Change"    = "#e74c3c",  # red (missense / in-frame indel)
                        "Frameshift"   = "#8e44ad",  # purple (out-of-frame indel)
                        "Unclassified" = "#e74c3c")  # fallback when no AA table

# Map an Amino_Acids$AA value to a category. AA is one of: "Synonymous",
# "Non-coding SNP"/"Non-coding", "Insertion/Deletion Not In-frame", or a
# "GENE:refPOSalt" string (missense or in-frame indel).
classify_aa <- function(aa) {
  dplyr::case_when(
    stringr::str_detect(aa, "Synonymous")    ~ "Synonymous",
    stringr::str_detect(aa, "Not In-frame")  ~ "Frameshift",
    stringr::str_detect(aa, "[Nn]on-coding") ~ "Non-coding",
    TRUE                                     ~ "AA Change"
  )
}

load(RDATA_INPUT)

array_info <- final_array
SNPS       <- final_snps

# Join per-position depth so plots can depth-filter; keep snp_distribution raw
# so each plot block can expand it inline (parse_snp_distribution destroys the column)
SNPS <- full_join(
  select(final_array, run, name, name_short, assay_name, position, depth),
  SNPS,
  by = c("run", "name", "name_short", "assay_name", "position" = "snp_position")
)
SNPS$snp_distribution[is.na(SNPS$snp_distribution)] <- "A=0, T=0, C=0, G=0, _=0"

# --- Handle Positions of Interest ---
if (!(POI_CSV %in% c("NA", "NULL", "", NA))) {
  genes_poi <- read.csv(POI_CSV)
  get_positions <- function(s, e) seq(min(s, e), max(s, e))
  genes_expanded <- genes_poi %>%
    rowwise() %>%
    reframe(
      X          = X,
      assay_name = seqnames,
      strand     = strand,
      type       = type,
      gene       = gene,
      position   = get_positions(start, end)
    )

  array_info <- genes_expanded %>%
    select(assay_name, position, gene) %>%
    inner_join(array_info, by = c("assay_name", "position")) %>%
    mutate(assay_name = paste(assay_name, "-", gene, sep = ""))

  SNPS <- genes_expanded %>%
    select(assay_name, position, gene) %>%
    inner_join(SNPS, by = c("assay_name", "position")) %>%
    mutate(assay_name = paste(assay_name, "-", gene, sep = ""))
}

# --- Plot: SNP Position Prevalence ---
safe_plot("SNP Prevalence", {
  # Expand snp_distribution ("A=5, T=3, ...") to one row per allele. Use
  # separate_longer_delim rather than separate()+pivot_longer(): the latter first
  # builds a (1 + max spaces)-wide frame for every row, which explodes to ~150M
  # transient rows (and >20GB) when indel-rich positions push the column count to 50+.
  SNPS_plot <- SNPS %>%
    separate_longer_delim(snp_distribution, delim = ", ") %>%
    separate_wider_delim(snp_distribution, delim = "=", names = c("Call", "n"),
                         too_many = "merge", too_few = "align_start") %>%
    mutate(snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))) %>%
    mutate(SNP = paste0(snp_reference, position, Call)) %>%
    filter(snp_reference != Call) %>%
    filter(!is.na(snp_proportion))

  SNPS_plot <- SNPS_plot %>% filter(depth >= MIN_DEPTH)
  SNPS_plot$snp_proportion[is.na(SNPS_plot$snp_proportion)] <- 0

  SNP_Plot_Data <- SNPS_plot %>%
    group_by(name, name_short, assay_name, position) %>%
    summarise(Max_SNP_proportion = max(snp_proportion), .groups = "drop")

  p_SNP <- SNP_Plot_Data %>%
    ggplot(aes(x = position, y = Max_SNP_proportion, col = name_short,
               text = paste0("Sample: ", name,
                             "<br>Position: ", position,
                             "<br>Max SNP Proportion: ", round(Max_SNP_proportion, 1), "%"))) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = SNP_THRESHOLD * 100, linetype = "dashed",
               color = "red", linewidth = 0.7, alpha = 0.8) +
    facet_wrap(~assay_name, scales = "free_x") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(y = "Max SNP Proportion (%)",
         x = "Reference Position (BP)",
         title = "SNP Locations across Reference Sequences",
         subtitle = paste0("Each point = one SNP call at that position. Red dashed line = calling threshold (",
                           SNP_THRESHOLD * 100, "%). Color = sample."),
         caption = wrap_cap("Maximum SNP proportion (%) per position per sample, for positions meeting the minimum depth threshold. Only non-reference variant calls are shown. Each dot marks a genomic position where a variant was detected.")) +
    theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

  dim_snp <- size_facet_grid(length(unique(SNP_Plot_Data$assay_name)))
  ggsave(paste0(PREFIX, "_SNP_position_prevalence.jpg"), plot = p_SNP,
         width = dim_snp$width, height = dim_snp$height, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_plot_SNP <- ggplotly(p_SNP, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_plot_SNP, paste0(PREFIX, "_SNP_position_prevalence.html"), selfcontained = TRUE)
  }

  # --- Plot: SNP Proportion Density ---
  p_snp_dist <- SNPS_plot %>%
    ggplot(aes(x = snp_proportion)) +
    geom_density(fill = "#3498db", color = "#2980b9", alpha = 0.5) +
    geom_vline(xintercept = SNP_THRESHOLD * 100, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    facet_wrap(~assay_name, scales = "free_y") +
    theme_bw() +
    labs(
      title    = "SNP Proportion Density",
      subtitle = paste0("Red dashed line = calling threshold (", SNP_THRESHOLD * 100, "%)."),
      caption  = wrap_cap("Kernel density estimate of non-reference variant proportions (%) at positions meeting the minimum depth threshold. The red dashed line marks the minimum proportion threshold used for SNP calling."),
      x = "SNP Proportion (%)", y = "Density"
    ) +
    theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

  dim_dens <- size_facet_grid(length(unique(SNPS_plot$assay_name)))
  ggsave(paste0(PREFIX, "_SNP_proportion_density.jpg"), plot = p_snp_dist,
         width = dim_dens$width, height = dim_dens$height, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_snp_dist <- ggplotly(p_snp_dist, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_snp_dist, paste0(PREFIX, "_SNP_proportion_density.html"), selfcontained = TRUE)
  }
})

# --- Plot: Strand Bias ---
safe_plot("Strand Bias", {
  if (all(c("snp_call_R1", "snp_call_R2") %in% names(SNPS))) {
    strand_data <- SNPS %>%
      filter(!is.na(snp_call_R1), !is.na(snp_call_R2)) %>%
      mutate(
        total_strand   = as.numeric(snp_call_R1) + as.numeric(snp_call_R2),
        strand_ratio   = as.numeric(snp_call_R1) / total_strand,
        # reads supporting the call (both strands) over total depth; strand_data has
        # no expanded per-allele "n" column, so the old as.numeric(n) picked up dplyr::n()
        snp_proportion = 100 * (total_strand / as.numeric(location_depth))
      ) %>%
      filter(total_strand > 0, !is.na(strand_ratio))

    p_strand <- strand_data %>%
      ggplot(aes(x = snp_proportion, y = strand_ratio, col = assay_name,
                 text = paste0("Sample: ", name,
                               "<br>Position: ", position,
                               "<br>SNP: ", snp_reference, ">", snp_call,
                               "<br>Proportion: ", round(snp_proportion, 2), "%",
                               "<br>Strand Ratio (R1/total): ", round(strand_ratio, 3)))) +
      geom_point(alpha = 0.5, size = 1) +
      geom_hline(yintercept = 0.5, col = "black", lty = "dashed", alpha = 0.7) +
      facet_wrap(~assay_name, scales = "free_x") +
      theme_bw() +
      theme(legend.position = "none") +
      labs(title = "Strand Bias Assessment",
           subtitle = "True variants cluster near 0.5 (dashed); values near 0 or 1 suggest strand-specific artifacts.",
           caption  = wrap_cap("Strand ratio (R1 / R1+R2 supporting reads) vs. SNP proportion (%) for each called variant. The dashed line at 0.5 represents equal contribution from both strands. Each point is one variant call; color indicates the amplicon."),
           x = "SNP Proportion (%)", y = "Strand Ratio (R1 / R1+R2)") +
      theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

    dim_strand <- size_facet_grid(length(unique(strand_data$assay_name)))
    ggsave(paste0(PREFIX, "_SNP_strand_bias.jpg"), plot = p_strand,
           width = dim_strand$width, height = dim_strand$height, dpi = 300)
    if (EXPORT_INTERACTIVE) {
      interactive_strand <- ggplotly(p_strand, tooltip = "text") %>% partial_bundle()
      saveWidget(interactive_strand, paste0(PREFIX, "_SNP_strand_bias.html"), selfcontained = TRUE)
    }
  } else {
    message("[INFO] Strand Bias plot skipped: snp_call_R1/snp_call_R2 columns not present.")
  }
})

# --- Plot: Base Quality ---
safe_plot("Base Quality", {
  if (all(c("snp_call_qual_mean", "snp_ref_qual_mean") %in% names(SNPS))) {
    qual_summary <- SNPS %>%
      filter(!is.na(snp_call_qual_mean), !is.na(snp_ref_qual_mean),
             snp_reference != snp_call) %>%
      mutate(snp_call_qual_mean = as.numeric(snp_call_qual_mean),
             snp_ref_qual_mean  = as.numeric(snp_ref_qual_mean)) %>%
      group_by(name, name_short, assay_name) %>%
      summarise(
        `Called SNP`     = mean(snp_call_qual_mean, na.rm = TRUE),
        `Reference Base` = mean(snp_ref_qual_mean,  na.rm = TRUE),
        .groups = "drop"
      )

    qual_long <- qual_summary %>%
      pivot_longer(c(`Called SNP`, `Reference Base`),
                   names_to = "Type", values_to = "Mean_Quality") %>%
      mutate(Type = factor(Type, levels = c("Called SNP", "Reference Base")))

    p_qual <- ggplot() +
      geom_line(data = qual_long,
                aes(x = Type, y = Mean_Quality, group = name_short),
                linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_point(data = qual_long,
                 aes(x = Type, y = Mean_Quality, color = Type,
                     text = paste0("Sample: ", name,
                                   "<br>Assay: ", assay_name,
                                   "<br>", Type, ": ", round(Mean_Quality, 1))),
                 size = 3) +
      facet_wrap(~assay_name) +
      scale_color_manual(values = c("Called SNP" = "#e74c3c", "Reference Base" = "#3498db")) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(title = "Base Quality: Called SNP vs Reference Base (per sample mean)",
           subtitle = "Red = Called SNP mean quality, Blue = Reference Base mean quality. Each dashed line = one sample.",
           caption  = wrap_cap("Mean base quality score per sample per amplicon for two categories: Called SNP (red) = mean quality of bases at positions where a variant was called; Reference Base (blue) = mean quality of bases at the same positions supporting the reference allele. Each dashed line connects the two values for one sample."),
           x = NULL, y = "Mean Base Quality Score") +
      theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

    # x is only two categories per facet, so facets can be narrower than the default
    dim_qual <- size_facet_grid(length(unique(qual_long$assay_name)), w_per = 3.5)
    ggsave(paste0(PREFIX, "_SNP_base_quality.jpg"), plot = p_qual,
           width = dim_qual$width, height = dim_qual$height, dpi = 300)
    if (EXPORT_INTERACTIVE) {
      interactive_qual <- ggplotly(p_qual, tooltip = "text") %>% partial_bundle()
      saveWidget(interactive_qual, paste0(PREFIX, "_SNP_base_quality.html"), selfcontained = TRUE)
    }
  } else {
    message("[INFO] Base Quality plot skipped: snp_call_qual_mean/snp_ref_qual_mean columns not present.")
  }
})

# --- Plot: Genome Track (A: genes / B: SNP density / C: sample heatmap) per GenBank file ---
safe_plot("Genome Track", {
  if (length(GB_FILES) == 0) {
    message("[INFO] Genome Track skipped: no GenBank files provided.")
    return(invisible(NULL))
  }

  # See SNP Prevalence block: separate_longer_delim avoids the ~150M-row / >20GB
  # blowup that separate()+pivot_longer() causes on indel-rich distributions.
  snp_local <- SNPS %>%
    separate_longer_delim(snp_distribution, delim = ", ") %>%
    separate_wider_delim(snp_distribution, delim = "=", names = c("Call", "n"),
                         too_many = "merge", too_few = "align_start") %>%
    mutate(snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))) %>%
    filter(snp_reference != Call, !is.na(snp_proportion), depth >= MIN_DEPTH)

  snp_pos_data <- snp_local %>%
    group_by(name, name_short, assay_name, position) %>%
    summarise(Max_SNP_proportion = max(snp_proportion), .groups = "drop")

  # Per-(assay, position) SNP effect category from the amino-acid table, used to
  # color the panel-C tiles. Combo (codon-merge) tokens contribute every position
  # they contain; a position keeps its highest-impact category across variants.
  aa_category_map <- NULL
  if (!is.na(AA_RDATA) && file.exists(AA_RDATA)) {
    aa_env <- new.env()
    load(AA_RDATA, envir = aa_env)
    if (exists("Amino_Acids", envir = aa_env) && nrow(aa_env$Amino_Acids) > 0) {
      aa_category_map <- aa_env$Amino_Acids %>%
        transmute(assay_name,
                  category = classify_aa(AA),
                  position = stringr::str_extract_all(SNP, "[0-9]+")) %>%
        tidyr::unnest(position) %>%
        mutate(position = as.integer(position),
               rank     = match(category, AA_CATEGORY_LEVELS)) %>%
        group_by(assay_name, position) %>%
        summarise(aa_category = AA_CATEGORY_LEVELS[max(rank)], .groups = "drop")
      message(sprintf("[INFO] Genome Track: loaded %d AA-annotated positions for tile coloring.",
                      nrow(aa_category_map)))
    }
  } else {
    message("[INFO] Genome Track: no amino-acid table supplied; SNP tiles not colored by effect.")
  }

  MAX_GENOME_LEN <- 2e6  # bp; this panel is designed for viral/phage-scale references

  for (gb_file in GB_FILES) {
    # Split multi-record references into per-contig files (genbankr cannot read a
    # multi-LOCUS file). Single-contig files yield exactly one record.
    records <- tryCatch(split_genbank_records(gb_file), error = function(e) {
      message(sprintf("[WARN] Genome Track: could not split %s: %s",
                      basename(gb_file), conditionMessage(e)))
      NULL
    })
    if (is.null(records)) next

    for (rec_i in seq_len(nrow(records))) {
      contig_path <- records$path[rec_i]
    tryCatch({
      ref    <- suppressWarnings(genbankr::readGenBank(contig_path))
      acc_id <- ref@accession
      fb     <- sub("\\.[^.]+$", "", basename(gb_file))

      gene_df <- extract_gene_table(ref) %>%
        mutate(
          gene_label  = as.character(gene),
          label_angle = ifelse((end - start) < nchar(as.character(ref@sequence)) / 20, 90, 0)
        ) %>%
        filter(!is.na(start), !is.na(end)) %>%
        mutate(gene_color = scales::hue_pal()(n()))

      genome_len  <- nchar(as.character(ref@sequence))

      if (genome_len > MAX_GENOME_LEN) {
        message(sprintf(
          "[INFO] Genome Track: %s is %s bp (exceeds %s bp limit), skipping — panel is not designed for bacterial/eukaryotic-scale genomes.",
          acc_id, format(genome_len, big.mark = ","), format(MAX_GENOME_LEN, big.mark = ",")
        ))
        next
      }

      gene_bounds <- sort(unique(c(gene_df$start, gene_df$end)))
      x_lim       <- c(0, genome_len)

      snp_match <- snp_pos_data %>%
        filter(grepl(acc_id, assay_name, fixed = TRUE) |
               grepl(fb,     assay_name, fixed = TRUE))

      cov_match <- array_info %>%
        filter(grepl(acc_id, assay_name, fixed = TRUE) |
               grepl(fb,     assay_name, fixed = TRUE)) %>%
        select(name, name_short, position, depth)

      if (nrow(gene_df) == 0 || nrow(snp_match) == 0) {
        message(sprintf("[INFO] Genome Track: insufficient data for %s, skipping.", acc_id))
        next
      }

      cov_tiles <- if (nrow(cov_match) > 0) {
        cov_match %>%
          filter(depth > 0) %>%
          arrange(name, position) %>%
          group_by(name) %>%
          mutate(run_id = cumsum(c(1L, diff(position)) > 1L)) %>%
          group_by(name, name_short, run_id) %>%
          summarise(
            cov_xmid  = (min(position) + max(position)) / 2,
            cov_width = max(position) - min(position) + 1L,
            .groups = "drop"
          )
      } else {
        tibble(name = character(), name_short = character(), cov_xmid = numeric(), cov_width = integer())
      }

      snp_width <- max(1, genome_len / 1000)

      p_genes <- gene_df %>%
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = gene_label),
                  color = "white", linewidth = 0.4, alpha = 0.55) +
        geom_text(aes(x = (start + end) / 2, y = 0.5,
                      label = gene_label, angle = label_angle),
                  size = 2.8, color = "white", fontface = "bold") +
        scale_fill_manual(values = setNames(gene_df$gene_color, gene_df$gene_label)) +
        scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = scales::comma) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(b = 2)) +
        labs(title = paste0("Genome Track: ", acc_id))

      safe_name <- gsub("[^A-Za-z0-9]", "_", acc_id)

      # SNP effect per position for THIS contig (matched like snp_match above,
      # then keyed on position). NULL when no amino-acid table was supplied.
      aa_cat_contig <- if (!is.null(aa_category_map)) {
        aa_category_map %>%
          filter(grepl(acc_id, assay_name, fixed = TRUE) |
                 grepl(fb,     assay_name, fixed = TRUE)) %>%
          group_by(position) %>%
          summarise(aa_category = AA_CATEGORY_LEVELS[max(match(aa_category, AA_CATEGORY_LEVELS))],
                    .groups = "drop") %>%
          mutate(position = as.numeric(position))
      } else NULL

      render_genome_track <- function(snp_match_i, cov_tiles_i, suffix, caption_extra = "") {
        if (nrow(snp_match_i) == 0) {
          message(sprintf("[INFO] Genome Track (%s): no samples remain for %s, skipping.", suffix, acc_id))
          return(invisible(NULL))
        }

        sample_order_i <- sort(unique(as.character(snp_match_i$name_short)))
        n_samples_i    <- length(sample_order_i)
        snp_match_i    <- mutate(snp_match_i, name_short = factor(as.character(name_short), levels = sample_order_i))
        cov_tiles_i    <- cov_tiles_i %>%
          filter(as.character(name_short) %in% sample_order_i) %>%
          mutate(name_short = factor(as.character(name_short), levels = sample_order_i))

        # Attach the SNP effect category used to fill panel-C tiles. Positions
        # with no AA annotation default to Non-coding; with no table at all,
        # every tile is a single "Unclassified" colour (legacy red).
        if (!is.null(aa_cat_contig)) {
          snp_match_i <- snp_match_i %>%
            mutate(position = as.numeric(position)) %>%
            left_join(aa_cat_contig, by = "position")
          snp_match_i$aa_category[is.na(snp_match_i$aa_category)] <- "Non-coding"
        } else {
          snp_match_i$aa_category <- "Unclassified"
        }
        snp_match_i$aa_category <- factor(snp_match_i$aa_category,
                                          levels = names(AA_CATEGORY_COLORS))

        p_density_i <- snp_match_i %>%
          ggplot(aes(x = as.numeric(position))) +
          geom_density(fill = "grey80", alpha = 0.5, color = "black",
                       bw = max(50, genome_len / 500)) +
          geom_vline(xintercept = gene_bounds, color = "grey60", alpha = 0.35, linewidth = 0.3) +
          scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = scales::comma) +
          theme_minimal() +
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                axis.ticks.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks.y = element_blank(), panel.grid = element_blank(),
                plot.margin = margin(t = 0, b = 0)) +
          labs(y = "SNP Density")

        p_heat_i <- ggplot() +
          geom_tile(data = cov_tiles_i,
                    aes(x = cov_xmid, y = name_short, width = cov_width, height = 0.85),
                    fill = "grey25", alpha = 0.35) +
          geom_vline(xintercept = gene_bounds, color = "grey60", alpha = 0.35, linewidth = 0.3) +
          geom_tile(data = snp_match_i,
                    aes(x = as.numeric(position), y = name_short,
                        fill = aa_category, alpha = Max_SNP_proportion,
                        text = paste0("Sample: ", name,
                                      "<br>Position: ", position,
                                      "<br>Max SNP %: ", round(Max_SNP_proportion, 1), "%",
                                      "<br>Effect: ", aa_category)),
                    width = snp_width, height = 0.85) +
          scale_fill_manual(values = AA_CATEGORY_COLORS, drop = TRUE,
                            name = "SNP effect", na.value = "grey70") +
          scale_alpha_continuous(range = c(0.45, 1), limits = c(0, 100),
                                 name = "Max SNP %") +
          scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = scales::comma) +
          theme_minimal() +
          theme(axis.text.x  = element_text(angle = 45, hjust = 1),
                axis.text.y  = element_text(size = 7),
                panel.grid.major.x = element_line(color = "grey95", linewidth = 0.1),
                panel.grid.minor   = element_blank(),
                panel.grid.major.y = element_blank(),
                plot.margin = margin(t = 0)) +
          labs(x = "Reference Position (BP)",
               y = paste0("Samples (n=", n_samples_i, ")"),
               caption = wrap_cap(paste0(
                 "(A) CDS gene annotations parsed from ", basename(gb_file), ". ",
                 "(B) Kernel density estimate of SNP positions across the genome. ",
                 "(C) Per-sample SNP locations colored by predicted effect ",
                 "(red = AA change, purple = frameshift, blue = synonymous, grey = non-coding), ",
                 "overlaid on covered regions (grey shading). Tile opacity scales with SNP proportion. ",
                 "Grey shading = positions with depth > 0. Only positions meeting minimum depth are eligible for SNP calling.",
                 caption_extra
               ))) +
          theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

        combined_track_i <- p_genes / p_density_i / p_heat_i +
          plot_layout(heights = c(1, 2, max(3, n_samples_i * 0.3))) +
          plot_annotation(
            tag_levels = "A", tag_suffix = ".",
            theme = theme(
              plot.tag    = element_text(face = "bold", size = 14),
              plot.margin = margin(l = 5, r = 5, t = 5, b = 5)
            )
          )

        # cap under ggplot2's 50-inch ggsave limit for very high sample counts
        track_height_i <- min(49, max(8, 3.5 + n_samples_i * 0.3))

        ggsave(paste0(PREFIX, "_SNP_", safe_name, suffix, ".jpg"),
               plot = combined_track_i, width = 18, height = track_height_i, dpi = 300)

        if (EXPORT_INTERACTIVE) tryCatch({
          interactive_track_i <- ggplotly(p_heat_i, tooltip = "text") %>% partial_bundle()
          saveWidget(interactive_track_i,
                     paste0(PREFIX, "_SNP_", safe_name, suffix, ".html"),
                     selfcontained = TRUE)
        }, error = function(e) {
          message(sprintf("[WARN] Genome Track (%s) HTML export failed for %s: %s",
                          suffix, acc_id, conditionMessage(e)))
        })
      }

      render_genome_track(snp_match, cov_tiles, "_genome_track")

      breadth_match <- final_asap %>%
        filter(grepl(acc_id, assay_name, fixed = TRUE) |
               grepl(fb,     assay_name, fixed = TRUE)) %>%
        mutate(breadth = as.numeric(breadth))

      passing_samples <- breadth_match %>%
        filter(breadth > BREADTH_THRESHOLD_PCT) %>%
        pull(name) %>%
        unique()

      render_genome_track(
        filter(snp_match, as.character(name) %in% passing_samples),
        filter(cov_tiles, as.character(name) %in% passing_samples),
        "_genome_track_filtered",
        caption_extra = paste0(
          " Samples restricted to those with breadth of coverage > ",
          round(BREADTH_THRESHOLD_PCT, 1), "%."
        )
      )

    }, error = function(e) {
      message(sprintf("[WARN] Genome Track for %s failed: %s",
                      basename(gb_file), conditionMessage(e)))
    })
    }  # end per-contig record loop
  }
})
