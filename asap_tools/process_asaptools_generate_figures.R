#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(plotly)
  library(htmlwidgets)
  library(zoo)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: process_asaptools_generate_figures.R <rdata> <prefix> <poi_csv> [<snp_threshold>] [<snp_depth>] [<interactive>]")
}

RDATA_INPUT   <- args[1]
PREFIX        <- args[2]
POI_CSV       <- args[3]
SNP_THRESHOLD <- if (length(args) >= 4 && !args[4] %in% c("NULL", "NA", "")) as.numeric(args[4]) else 0.03
MIN_DEPTH     <- if (length(args) >= 5 && !args[5] %in% c("NULL", "NA", "")) as.numeric(args[5]) else 100
# Interactive HTML widgets (selfcontained ggplotly) are expensive; off by default,
# exported only when arg 6 is TRUE. JPGs are always written.
EXPORT_INTERACTIVE <- length(args) >= 6 && toupper(args[6]) == "TRUE"

# Wraps a plot block so a single failure doesn't abort all plots
safe_plot <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[WARN] %s failed: %s", label, conditionMessage(e)))
  })
}

# Word-wrap caption text; preserves existing \n paragraph breaks
wrap_cap <- function(txt, w = 110) {
  paste(sapply(strsplit(txt, "\n")[[1]], stringr::str_wrap, width = w), collapse = "\n")
}

# --- Automatic figure sizing -------------------------------------------------
# All ggsave() dimensions below are derived from the number of facets (assays)
# and samples actually present in each plot, so figures stay legible whether a
# run has 2 samples or 200 and 1 amplicon or 50.

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

load(RDATA_INPUT)

array_info <- final_array

# Join metadata
array_info <- left_join(array_info, select(final_asap, run, assay_name, name, name_short, amplicon_reads, avg_depth, breadth))

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
}

# Dynamic plot data reduction
target_points <- 10000
total_span    <- length(unique(array_info$position))
dynamic_k     <- max(1, floor(total_span / target_points))

# Rolling average dataset for coverage/N-read plots
array_avg <- array_info %>%
  group_by(name, assay_name) %>%
  arrange(position) %>%
  mutate(
    depth_avg    = rollmean(depth, k = dynamic_k, fill = NA),
    n_reads_prop = rollmean(100 * (n_reads / depth), k = dynamic_k, fill = NA)
  ) %>%
  filter(!is.na(depth_avg)) %>%
  filter(row_number() %% dynamic_k == 0)

# --- Plot 1: Coverage Depth ---
safe_plot("Coverage Depth", {
  p_cov <- array_avg %>%
    ggplot(aes(x = position, y = depth_avg, col = name_short, group = name_short,
               text = paste0("Sample: ", name,
                             "<br>~Position: ", position,
                             "<br>Mean Depth (10bp): ", round(depth_avg, 1)))) +
      geom_line(alpha = 0.7) +
      geom_hline(yintercept = MIN_DEPTH, col = "Red", lty = "dashed", alpha = 0.6) +
      facet_wrap(~assay_name, scales = "free_x", ncol = 1) +
      scale_y_log10() +
      theme_bw() +
      theme(legend.position = "none") +
      labs(
        title    = "Reference Depth of Coverage",
        subtitle = paste0("Red dashed line = minimum depth (", MIN_DEPTH, "x) required for SNP calls. Each line = one sample."),
        caption  = wrap_cap("Mean sequencing depth (log10 scale) per reference position, smoothed over a rolling window. Each line represents one sample. Depth is derived from reads aligned to the reference. Positions below the red dashed line do not meet the minimum depth threshold for SNP calling."),
        y = paste0("Mean Coverage Depth (", dynamic_k, "bp window, log10 scaled)"),
        x = "Reference Position (BP)"
      ) +
      theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

  # ncol=1 stack: one facet row per assay -> height grows with facet count
  n_facets_cov <- length(unique(array_avg$assay_name))
  cov_h <- clamp(3 + 2.3 * n_facets_cov, 6, 49)
  ggsave(paste0(PREFIX, "_QC_coverage_depth.jpg"), plot = p_cov, width = 12, height = cov_h, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_plot_coverage <- ggplotly(p_cov, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_plot_coverage, paste0(PREFIX, "_QC_coverage_depth.html"), selfcontained = TRUE)
  }
})

# --- Plot 2: N Read Proportion ---
safe_plot("N Read Proportion", {
  p_n <- array_avg %>%
    ggplot(aes(x = position, y = n_reads_prop, col = name_short, group = name_short,
               text = paste0("Sample: ", name,
                             "<br>~Position: ", position,
                             "<br>Proporion 'N' Reads (10bp window): ", round(n_reads_prop, 1)))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~assay_name, scales = "free_x", ncol = 1) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(y = paste0("Proportion 'N' Reads (", dynamic_k, "bp window)"),
         x = "Reference Position (BP)",
         title = "Percent Reads with 'N's across Reference",
         subtitle = "'N' bases accumulate from base quality masking, primer masking, and SMOR deduplication.",
         caption = wrap_cap("Percentage of reads carrying 'N' bases at each reference position, smoothed over a rolling window. 'N' bases are contributed by base quality masking, primer sequence masking, and SMOR deduplication masking.")) +
    theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

  # ncol=1 stack: one facet row per assay -> height grows with facet count
  n_facets_n <- length(unique(array_avg$assay_name))
  n_h <- clamp(3 + 2.3 * n_facets_n, 6, 49)
  ggsave(paste0(PREFIX, "_QC_n_reads_proportion.jpg"), plot = p_n, width = 12, height = n_h, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_plot_n_reads <- ggplotly(p_n, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_plot_n_reads, paste0(PREFIX, "_QC_n_reads_proportion.html"), selfcontained = TRUE)
  }
})

# --- Plot 3: Breadth of Coverage Heatmap ---
safe_plot("Breadth Heatmap", {
  breadth_data <- final_asap %>%
    select(name, name_short, assay_name, breadth) %>%
    mutate(breadth = as.numeric(breadth))

  p_breadth <- breadth_data %>%
    ggplot(aes(x = assay_name, y = name_short, fill = breadth,
               text = paste0("Sample: ", name,
                             "<br>Assay: ", assay_name,
                             "<br>Breadth: ", round(breadth, 1), "%"))) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "#2ecc71", limits = c(0, 100),
                        name = paste0("Breadth\n(%, ≥", MIN_DEPTH, "x)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Breadth of Coverage per Sample and Amplicon",
         subtitle = paste0("Proportion of amplicon positions with ASAP depth ≥ ", MIN_DEPTH,
                           "x. White = no coverage, green = full coverage."),
         caption  = wrap_cap(paste0(
           "Breadth of coverage per sample and amplicon as computed by ASAP, reflecting positions ",
           "with at least ", MIN_DEPTH, "x depth (--depth) from reads passing all pipeline filters ",
           "(primer matching, identity threshold, SMOR deduplication). ",
           "Values range from 0% (no positions meet threshold) to 100% (all positions meet threshold)."
         )),
         x = "Assay", y = "Sample") +
    theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

  # Heatmap: assays on x, samples on y -> width per assay tile, height per sample row
  n_assays_br  <- length(unique(breadth_data$assay_name))
  n_samples_br <- length(unique(breadth_data$name_short))
  br_w <- clamp(4 + 0.55 * n_assays_br,  8, 40)
  br_h <- clamp(3 + 0.35 * n_samples_br, 6, 40)
  ggsave(paste0(PREFIX, "_Breadth_coverage_heatmap.jpg"), plot = p_breadth, width = br_w, height = br_h, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_breadth <- ggplotly(p_breadth, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_breadth, paste0(PREFIX, "_Breadth_coverage_heatmap.html"), selfcontained = TRUE)
  }
})

# --- Plot 4: Alignment Summary (counts + percentage panels) ---
safe_plot("Alignment Summary", {
  align_data <- final_asap %>%
    select(name, name_short, total_reads, trimmed_reads, mapped_reads, unassigned_reads, unmapped_reads) %>%
    distinct() %>%
    mutate(
      across(c(total_reads, trimmed_reads, mapped_reads, unassigned_reads, unmapped_reads), as.numeric),
      # unmapped_reads (pysam .unmapped, from newBamProcessor.py) is the TOTAL unmapped
      # count; unassigned_reads (pysam .nocoordinate) is the subset with no alignment
      # coordinate at all (neither mate mapped anywhere) -- already included inside
      # unmapped_reads, not a separate pool. Subtract it back out so the two are
      # mutually exclusive and the stacked categories actually sum to total_reads.
      unmapped_only = pmax(0, unmapped_reads - unassigned_reads),
      lost_fastp    = pmax(0, total_reads - trimmed_reads),
      lost_other    = pmax(0, trimmed_reads - mapped_reads - unmapped_reads)
    ) %>%
    select(name, name_short, total_reads,
           "Aligned"          = mapped_reads,
           "Unassigned"       = unassigned_reads,
           "Unmapped"         = unmapped_only,
           "Other"            = lost_other,
           "Removed by FastP" = lost_fastp) %>%
    pivot_longer(-c(name, name_short, total_reads), names_to = "Category", values_to = "Reads") %>%
    mutate(
      Category = factor(Category, levels = c("Aligned", "Unassigned", "Unmapped",
                                             "Other", "Removed by FastP")),
      Percent  = 100 * Reads / total_reads
    ) %>%
    filter(!is.na(Reads), Reads > 0)

  align_colours <- c("Aligned"          = "#2ecc71",
                     "Unassigned"       = "#f1c40f",
                     "Unmapped"         = "#e74c3c",
                     "Other"            = "#95a5a6",
                     "Removed by FastP" = "#bdc3c7")

  align_base_theme <- theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  p_align <- align_data %>%
    ggplot(aes(x = name_short, y = Reads, fill = Category,
               text = paste0("Sample: ", name,
                             "<br>Category: ", Category,
                             "<br>Reads: ", scales::comma(Reads)))) +
    geom_col() +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = align_colours) +
    align_base_theme +
    labs(title = "Read Counts", x = "Sample", y = "Read Count", fill = NULL)

  p_align_pct <- align_data %>%
    ggplot(aes(x = name_short, y = Percent, fill = Category,
               text = paste0("Sample: ", name,
                             "<br>Category: ", Category,
                             "<br>Percent: ", round(Percent, 1), "%"))) +
    geom_col() +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
    scale_fill_manual(values = align_colours) +
    align_base_theme +
    labs(title = "Percentage of Reads", x = "Sample", y = "Percentage (%)", fill = NULL)

  p_align_combined <- (p_align / p_align_pct) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  p_align_combined <- p_align_combined + plot_annotation(
    title      = "Alignment Summary",
    subtitle   = "Bar height = total reads (pre-trim). Colors show read fate through FastP trimming and alignment.",
    tag_levels = "A",
    caption    = wrap_cap(paste(
      "(A) Total read count per sample stacked by alignment outcome. Bar height = total reads before FastP trimming. Categories: Aligned = primary reads mapped to the reference; Unassigned = reads with no alignment coordinate at all (neither mate mapped anywhere); Unmapped = reads flagged unmapped but still assigned a coordinate because their mate mapped nearby (mate-rescued placement) — a subset of the total unmapped count not already captured by Unassigned; Removed by FastP = discarded during adapter/quality trimming; Other = reads not accounted for by the above categories.",
      "(B) Same data expressed as a percentage of total pre-trim reads per sample.",
      sep = "\n"
    ), w = 150),
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.4),
      plot.tag     = element_text(face = "bold", size = 12)
    )
  )

  # Two stacked panels, samples on x -> width grows with sample count;
  # height is two panels of bars plus room for the shared legend + caption.
  n_samples_al <- length(unique(align_data$name_short))
  al_w <- clamp(6 + 0.45 * n_samples_al, 12, 40)
  al_h <- clamp(al_w * 0.85,             14, 30)
  ggsave(paste0(PREFIX, "_QC_alignment_summary.jpg"), plot = p_align_combined,
         width = al_w, height = al_h, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_align <- subplot(
      ggplotly(p_align,     tooltip = "text") %>% partial_bundle(),
      ggplotly(p_align_pct, tooltip = "text") %>% partial_bundle(),
      nrows = 1, shareY = FALSE, titleX = TRUE, titleY = TRUE
    ) %>% layout(title = "Alignment Summary")
    saveWidget(interactive_align, paste0(PREFIX, "_QC_alignment_summary.html"), selfcontained = TRUE)
  }
})

# --- Plot 5: Read Funnel (counts + percentage panels) ---
safe_plot("Read Funnel", {
  # NOTE: no_primer_reads is intentionally NOT a loss category. With
  # primer_only=false (the default), reads where no primer was detected are only
  # primer-masked (bases -> N); they still flow through to amplicon_reads. The
  # true conservation is aligned_reads = amplicon_reads + identity_discarded +
  # smor_pairs_dropped + residual, so counting no_primer_reads here would
  # double-count them (once as "lost", once inside Final Reads). SMOR columns are
  # NA when SMOR is off, so coalesce to 0 to keep the arithmetic well-defined.
  funnel_data <- final_asap %>%
    select(name, name_short, assay_name,
           aligned_reads,
           identity_discarded,
           smor_pairs_dropped,
           amplicon_reads) %>%
    mutate(
      across(c(aligned_reads, identity_discarded,
               smor_pairs_dropped, amplicon_reads), as.numeric),
      lost_identity = coalesce(identity_discarded, 0),
      lost_smor     = coalesce(smor_pairs_dropped, 0),
      kept          = coalesce(amplicon_reads, 0),
      lost_other    = pmax(0, coalesce(aligned_reads, 0) - kept - lost_smor - lost_identity)
    )

  fate_levels  <- c("Lost: Identity Filter", "Lost: SMOR",
                    "Lost: Other", "Final Reads")
  fate_colours <- c("Lost: Identity Filter"   = "#e67e22",
                    "Lost: SMOR"              = "#f1c40f",
                    "Lost: Other"             = "#95a5a6",
                    "Final Reads"             = "#2ecc71")

  funnel_long <- funnel_data %>%
    select(name, name_short, assay_name,
           "Lost: Identity Filter" = lost_identity,
           "Lost: SMOR"            = lost_smor,
           "Lost: Other"           = lost_other,
           "Final Reads"           = kept) %>%
    pivot_longer(-c(name, name_short, assay_name), names_to = "Fate", values_to = "Reads") %>%
    mutate(Fate = factor(Fate, levels = fate_levels)) %>%
    filter(!is.na(Reads))

  funnel_totals <- funnel_long %>%
    group_by(name, assay_name) %>%
    summarise(total_aligned = sum(Reads, na.rm = TRUE), .groups = "drop")

  funnel_long <- funnel_long %>%
    left_join(funnel_totals, by = c("name", "assay_name")) %>%
    mutate(Percent = if_else(total_aligned == 0, 0, 100 * Reads / total_aligned))

  funnel_base_theme <- theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  p_funnel <- funnel_long %>%
    filter(Reads > 0) %>%
    ggplot(aes(x = name_short, y = Reads, fill = Fate,
               text = paste0("Sample: ", name,
                             "<br>Fate: ", Fate,
                             "<br>Reads: ", scales::comma(Reads)))) +
    geom_col() +
    facet_wrap(~assay_name, scales = "free_x") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = fate_colours) +
    funnel_base_theme +
    labs(title = "Read Counts", x = "Sample", y = "Read Count", fill = NULL)

  p_funnel_pct <- funnel_long %>%
    ggplot(aes(x = name_short, y = Percent, fill = Fate,
               text = paste0("Sample: ", name,
                             "<br>Fate: ", Fate,
                             "<br>Percent: ", round(Percent, 1), "%"))) +
    geom_col() +
    facet_wrap(~assay_name, scales = "free_x") +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
    scale_fill_manual(values = fate_colours) +
    funnel_base_theme +
    labs(title = "Percentage of Aligned Reads", x = "Sample", y = "Percentage (%)", fill = NULL)

  p_funnel_combined <- (p_funnel / p_funnel_pct) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  p_funnel_combined <- p_funnel_combined + plot_annotation(
    title      = "Read Fate per Amplicon",
    subtitle   = "Bar height = reads aligned to amplicon. Colors show where reads were removed or retained.",
    tag_levels = "A",
    caption    = wrap_cap(paste(
      "(A) Read counts per amplicon per sample stacked by filtering fate. Bar height = total reads aligned to the amplicon.",
      "(B) Same data as percentage of aligned reads per amplicon.",
      "Categories — Lost:Identity: read pair did not meet percent-identity threshold; Lost:SMOR: duplicate pair removed by SMOR deduplication; Lost:Other: aligned reads not accounted for by the above filters; Final Reads: reads passing all filters and counted as amplicon_reads. Reads with no detected primer are primer-masked but still retained (they pass through to Final Reads), so they are not shown as a loss.",
      sep = "\n"
    ), w = 120),
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.4),
      plot.tag     = element_text(face = "bold", size = 12)
    )
  )

  # facet_wrap grid stacked twice (counts over percentage). Width scales with the
  # facet columns and the samples shown per facet; height with the facet rows,
  # doubled for the two patchwork panels, plus room for the legend + caption.
  n_samples_fn <- length(unique(funnel_long$name_short))
  n_assays_fn  <- length(unique(funnel_long$assay_name))
  fdim <- facet_grid_dims(n_assays_fn)
  fn_w <- clamp(fdim$ncol * (2.5 + 0.4 * n_samples_fn), 12, 49)
  fn_h <- clamp(2 * fdim$nrow * 3.2 + 3,                12, 49)
  ggsave(paste0(PREFIX, "_QC_read_funnel.jpg"), plot = p_funnel_combined,
         width = fn_w, height = fn_h, dpi = 300)
  if (EXPORT_INTERACTIVE) {
    interactive_funnel <- subplot(
      ggplotly(p_funnel,     tooltip = "text") %>% partial_bundle(),
      ggplotly(p_funnel_pct, tooltip = "text") %>% partial_bundle(),
      nrows = 2, shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE
    ) %>% layout(title = "Read Fate per Amplicon")
    saveWidget(interactive_funnel, paste0(PREFIX, "_QC_read_funnel.html"), selfcontained = TRUE)
  }
})
