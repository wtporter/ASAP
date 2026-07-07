#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)
library(plotly)
library(htmlwidgets)
library(patchwork)

# Resolve path to local function files relative to this script
.script_path   <- normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
.functions_dir <- file.path(dirname(.script_path), "asap_tools_functions")
source(file.path(.functions_dir, "_shorten_sample_names.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: process_fastp_panel.R <prefix> <json1> [<json2> ...]")
}

PREFIX    <- args[1]
json_files <- args[-1]

# Wraps a plot block so a single panel failure doesn't abort all others
safe_plot <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[WARN] fastp panel '%s' failed: %s", label, conditionMessage(e)))
    NULL
  })
}

# --- Parse JSONs ---
parse_fastp_json <- function(path) {
  sample_id <- sub("\\.fastplong\\.json$", "", sub("\\.fastp\\.json$", "", basename(path)))
  d <- fromJSON(path, simplifyVector = TRUE)

  bf  <- d$summary$before_filtering
  af  <- d$summary$after_filtering
  fr  <- d$filtering_result
  dup <- d$duplication$rate
  ins <- d$insert_size$peak
  ac  <- d$adapter_cutting

  adapter_trimmed_reads <- ac$adapter_trimmed_reads %||% 0L
  adapter_trimmed_bases <- ac$adapter_trimmed_bases %||% 0L

  tibble(
    sample                    = sample_id,
    total_reads_before        = bf$total_reads,
    total_reads_after         = af$total_reads,
    reads_removed             = bf$total_reads - af$total_reads,
    q20_before                = bf$q20_rate * 100,
    q20_after                 = af$q20_rate  * 100,
    q30_before                = bf$q30_rate  * 100,
    q30_after                 = af$q30_rate  * 100,
    gc_before                 = bf$gc_content * 100,
    gc_after                  = af$gc_content * 100,
    read_len_before           = bf$read1_mean_length,
    read_len_after            = af$read1_mean_length,
    low_quality_reads         = fr$low_quality_reads %||% 0L,
    too_short_reads           = fr$too_short_reads   %||% 0L,
    too_long_reads            = fr$too_long_reads    %||% 0L,
    too_many_N_reads          = fr$too_many_N_reads  %||% 0L,
    duplication_rate          = dup * 100,
    insert_size_peak          = if (!is.null(ins)) ins else NA_real_,
    pct_reads_adapter_trimmed = adapter_trimmed_reads / bf$total_reads * 100,
    pct_bases_adapter_trimmed = adapter_trimmed_bases / bf$total_bases * 100
  )
}

# %||% operator (NULL coalescing)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Extract per-cycle mean quality for R1 (and R2 if paired) before/after filtering
parse_fastp_cycles <- function(path) {
  sample_id <- sub("\\.fastplong\\.json$", "", sub("\\.fastp\\.json$", "", basename(path)))
  d <- fromJSON(path, simplifyVector = TRUE)

  reads <- list(
    list(key = "read1_before_filtering", read = "R1", timing = "Before"),
    list(key = "read1_after_filtering",  read = "R1", timing = "After")
  )
  if (!is.null(d[["read2_before_filtering"]])) {
    reads <- c(reads,
      list(list(key = "read2_before_filtering", read = "R2", timing = "Before")),
      list(list(key = "read2_after_filtering",  read = "R2", timing = "After"))
    )
  }

  map_dfr(reads, function(r) {
    quals <- d[[r$key]]$quality_curves$mean
    if (is.null(quals)) return(NULL)
    tibble(
      sample       = sample_id,
      read         = r$read,
      timing       = r$timing,
      cycle        = seq_along(quals),
      mean_quality = quals
    )
  })
}

qc <- map_dfr(json_files, safely(parse_fastp_json, otherwise = NULL)) %>%
  { bind_rows(.$result) }

if (nrow(qc) == 0) stop("No fastp JSON files could be parsed.")

qc$sample_short <- shorten_sample_names(qc$sample)

n_samples <- nrow(qc)

# Each 2-column tier (and each facet within panel C) only gets half the total
# figure width, so width must scale with n_samples to keep rotated labels from
# overlapping. Grows linearly at 18in per 100 samples so massive cohorts keep
# getting wider instead of plateauing early; floored at 18in (today's fixed
# width) so cohorts under ~100 samples aren't shrunk below the old baseline.
plot_width <- max(18, 18 * (n_samples / 100))

axis_text_size <- case_when(
  n_samples <= 15 ~ 9,
  n_samples <= 40 ~ 7,
  n_samples <= 80 ~ 5.5,
  TRUE             ~ 4
)
qc_axis_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size))

# Zooms each panel's y-axis to its own observed range (padded 2.5% beyond the
# min/max) instead of the ggplot default (bars from 0, lines auto-expanded),
# so cross-sample differences are easy to see even when every sample is
# clustered in a narrow band -- the tradeoff being bars no longer start at a
# true zero baseline.
scaled_ylim <- function(values) {
  rng <- range(values, na.rm = TRUE)
  c(rng[1] * 0.975, rng[2] * 1.025)
}

cycles <- map_dfr(json_files, safely(parse_fastp_cycles, otherwise = NULL)) %>%
  { bind_rows(.$result) }

reads_data <- qc %>%
  select(sample, sample_short, Before = total_reads_before, After = total_reads_after) %>%
  pivot_longer(-c(sample, sample_short), names_to = "Timing", values_to = "Reads") %>%
  mutate(Timing = factor(Timing, levels = c("Before", "After")))

# --- Panel 1: Read Counts Before vs After (log10 scale) ---
p_reads <- safe_plot("Read Counts", {
  reads_data %>%
    ggplot(aes(x = sample_short, y = Reads, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>Timing: ", Timing,
                             "<br>Reads: ", scales::comma(Reads)))) +
    geom_col(position = "dodge") +
    # No coord_cartesian ylim: geom_col always draws bars from 0, which floors
    # a log10 axis at raw-value 1 -- cropping further would just clip bars to
    # solid blocks (verified empirically). Left as the plain log10 scale.
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#2ecc71")) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Read Counts Before and After Trimming (log10 scale)",
         x = NULL, y = "Total Reads", fill = NULL)
})

# --- Panel 2: Read Counts Before vs After (linear scale) ---
p_reads_linear <- safe_plot("Read Counts (linear)", {
  reads_data %>%
    ggplot(aes(x = sample_short, y = Reads, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>Timing: ", Timing,
                             "<br>Reads: ", scales::comma(Reads)))) +
    geom_col(position = "dodge") +
    coord_cartesian(ylim = scaled_ylim(reads_data$Reads)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#2ecc71")) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Read Counts Before and After Trimming (linear scale)",
         x = NULL, y = "Total Reads", fill = NULL)
})

# --- Panel 3: Filtering Breakdown ---
p_filter <- safe_plot("Filtering Breakdown", {
  qc %>%
    select(sample, sample_short,
           "Low Quality"  = low_quality_reads,
           "Too Short"    = too_short_reads,
           "Too Long"     = too_long_reads,
           "Too Many N's" = too_many_N_reads) %>%
    pivot_longer(-c(sample, sample_short), names_to = "Reason", values_to = "Reads") %>%
    ggplot(aes(x = sample_short, y = Reads, fill = Reason,
               text = paste0("Sample: ", sample,
                             "<br>Reason: ", Reason,
                             "<br>Reads removed: ", scales::comma(Reads)))) +
    geom_col(position = "stack") +
    # Range is based on the per-sample stack *total* (reads_removed), not the
    # individual reason segments, since the stack total is the bar height a
    # viewer actually sees.
    coord_cartesian(ylim = scaled_ylim(qc$reads_removed)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_brewer(palette = "Reds") +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Reads Removed by Filter Category",
         x = NULL, y = "Reads Removed", fill = "Filter Reason")
})

# --- Panel 3: Q20 / Q30 Rates ---
p_quality <- safe_plot("Q20/Q30 Rates", {
  quality_data <- qc %>%
    select(sample, sample_short,
           "Q20 Before" = q20_before, "Q20 After" = q20_after,
           "Q30 Before" = q30_before, "Q30 After" = q30_after) %>%
    pivot_longer(-c(sample, sample_short), names_to = "Metric", values_to = "Rate") %>%
    mutate(
      Score  = if_else(str_starts(Metric, "Q20"), "Q20", "Q30"),
      Timing = if_else(str_ends(Metric, "Before"), "Before", "After"),
      Timing = factor(Timing, levels = c("Before", "After"))
    )

  # Shared y-limits across both facets (not per-facet free scales), padded 2.5%
  # beyond the observed range and clamped to [0, 100], so Q20 and Q30 stay
  # visually comparable while still keeping both fully visible -- unlike the
  # old hardcoded c(90, 100), which clipped Q30 off-plot on lower-quality runs.
  y_limits <- pmax(0, pmin(100, scaled_ylim(quality_data$Rate)))

  quality_data %>%
    ggplot(aes(x = sample_short, y = Rate, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>", Score, " ", Timing, ": ", round(Rate, 2), "%"))) +
    geom_col(position = "dodge") +
    facet_wrap(~Score) +
    coord_cartesian(ylim = y_limits) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#3498db")) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Q20 / Q30 Rates Before and After Trimming",
         x = NULL, y = "Rate (%)", fill = NULL)
})

# --- Panel 4: GC Content ---
p_gc <- safe_plot("GC Content", {
  gc_data <- qc %>%
    select(sample, sample_short, Before = gc_before, After = gc_after) %>%
    pivot_longer(-c(sample, sample_short), names_to = "Timing", values_to = "GC") %>%
    mutate(Timing = factor(Timing, levels = c("Before", "After")))

  gc_data %>%
    ggplot(aes(x = sample_short, y = GC, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>Timing: ", Timing,
                             "<br>GC: ", round(GC, 2), "%"))) +
    geom_col(position = "dodge") +
    coord_cartesian(ylim = scaled_ylim(gc_data$GC)) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#e67e22")) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "GC Content Before and After Trimming",
         x = NULL, y = "GC Content (%)", fill = NULL)
})

# --- Panel 5: Mean Read Length ---
p_readlen <- safe_plot("Read Length", {
  readlen_data <- qc %>%
    select(sample, sample_short, Before = read_len_before, After = read_len_after) %>%
    pivot_longer(-c(sample, sample_short), names_to = "Timing", values_to = "Length") %>%
    mutate(Timing = factor(Timing, levels = c("Before", "After")))

  readlen_data %>%
    ggplot(aes(x = sample_short, y = Length, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>Timing: ", Timing,
                             "<br>Mean Length: ", round(Length, 1), " bp"))) +
    geom_col(position = "dodge") +
    coord_cartesian(ylim = scaled_ylim(readlen_data$Length)) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#f1c40f")) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Mean Read Length Before and After Trimming",
         x = NULL, y = "Mean Read Length (bp)", fill = NULL)
})

# --- Panel 6: Duplication Rate ---
p_dup <- safe_plot("Duplication Rate", {
  qc %>%
    ggplot(aes(x = sample_short, y = duplication_rate,
               text = paste0("Sample: ", sample,
                             "<br>Duplication Rate: ", round(duplication_rate, 2), "%"))) +
    geom_col(fill = "#9b59b6", alpha = 0.8) +
    coord_cartesian(ylim = scaled_ylim(qc$duplication_rate)) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Duplication Rate per Sample",
         x = NULL, y = "Duplication Rate (%)")
})

# --- Panel 7: Insert Size Peak (paired-end only) ---
p_insert <- safe_plot("Insert Size", {
  if (any(!is.na(qc$insert_size_peak))) {
    qc %>%
      filter(!is.na(insert_size_peak)) %>%
      ggplot(aes(x = sample_short, y = insert_size_peak,
                 text = paste0("Sample: ", sample,
                               "<br>Insert Size Peak: ", insert_size_peak, " bp"))) +
      geom_col(fill = "#1abc9c", alpha = 0.8) +
      coord_cartesian(ylim = scaled_ylim(qc$insert_size_peak)) +
      theme_bw() +
      qc_axis_theme +
      labs(title = "Insert Size Peak",
           x = NULL, y = "Insert Size (bp)")
  } else {
    NULL
  }
})

# --- Panel 8: Adapter Trimming ---
p_adapter <- safe_plot("Adapter Trimming", {
  adapter_data <- qc %>%
    select(sample, sample_short,
           "% Reads"  = pct_reads_adapter_trimmed,
           "% Bases"  = pct_bases_adapter_trimmed) %>%
    pivot_longer(-c(sample, sample_short), names_to = "Metric", values_to = "Percent")

  adapter_data %>%
    ggplot(aes(x = sample_short, y = Percent, fill = Metric,
               text = paste0("Sample: ", sample,
                             "<br>", Metric, " Adapter-Trimmed: ", round(Percent, 2), "%"))) +
    geom_col(position = "dodge") +
    coord_cartesian(ylim = scaled_ylim(adapter_data$Percent)) +
    scale_fill_manual(values = c("% Reads" = "#34495e", "% Bases" = "#16a085")) +
    theme_bw() +
    qc_axis_theme +
    labs(title = "Adapter Trimming per Sample",
         x = NULL, y = "Adapter-Trimmed (%)", fill = NULL)
})

# --- Panel 9: Mean Quality per Sequencing Cycle ---
p_cycles <- safe_plot("Quality per Cycle", {
  if (nrow(cycles) == 0) return(NULL)

  cycle_summary <- cycles %>%
    mutate(timing = factor(timing, levels = c("Before", "After")))

  # Population mean line across all samples per cycle/read/timing
  cycle_mean <- cycle_summary %>%
    group_by(read, timing, cycle) %>%
    summarise(mean_quality = mean(mean_quality, na.rm = TRUE), .groups = "drop")

  ggplot() +
    geom_line(data = cycle_summary,
              aes(x = cycle, y = mean_quality, group = sample,
                  text = paste0("Sample: ", sample,
                                "<br>Cycle: ", cycle,
                                "<br>Mean Quality: ", round(mean_quality, 2))),
              alpha = 0.35, linewidth = 0.4, color = "#3498db") +
    geom_line(data = cycle_mean,
              aes(x = cycle, y = mean_quality),
              color = "#e74c3c", linewidth = 0.9, linetype = "solid") +
    geom_hline(yintercept = 30, linetype = "dashed", color = "black", alpha = 0.5) +
    facet_grid(read ~ timing) +
    coord_cartesian(ylim = scaled_ylim(cycle_summary$mean_quality)) +
    theme_bw() +
    labs(title = "Mean Base Quality per Sequencing Cycle",
         subtitle = "Blue lines = individual samples. Red line = cohort mean. Dashed = Q30.",
         x = "Cycle", y = "Mean Quality Score")
})

# --- Panel tiers: single source of truth for both static and interactive layouts ---
# Each tier is a row; ncol = 1 spans the tier's plot(s) across the full figure
# width (e.g. panel I's Q20/Q30 facets, panel J's per-cycle lines), ncol = 2
# lays same-tier plots side by side in a grid. `height` (inches) lets some
# rows run shorter than the busier ones -- A-H all share the same row height
# so every simple per-sample bar panel reads at the same scale.
tiers <- list(
  list(plots = list(p_reads, p_reads_linear), ncol = 2, height = 4),  # A, B
  list(plots = list(p_filter, p_readlen),      ncol = 2, height = 4),  # C, D
  list(plots = list(p_gc, p_dup),              ncol = 2, height = 4),  # E, F
  list(plots = list(p_insert, p_adapter),      ncol = 2, height = 4),  # G, H
  list(plots = list(p_quality),                ncol = 1, height = 4),  # I -- full width, just above J
  list(plots = list(p_cycles),                 ncol = 1, height = 12)   # J -- full width
)
tiers <- lapply(tiers, function(t) { t$plots <- Filter(Negate(is.null), t$plots); t })
tiers <- Filter(function(t) length(t$plots) > 0, tiers)
tier_heights <- sapply(tiers, function(t) ceiling(length(t$plots) / t$ncol) * t$height)

# --- Static Output (patchwork combined panel) ---
if (length(tiers) > 0) {
  tryCatch({
    annotation <- plot_annotation(
      title      = paste0("FastP QC Panel: ", PREFIX),
      subtitle   = paste0(nrow(qc), " samples"),
      tag_levels = "A",
      caption    = {
        lines <- c(
          "(A) Total read count per sample before (gray) and after (green) FastP adapter and quality trimming, shown as grouped bars on a log10-scaled y-axis.",
          "(B) The same read counts as (A) on a linear y-axis, zoomed to the observed data range (± 2.5%) for direct sample-to-sample comparison.",
          "(C) Read counts discarded during trimming, stacked by filter category: Low Quality, Too Short, Too Long, Too Many N's.",
          "(D) Mean read length (bp) per sample before and after trimming, shown as grouped bars.",
          "(E) GC content (%) per sample before and after trimming, shown as grouped bars.",
          "(F) Estimated duplicate rate (%) per sample as reported by FastP, reflecting PCR and optical duplicates. FastP reports a single rate estimated from the raw input reads -- this pipeline does not run FastP's --dedup flag, so no reads are actually removed for duplication and there is no before/after split to show.",
          "(G) Peak insert size (bp) per sample estimated from paired-end read overlap, where available. As with duplication rate, FastP reports only a single estimate here, not a before/after comparison.",
          "(H) Percentage of reads and bases removed by adapter trimming per sample.",
          "(I) Percentage of bases meeting Q20 (≥99% accuracy) and Q30 (≥99.9% accuracy) quality thresholds, before and after trimming. Y-axis is shared across both facets and scaled to the observed data range (± 2.5%, clamped to 0-100%) so Q20 and Q30 stay directly comparable.",
          "(J) Mean base quality score at each sequencing cycle for R1 and R2, before and after trimming. Blue lines = individual samples; red line = cohort mean; dashed line = Q30 threshold."
        )
        paste(sapply(lines, stringr::str_wrap, width = 160), collapse = "\n")
      },
      theme = theme(
        plot.title   = element_text(size = 16, face = "bold"),
        plot.caption = element_text(size = 8,  hjust = 0, lineheight = 1.4)
      )
    )

    total_height <- sum(tier_heights)
    tier_plots   <- lapply(tiers, function(t) wrap_plots(t$plots, ncol = t$ncol))
    # Reduce(`/`, ...) alone gives every stacked tier equal vertical space --
    # plot_layout(heights=) is what actually makes each tier's row occupy
    # space proportional to its `height` field.
    combined     <- Reduce(`/`, tier_plots) +
      plot_layout(heights = tier_heights) +
      annotation

    ggsave(paste0(PREFIX, "_QC_fastp_panel.jpg"), plot = combined,
           width = plot_width, height = total_height, dpi = 300)
  }, error = function(e) {
    message(sprintf("[WARN] Static fastp panel save failed: %s", conditionMessage(e)))
  })
}

# --- Interactive Output (plotly subplots) ---
tryCatch({
  tier_interactive <- Filter(Negate(is.null), lapply(tiers, function(t) {
    ps <- lapply(t$plots, function(p) ggplotly(p, tooltip = "text") %>% partial_bundle())
    if (length(ps) == 1) return(ps[[1]])
    subplot(ps, nrows = ceiling(length(ps) / t$ncol),
            shareX = FALSE, shareY = FALSE, titleX = TRUE, titleY = TRUE)
  }))

  if (length(tier_interactive) > 0) {
    # tiers (and therefore tier_heights) were already filtered to non-empty
    # entries before tier_interactive was built, so the two stay aligned --
    # heights must sum to 1 for plotly::subplot.
    combined_interactive <- subplot(tier_interactive,
                                    nrows   = length(tier_interactive),
                                    heights = tier_heights / sum(tier_heights),
                                    shareX  = FALSE,
                                    shareY  = FALSE,
                                    titleX  = TRUE,
                                    titleY  = TRUE) %>%
      layout(title = paste0("FastP QC Panel: ", PREFIX))
    saveWidget(combined_interactive, paste0(PREFIX, "_QC_fastp_panel.html"), selfcontained = TRUE)
  }
}, error = function(e) {
  message(sprintf("[WARN] Interactive fastp panel failed: %s", conditionMessage(e)))
})
