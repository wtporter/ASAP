#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)
library(plotly)
library(htmlwidgets)
library(patchwork)

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

  tibble(
    sample              = sample_id,
    total_reads_before  = bf$total_reads,
    total_reads_after   = af$total_reads,
    reads_removed       = bf$total_reads - af$total_reads,
    q20_before          = bf$q20_rate * 100,
    q20_after           = af$q20_rate  * 100,
    q30_before          = bf$q30_rate  * 100,
    q30_after           = af$q30_rate  * 100,
    gc_before           = bf$gc_content * 100,
    gc_after            = af$gc_content * 100,
    read_len_before     = bf$read1_mean_length,
    read_len_after      = af$read1_mean_length,
    low_quality_reads   = fr$low_quality_reads %||% 0L,
    too_short_reads     = fr$too_short_reads   %||% 0L,
    too_long_reads      = fr$too_long_reads    %||% 0L,
    too_many_N_reads    = fr$too_many_N_reads  %||% 0L,
    duplication_rate    = dup * 100,
    insert_size_peak    = if (!is.null(ins)) ins else NA_real_
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

cycles <- map_dfr(json_files, safely(parse_fastp_cycles, otherwise = NULL)) %>%
  { bind_rows(.$result) }

# --- Panel 1: Read Counts Before vs After ---
p_reads <- safe_plot("Read Counts", {
  qc %>%
    select(sample, Before = total_reads_before, After = total_reads_after) %>%
    pivot_longer(-sample, names_to = "Timing", values_to = "Reads") %>%
    mutate(Timing = factor(Timing, levels = c("Before", "After"))) %>%
    ggplot(aes(x = sample, y = Reads, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>Timing: ", Timing,
                             "<br>Reads: ", scales::comma(Reads)))) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#2ecc71")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Read Counts Before and After Trimming",
         x = NULL, y = "Total Reads", fill = NULL)
})

# --- Panel 2: Filtering Breakdown ---
p_filter <- safe_plot("Filtering Breakdown", {
  qc %>%
    select(sample,
           "Low Quality"  = low_quality_reads,
           "Too Short"    = too_short_reads,
           "Too Long"     = too_long_reads,
           "Too Many N's" = too_many_N_reads) %>%
    pivot_longer(-sample, names_to = "Reason", values_to = "Reads") %>%
    ggplot(aes(x = sample, y = Reads, fill = Reason,
               text = paste0("Sample: ", sample,
                             "<br>Reason: ", Reason,
                             "<br>Reads removed: ", scales::comma(Reads)))) +
    geom_col(position = "stack") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_brewer(palette = "Reds") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Reads Removed by Filter Category",
         x = NULL, y = "Reads Removed", fill = "Filter Reason")
})

# --- Panel 3: Q20 / Q30 Rates ---
p_quality <- safe_plot("Q20/Q30 Rates", {
  qc %>%
    select(sample,
           "Q20 Before" = q20_before, "Q20 After" = q20_after,
           "Q30 Before" = q30_before, "Q30 After" = q30_after) %>%
    pivot_longer(-sample, names_to = "Metric", values_to = "Rate") %>%
    mutate(
      Score  = if_else(str_starts(Metric, "Q20"), "Q20", "Q30"),
      Timing = if_else(str_ends(Metric, "Before"), "Before", "After"),
      Timing = factor(Timing, levels = c("Before", "After"))
    ) %>%
    ggplot(aes(x = sample, y = Rate, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>", Score, " ", Timing, ": ", round(Rate, 2), "%"))) +
    geom_col(position = "dodge") +
    facet_wrap(~Score) +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#3498db")) +
    coord_cartesian(ylim = c(90, 100)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Q20 / Q30 Rates Before and After Trimming",
         x = NULL, y = "Rate (%)", fill = NULL)
})

# --- Panel 4: GC Content ---
p_gc <- safe_plot("GC Content", {
  qc %>%
    select(sample, Before = gc_before, After = gc_after) %>%
    pivot_longer(-sample, names_to = "Timing", values_to = "GC") %>%
    mutate(Timing = factor(Timing, levels = c("Before", "After"))) %>%
    ggplot(aes(x = sample, y = GC, fill = Timing,
               text = paste0("Sample: ", sample,
                             "<br>Timing: ", Timing,
                             "<br>GC: ", round(GC, 2), "%"))) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c(Before = "#95a5a6", After = "#e67e22")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "GC Content Before and After Trimming",
         x = NULL, y = "GC Content (%)", fill = NULL)
})

# --- Panel 5: Duplication Rate ---
p_dup <- safe_plot("Duplication Rate", {
  qc %>%
    ggplot(aes(x = sample, y = duplication_rate,
               text = paste0("Sample: ", sample,
                             "<br>Duplication Rate: ", round(duplication_rate, 2), "%"))) +
    geom_col(fill = "#9b59b6", alpha = 0.8) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Duplication Rate per Sample",
         x = NULL, y = "Duplication Rate (%)")
})

# --- Panel 6: Insert Size Peak (paired-end only) ---
p_insert <- safe_plot("Insert Size", {
  if (any(!is.na(qc$insert_size_peak))) {
    qc %>%
      filter(!is.na(insert_size_peak)) %>%
      ggplot(aes(x = sample, y = insert_size_peak,
                 text = paste0("Sample: ", sample,
                               "<br>Insert Size Peak: ", insert_size_peak, " bp"))) +
      geom_col(fill = "#1abc9c", alpha = 0.8) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Insert Size Peak",
           x = NULL, y = "Insert Size (bp)")
  } else {
    NULL
  }
})

# --- Panel 7: Mean Quality per Sequencing Cycle ---
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
    theme_bw() +
    labs(title = "Mean Base Quality per Sequencing Cycle",
         subtitle = "Blue lines = individual samples. Red line = cohort mean. Dashed = Q30.",
         x = "Cycle", y = "Mean Quality Score")
})

# --- Static Output (patchwork combined panel) ---
# p_cycles spans full width at the bottom; all other plots are in a 2-column grid above it.
grid_plots      <- Filter(Negate(is.null), list(p_reads, p_filter, p_quality, p_gc, p_dup, p_insert))
fullwidth_plots <- Filter(Negate(is.null), list(p_cycles))

if (length(grid_plots) > 0 || length(fullwidth_plots) > 0) {
  tryCatch({
    annotation <- plot_annotation(
      title      = paste0("FastP QC Panel: ", PREFIX),
      subtitle   = paste0(nrow(qc), " samples"),
      tag_levels = "A",
      caption    = {
        lines <- c(
          "(A) Total read count per sample before (gray) and after (green) FastP adapter and quality trimming, shown as grouped bars.",
          "(B) Read counts discarded during trimming, stacked by filter category: Low Quality, Too Short, Too Long, Too Many N's.",
          "(C) Percentage of bases meeting Q20 (≥99% accuracy) and Q30 (≥99.9% accuracy) quality thresholds, before and after trimming. Y-axis starts at 90%.",
          "(D) GC content (%) per sample before and after trimming, shown as grouped bars.",
          "(E) Estimated duplicate rate (%) per sample as reported by FastP, reflecting PCR and optical duplicates.",
          "(F) Peak insert size (bp) per sample estimated from paired-end read overlap, where available.",
          "(G) Mean base quality score at each sequencing cycle for R1 and R2, before and after trimming. Blue lines = individual samples; red line = cohort mean; dashed line = Q30 threshold."
        )
        paste(sapply(lines, stringr::str_wrap, width = 160), collapse = "\n")
      },
      theme = theme(
        plot.title   = element_text(size = 16, face = "bold"),
        plot.caption = element_text(size = 8,  hjust = 0, lineheight = 1.4)
      )
    )

    n_grid_rows     <- ceiling(length(grid_plots) / 2)
    n_fullwidth_rows <- length(fullwidth_plots)
    total_height    <- (n_grid_rows + n_fullwidth_rows) * 6

    if (length(grid_plots) > 0 && length(fullwidth_plots) > 0) {
      combined <- wrap_plots(grid_plots, ncol = 2) /
                  wrap_plots(fullwidth_plots, ncol = 1)
    } else if (length(grid_plots) > 0) {
      combined <- wrap_plots(grid_plots, ncol = 2)
    } else {
      combined <- wrap_plots(fullwidth_plots, ncol = 1)
    }

    combined <- combined + annotation
    ggsave(paste0(PREFIX, "_QC_fastp_panel.jpg"), plot = combined,
           width = 18, height = total_height, dpi = 300)
  }, error = function(e) {
    message(sprintf("[WARN] Static fastp panel save failed: %s", conditionMessage(e)))
  })
}

# --- Interactive Output (plotly subplots) ---
tryCatch({
  interactive_panels <- Filter(Negate(is.null), list(
    if (!is.null(p_reads))   ggplotly(p_reads,   tooltip = "text") %>% partial_bundle(),
    if (!is.null(p_filter))  ggplotly(p_filter,  tooltip = "text") %>% partial_bundle(),
    if (!is.null(p_quality)) ggplotly(p_quality, tooltip = "text") %>% partial_bundle(),
    if (!is.null(p_gc))      ggplotly(p_gc,      tooltip = "text") %>% partial_bundle(),
    if (!is.null(p_dup))     ggplotly(p_dup,     tooltip = "text") %>% partial_bundle(),
    if (!is.null(p_insert))  ggplotly(p_insert,  tooltip = "text") %>% partial_bundle(),
    if (!is.null(p_cycles))  ggplotly(p_cycles,  tooltip = "text") %>% partial_bundle()
  ))

  if (length(interactive_panels) > 0) {
    n_cols <- 2
    n_rows <- ceiling(length(interactive_panels) / n_cols)
    combined_interactive <- subplot(interactive_panels,
                                    nrows   = n_rows,
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
