#!/usr/bin/env Rscript

library(tidyverse)
library(plotly)
library(htmlwidgets)
library(patchwork)
library(genbankr)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: process_asaptools_snp_figures.R <rdata> <prefix> <poi_csv> [<snp_threshold>] [<snp_depth>] [<gb1> ...]")
}

RDATA_INPUT   <- args[1]
PREFIX        <- args[2]
POI_CSV       <- args[3]
SNP_THRESHOLD <- if (length(args) >= 4 && !args[4] %in% c("NULL", "NA", "")) as.numeric(args[4]) else 0.03
MIN_DEPTH     <- if (length(args) >= 5 && !args[5] %in% c("NULL", "NA", "")) as.numeric(args[5]) else 100
GB_FILES      <- if (length(args) >= 6) args[6:length(args)] else character(0)

safe_plot <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[WARN] %s failed: %s", label, conditionMessage(e)))
  })
}

wrap_cap <- function(txt, w = 110) {
  paste(sapply(strsplit(txt, "\n")[[1]], stringr::str_wrap, width = w), collapse = "\n")
}

load(RDATA_INPUT)

array_info <- final_array
SNPS       <- final_snps

# Join per-position depth so plots can depth-filter; keep snp_distribution raw
# so each plot block can expand it inline (parse_snp_distribution destroys the column)
SNPS <- full_join(
  select(final_array, run, name, assay_name, position, depth),
  SNPS,
  by = c("run", "name", "assay_name", "position" = "snp_position")
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
  SNPS_plot <- SNPS %>% mutate(space_count = str_count(snp_distribution, " "))
  max_spaces <- max(SNPS_plot$space_count, na.rm = TRUE)

  SNPS_plot <- SNPS_plot %>%
    relocate(snp_distribution, .after = last_col()) %>%
    separate(snp_distribution, into = paste0("Dist", 1:(1 + max_spaces)), sep = ", ", fill = "right") %>%
    pivot_longer(starts_with("Dist"), names_to = "Temp", values_to = "Dist") %>%
    select(-Temp) %>%
    filter(!is.na(Dist)) %>%
    separate(Dist, into = c("Call", "n"), sep = "=") %>%
    mutate(snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))) %>%
    mutate(SNP = paste0(snp_reference, position, Call)) %>%
    filter(snp_reference != Call) %>%
    filter(!is.na(snp_proportion))

  SNPS_plot <- SNPS_plot %>% filter(depth >= MIN_DEPTH)
  SNPS_plot$snp_proportion[is.na(SNPS_plot$snp_proportion)] <- 0

  SNP_Plot_Data <- SNPS_plot %>%
    group_by(name, assay_name, position) %>%
    summarise(Max_SNP_proportion = max(snp_proportion), .groups = "drop")

  p_SNP <- SNP_Plot_Data %>%
    ggplot(aes(x = position, y = Max_SNP_proportion, col = name,
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

  ggsave(paste0(PREFIX, "_SNP_position_prevalence.jpg"), plot = p_SNP, width = 12, height = 8, dpi = 300)
  interactive_plot_SNP <- ggplotly(p_SNP, tooltip = "text") %>% partial_bundle()
  saveWidget(interactive_plot_SNP, paste0(PREFIX, "_SNP_position_prevalence.html"), selfcontained = TRUE)

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

  ggsave(paste0(PREFIX, "_SNP_proportion_density.jpg"), plot = p_snp_dist, width = 12, height = 8, dpi = 300)
  interactive_snp_dist <- ggplotly(p_snp_dist, tooltip = "text") %>% partial_bundle()
  saveWidget(interactive_snp_dist, paste0(PREFIX, "_SNP_proportion_density.html"), selfcontained = TRUE)
})

# --- Plot: Strand Bias ---
safe_plot("Strand Bias", {
  if (all(c("snp_call_R1", "snp_call_R2") %in% names(SNPS))) {
    strand_data <- SNPS %>%
      filter(!is.na(snp_call_R1), !is.na(snp_call_R2)) %>%
      mutate(
        total_strand   = as.numeric(snp_call_R1) + as.numeric(snp_call_R2),
        strand_ratio   = as.numeric(snp_call_R1) / total_strand,
        snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))
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

    ggsave(paste0(PREFIX, "_SNP_strand_bias.jpg"), plot = p_strand, width = 12, height = 8, dpi = 300)
    interactive_strand <- ggplotly(p_strand, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_strand, paste0(PREFIX, "_SNP_strand_bias.html"), selfcontained = TRUE)
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
      group_by(name, assay_name) %>%
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
                aes(x = Type, y = Mean_Quality, group = name),
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

    ggsave(paste0(PREFIX, "_SNP_base_quality.jpg"), plot = p_qual, width = 12, height = 8, dpi = 300)
    interactive_qual <- ggplotly(p_qual, tooltip = "text") %>% partial_bundle()
    saveWidget(interactive_qual, paste0(PREFIX, "_SNP_base_quality.html"), selfcontained = TRUE)
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

  snp_local <- SNPS %>%
    mutate(space_count = str_count(snp_distribution, " "))
  max_sp <- max(snp_local$space_count, na.rm = TRUE)

  snp_local <- snp_local %>%
    relocate(snp_distribution, .after = last_col()) %>%
    separate(snp_distribution, into = paste0("Dist", 1:(1 + max_sp)), sep = ", ", fill = "right") %>%
    pivot_longer(starts_with("Dist"), names_to = "Temp", values_to = "Dist") %>%
    select(-Temp) %>%
    filter(!is.na(Dist)) %>%
    separate(Dist, into = c("Call", "n"), sep = "=") %>%
    mutate(snp_proportion = 100 * (as.numeric(n) / as.numeric(location_depth))) %>%
    filter(snp_reference != Call, !is.na(snp_proportion), depth >= MIN_DEPTH)

  snp_pos_data <- snp_local %>%
    group_by(name, assay_name, position) %>%
    summarise(Max_SNP_proportion = max(snp_proportion), .groups = "drop")

  gene_pal <- c("#4E79A7","#59A14F","#F28E2B","#E15759","#B07AA1","#9C755F",
                "#FF9DA7","#76B7B2","#EDC948","#BAB0AC","#D37295","#FABFD2")

  for (gb_file in GB_FILES) {
    tryCatch({
      ref    <- suppressWarnings(genbankr::readGenBank(gb_file))
      acc_id <- ref@accession
      fb     <- sub("\\.[^.]+$", "", basename(gb_file))

      gene_df <- left_join(
        data.frame(ref@genes),
        data.frame(ref@cds) %>% select(locus_tag, product)
      ) %>%
        mutate(
          gene_label  = ifelse(is.na(gene), as.character(locus_tag), as.character(gene)),
          label_angle = ifelse((end - start) < nchar(as.character(ref@sequence)) / 20, 90, 0),
          gene_color  = gene_pal[(seq_len(n()) - 1) %% length(gene_pal) + 1]
        ) %>%
        filter(!is.na(start), !is.na(end))

      genome_len  <- nchar(as.character(ref@sequence))
      gene_bounds <- sort(unique(c(gene_df$start, gene_df$end)))
      x_lim       <- c(0, genome_len)

      snp_match <- snp_pos_data %>%
        filter(grepl(acc_id, assay_name, fixed = TRUE) |
               grepl(fb,     assay_name, fixed = TRUE))

      cov_match <- array_info %>%
        filter(grepl(acc_id, assay_name, fixed = TRUE) |
               grepl(fb,     assay_name, fixed = TRUE)) %>%
        select(name, position, depth)

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
          group_by(name, run_id) %>%
          summarise(
            cov_xmid  = (min(position) + max(position)) / 2,
            cov_width = max(position) - min(position) + 1L,
            .groups = "drop"
          )
      } else {
        tibble(name = character(), cov_xmid = numeric(), cov_width = integer())
      }

      sample_order <- sort(unique(snp_match$name))
      n_samples    <- length(sample_order)
      snp_match    <- mutate(snp_match, name = factor(name, levels = sample_order))
      cov_tiles    <- mutate(cov_tiles, name = factor(name, levels = sample_order))
      snp_width    <- max(1, genome_len / 400)

      p_genes <- gene_df %>%
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = gene_label),
                  color = "white", linewidth = 0.4, alpha = 0.9) +
        geom_text(aes(x = (start + end) / 2, y = 0.5,
                      label = gene_label, angle = label_angle),
                  size = 2.8, color = "white", fontface = "bold") +
        scale_fill_manual(values = setNames(gene_df$gene_color, gene_df$gene_label)) +
        scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = scales::comma) +
        theme_void() +
        theme(legend.position = "none", plot.margin = margin(b = 2)) +
        labs(title = paste0("Genome Track: ", acc_id))

      p_density <- snp_match %>%
        ggplot(aes(x = as.numeric(position))) +
        geom_density(fill = "grey80", alpha = 0.5, color = "black",
                     bw = max(50, genome_len / 100)) +
        geom_vline(xintercept = gene_bounds, color = "grey60", alpha = 0.35, linewidth = 0.3) +
        scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = scales::comma) +
        theme_minimal() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), panel.grid = element_blank(),
              plot.margin = margin(t = 0, b = 0)) +
        labs(y = "SNP Density")

      p_heat <- ggplot() +
        geom_tile(data = cov_tiles,
                  aes(x = cov_xmid, y = name, width = cov_width, height = 0.85),
                  fill = "grey50", alpha = 0.18) +
        geom_vline(xintercept = gene_bounds, color = "grey60", alpha = 0.35, linewidth = 0.3) +
        geom_tile(data = snp_match,
                  aes(x = as.numeric(position), y = name,
                      text = paste0("Sample: ", name,
                                    "<br>Position: ", position,
                                    "<br>Max SNP %: ", round(Max_SNP_proportion, 1), "%")),
                  width = snp_width, height = 0.85, fill = "red") +
        scale_x_continuous(limits = x_lim, expand = c(0, 0), labels = scales::comma) +
        theme_minimal() +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1),
              axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
              panel.grid.major.x = element_line(color = "grey95", linewidth = 0.1),
              panel.grid.minor   = element_blank(),
              panel.grid.major.y = element_blank(),
              plot.margin = margin(t = 0)) +
        labs(x = "Reference Position (BP)",
             y = paste0("Samples (n=", n_samples, ")"),
             caption = wrap_cap(paste0(
               "(A) CDS gene annotations parsed from ", basename(gb_file), ". ",
               "(B) Kernel density estimate of SNP positions across the genome. ",
               "(C) Per-sample SNP locations (red) overlaid on covered regions (grey). ",
               "Grey shading = positions with depth > 0. Only positions meeting minimum depth are eligible for SNP calling."
             ))) +
        theme(plot.caption = element_text(hjust = 0, size = 8, lineheight = 1.3))

      combined_track <- p_genes / p_density / p_heat +
        plot_layout(heights = c(1, 2, max(3, n_samples * 0.3))) +
        plot_annotation(
          tag_levels = "A", tag_suffix = ".",
          theme = theme(
            plot.tag    = element_text(face = "bold", size = 14),
            plot.margin = margin(l = 5, r = 5, t = 5, b = 5)
          )
        )

      safe_name    <- gsub("[^A-Za-z0-9]", "_", acc_id)
      track_height <- max(8, 3.5 + n_samples * 0.3)

      ggsave(paste0(PREFIX, "_SNP_", safe_name, "_genome_track.jpg"),
             plot = combined_track, width = 18, height = track_height, dpi = 300)

      interactive_track <- ggplotly(p_heat, tooltip = "text") %>% partial_bundle()
      saveWidget(interactive_track,
                 paste0(PREFIX, "_SNP_", safe_name, "_genome_track.html"),
                 selfcontained = TRUE)

    }, error = function(e) {
      message(sprintf("[WARN] Genome Track for %s failed: %s",
                      basename(gb_file), conditionMessage(e)))
    })
  }
})
