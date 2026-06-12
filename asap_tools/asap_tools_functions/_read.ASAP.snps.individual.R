read.ASAP.snps.individual <- function(XML) {
  library(xml2)
  library(tidyverse)

  Out      <- data.frame()
  xml_data <- read_xml(XML, options = "HUGE")
  Run_Info <- data.frame(run = "Individual_XML_processing")

  Sample_Node <- xml_data

  Sample_Info <- data.frame(
    name             = xml_attr(Sample_Node, "name"),
    Total_Reads      = xml_attr(Sample_Node, "total_reads"),
    Trimmed_Reads    = xml_attr(Sample_Node, "trimmed_reads"),
    Mapped_Reads     = xml_attr(Sample_Node, "mapped_reads"),
    unassigned_reads = xml_attr(Sample_Node, "unassigned_reads"),
    unmapped_reads   = xml_attr(Sample_Node, "unmapped_reads")
  )

  Assays <- xml_children(Sample_Node)

  for (i in seq_along(Assays)) {
    Assay_Node <- Assays[[i]]

    Assay_Info <- data.frame(
      assay_function = xml_attr(Assay_Node, "function"),
      assay_gene     = xml_attr(Assay_Node, "gene"),
      assay_name     = xml_attr(Assay_Node, "name"),
      assay_type     = xml_attr(Assay_Node, "type")
    )

    Amplicons <- xml_children(Assay_Node)

    for (j in seq_along(Amplicons)) {
      Amplicon_Node <- Amplicons[[j]]
      Amplicon_Info <- data.frame(amplicon_number = j)

      All_Children <- xml_children(Amplicon_Node)
      SNP_Nodes    <- All_Children[xml_name(All_Children) == "snp"]

      if (length(SNP_Nodes) > 0) {
        for (k in seq_along(SNP_Nodes)) {
          SNP_Node <- SNP_Nodes[[k]]

          Dist_Node <- xml_child(SNP_Node, "base_distribution")
          if (!inherits(Dist_Node, "xml_missing")) {
            attrs    <- xml_attrs(Dist_Node)
            SNP_Dist <- paste(names(attrs), attrs, sep = "=", collapse = ", ")
          } else {
            SNP_Dist <- "No reads matching SNP"
          }

          Call_Node <- xml_child(SNP_Node, "snp_call")
          if (!inherits(Call_Node, "xml_missing")) {
            snp_depth      <- xml_attr(Call_Node, "count")
            snp_proportion <- xml_attr(Call_Node, "percent")
            snp_call       <- xml_text(Call_Node)
          } else {
            snp_depth <- "0"; snp_proportion <- "0"; snp_call <- "N/A"
          }

          CodonMerge_Nodes <- xml_find_all(SNP_Node, "codon_merge")
          if (length(CodonMerge_Nodes) > 0) {
            cm_distributions <- character(length(CodonMerge_Nodes))
            for (m in seq_along(CodonMerge_Nodes)) {
              dist_node <- xml_find_first(CodonMerge_Nodes[[m]], "codon_distribution")
              if (!inherits(dist_node, "xml_missing")) {
                attrs <- xml_attrs(dist_node)
                cm_distributions[m] <- paste(names(attrs), attrs, sep = "=", collapse = ",")
              } else {
                cm_distributions[m] <- ""
              }
            }

            codon_merge_name         <- paste(xml_attr(CodonMerge_Nodes, "name"),      collapse = ";")
            codon_merge_region       <- paste(xml_attr(CodonMerge_Nodes, "region"),    collapse = ";")
            codon_merge_position     <- paste(xml_attr(CodonMerge_Nodes, "position"),  collapse = ";")
            codon_merge_codon_depth  <- paste(xml_attr(CodonMerge_Nodes, "codon_depth"), collapse = ";")
            codon_merge_reference    <- paste(xml_attr(CodonMerge_Nodes, "reference"), collapse = ";")
            codon_merge_distribution <- paste(cm_distributions, collapse = ";")
          } else {
            codon_merge_name         <- NA_character_
            codon_merge_region       <- NA_character_
            codon_merge_position     <- NA_character_
            codon_merge_codon_depth  <- NA_character_
            codon_merge_reference    <- NA_character_
            codon_merge_distribution <- NA_character_
          }

          LinkedSNP_Nodes <- xml_find_all(SNP_Node, "linked_snps/linked_snp")
          if (length(LinkedSNP_Nodes) > 0) {
            linked_snp_targets       <- paste(xml_attr(LinkedSNP_Nodes, "target_variant"),      collapse = ";")
            linked_snp_linkage_pcts  <- paste(xml_attr(LinkedSNP_Nodes, "linkage_pct"),         collapse = ";")
            linked_snp_co_counts     <- paste(xml_attr(LinkedSNP_Nodes, "co_occurring_count"),   collapse = ";")
            linked_snp_shared_depths <- paste(xml_attr(LinkedSNP_Nodes, "shared_read_depth"),    collapse = ";")
          } else {
            linked_snp_targets       <- NA_character_
            linked_snp_linkage_pcts  <- NA_character_
            linked_snp_co_counts     <- NA_character_
            linked_snp_shared_depths <- NA_character_
          }

          call_qual_node <- xml_find_first(SNP_Node,
            paste0("base_quality/qual[@base='", snp_call, "']"))
          if (!inherits(call_qual_node, "xml_missing")) {
            snp_call_qual_mean   <- xml_attr(call_qual_node, "mean")
            snp_call_qual_median <- xml_attr(call_qual_node, "median")
            snp_call_qual_min    <- xml_attr(call_qual_node, "min")
            snp_call_qual_max    <- xml_attr(call_qual_node, "max")
          } else {
            snp_call_qual_mean <- snp_call_qual_median <- snp_call_qual_min <- snp_call_qual_max <- NA_character_
          }

          strand_node <- xml_find_first(SNP_Node,
            paste0("base_strand_distribution/strand[@base='", snp_call, "']"))
          if (!inherits(strand_node, "xml_missing")) {
            snp_call_R1 <- xml_attr(strand_node, "R1")
            snp_call_R2 <- xml_attr(strand_node, "R2")
          } else {
            snp_call_R1 <- snp_call_R2 <- NA_character_
          }

          SNP_Info <- data.frame(
            location_depth   = xml_attr(SNP_Node, "depth"),
            snp_name         = xml_attr(SNP_Node, "name"),
            snp_position     = xml_attr(SNP_Node, "position"),
            snp_reference    = xml_attr(SNP_Node, "reference"),
            snp_depth        = snp_depth,
            snp_proportion   = snp_proportion,
            snp_call         = snp_call,
            snp_distribution = SNP_Dist,
            codon_merge_name         = codon_merge_name,
            codon_merge_region       = codon_merge_region,
            codon_merge_position     = codon_merge_position,
            codon_merge_codon_depth  = codon_merge_codon_depth,
            codon_merge_reference    = codon_merge_reference,
            codon_merge_distribution = codon_merge_distribution,
            linked_snp_targets       = linked_snp_targets,
            linked_snp_linkage_pcts  = linked_snp_linkage_pcts,
            linked_snp_co_counts     = linked_snp_co_counts,
            linked_snp_shared_depths = linked_snp_shared_depths,
            snp_call_qual_mean       = snp_call_qual_mean,
            snp_call_qual_median     = snp_call_qual_median,
            snp_call_qual_min        = snp_call_qual_min,
            snp_call_qual_max        = snp_call_qual_max,
            snp_call_R1              = snp_call_R1,
            snp_call_R2              = snp_call_R2
          )

          Temp <- cbind(Run_Info, Sample_Info, Assay_Info, Amplicon_Info, SNP_Info)
          Out  <- rbind(Out, Temp)
        }
      }
    }
  }

  print(paste("Processing complete for:", xml_attr(Sample_Node, "name")))

  if (nrow(Out) > 0) {
    row.names(Out) <- 1:nrow(Out)
    num_cols <- c("Total_Reads", "Trimmed_Reads", "Mapped_Reads", "unassigned_reads",
                  "unmapped_reads", "location_depth", "snp_position", "snp_depth", "snp_proportion",
                  "snp_call_qual_mean", "snp_call_qual_median", "snp_call_qual_min", "snp_call_qual_max",
                  "snp_call_R1", "snp_call_R2")
    num_cols <- intersect(num_cols, names(Out))
    Out[num_cols] <- lapply(Out[num_cols], function(x) as.numeric(as.character(x)))
  }

  return(Out)
}
