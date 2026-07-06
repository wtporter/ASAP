library(tidyverse)
library(genbankr)
library(Biostrings)
library(foreach)
library(doParallel)
library(parallelly)

# ---------------------------------------------------------------------------
# Pure translation helpers — take scalar gene metadata + SNP vectors, return
# a single-row data.frame. No file I/O; directly unit-testable.
# ---------------------------------------------------------------------------

.translate_snp <- function(position_vec, mutation_vec, genome_snp,
                            gene_seq, gene_start, gene_end, gene_strand,
                            gene_name, gene_product) {
  n_comp       <- length(position_vec)
  gene_seq_dna <- Biostrings::DNAString(gene_seq)
  snp_in_gene  <- (position_vec - gene_start) + 1

  gene_len <- nchar(gene_seq)
  if (any(snp_in_gene < 1L) || any(snp_in_gene > gene_len)) {
    return(data.frame(
      SNP                   = as.character(genome_snp),
      snp_position_genome   = paste(position_vec, collapse = "|"),
      snp_position_gene     = paste(snp_in_gene,  collapse = "|"),
      Theoretical_Reference = NA_character_,
      Gene                  = as.character(gene_name),
      Product               = as.character(gene_product),
      AA                    = "Non-coding",
      SNP_Gene              = NA_character_,
      stringsAsFactors      = FALSE
    ))
  }

  ref_bases <- character(n_comp)
  obs_seq   <- gene_seq_dna
  for (i in seq_len(n_comp)) {
    ref_bases[i]          <- as.character(gene_seq_dna[snp_in_gene[i]])
    obs_seq[snp_in_gene[i]] <- mutation_vec[i]
  }
  mut_out <- mutation_vec

  if (gene_strand == "-") {
    gene_seq_dna <- Biostrings::reverseComplement(gene_seq_dna)
    obs_seq      <- Biostrings::reverseComplement(obs_seq)
    for (i in seq_len(n_comp)) {
      ref_bases[i] <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_bases[i])))
      mut_out[i]   <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(mutation_vec[i])))
    }
    snp_in_gene <- (gene_end - position_vec) + 1
  }

  aa_ref <- Biostrings::translate(gene_seq_dna, if.fuzzy.codon = "solve")
  aa_obs <- Biostrings::translate(obs_seq,       if.fuzzy.codon = "solve")

  mutations <- Biostrings::pairwiseAlignment(aa_ref, aa_obs) %>%
    Biostrings::mismatchTable() %>%
    mutate(AA_Change = paste0(gene_name, ":", PatternSubstring, PatternStart, SubjectSubstring))

  data.frame(
    SNP                   = as.character(genome_snp),
    snp_position_genome   = paste(position_vec,  collapse = "|"),
    snp_position_gene     = paste(snp_in_gene,   collapse = "|"),
    Theoretical_Reference = paste(ref_bases,     collapse = "|"),
    Gene                  = as.character(gene_name),
    Product               = as.character(gene_product),
    AA = as.character(ifelse(
      length(as.character(unique(mutations$AA_Change))) > 0,
      as.character(unique(mutations$AA_Change)),
      "Synonymous"
    )),
    SNP_Gene = paste(paste0(ref_bases, snp_in_gene, mut_out), collapse = "|")
  )
}

.translate_insertion <- function(position, mutation, genome_snp,
                                  gene_seq, gene_start, gene_end, gene_strand,
                                  gene_name, gene_product) {
  snp_in_gene  <- (position - gene_start) + 1
  ref_gene_str <- gene_seq
  mut_out      <- mutation
  theo_ref     <- substr(ref_gene_str, snp_in_gene, snp_in_gene)

  if (gene_strand == "-") {
    ref_gene_str <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_gene_str)))
    snp_in_gene  <- (gene_end - position) + 1
    mut_out      <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(mutation)))
    theo_ref     <- substr(ref_gene_str, snp_in_gene, snp_in_gene)
  }

  obs_gene_str <- paste0(substr(ref_gene_str, 1, snp_in_gene - 1),
                          mut_out,
                          substr(ref_gene_str, snp_in_gene + 1, nchar(ref_gene_str)))

  aa_ref_str <- as.character(Biostrings::translate(Biostrings::DNAString(ref_gene_str), if.fuzzy.codon = "solve"))
  aa_obs_str <- as.character(Biostrings::translate(Biostrings::DNAString(obs_gene_str),  if.fuzzy.codon = "solve"))

  ins_bases <- substr(mut_out, 2, nchar(mut_out))

  if (aa_ref_str == aa_obs_str) {
    aa_change <- "Synonymous"
  } else {
    n_min      <- min(nchar(aa_ref_str), nchar(aa_obs_str))
    first_diff <- n_min + 1
    for (.i in seq_len(n_min)) {
      if (substr(aa_ref_str, .i, .i) != substr(aa_obs_str, .i, .i)) { first_diff <- .i; break }
    }
    ins_len_aa <- nchar(aa_obs_str) - nchar(aa_ref_str)
    ins_aa_str <- substr(aa_obs_str, first_diff, first_diff + ins_len_aa - 1)
    aa_change  <- paste0(gene_name, ":", first_diff, "ins", ins_aa_str)
  }

  data.frame(
    SNP                   = as.character(genome_snp),
    snp_position_genome   = as.character(position),
    snp_position_gene     = as.character(snp_in_gene),
    Theoretical_Reference = theo_ref,
    Gene                  = as.character(gene_name),
    Product               = as.character(gene_product),
    AA                    = aa_change,
    SNP_Gene              = paste0(theo_ref, snp_in_gene, "ins", ins_bases)
  )
}

.translate_deletion <- function(position_vec, genome_snp,
                                 gene_seq, gene_start, gene_end, gene_strand,
                                 gene_name, gene_product) {
  n_comp       <- length(position_vec)
  snp_in_gene  <- (position_vec - gene_start) + 1
  ref_gene_str <- gene_seq
  theo_refs    <- sapply(snp_in_gene, function(p) substr(ref_gene_str, p, p))

  if (gene_strand == "-") {
    ref_gene_str <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_gene_str)))
    snp_in_gene  <- sort((gene_end - position_vec) + 1)
    theo_refs    <- sapply(snp_in_gene, function(p) substr(ref_gene_str, p, p))
  }

  del_start    <- min(snp_in_gene)
  del_end      <- max(snp_in_gene)
  obs_gene_str <- paste0(substr(ref_gene_str, 1, del_start - 1),
                          substr(ref_gene_str, del_end + 1, nchar(ref_gene_str)))

  aa_ref_str <- as.character(Biostrings::translate(Biostrings::DNAString(ref_gene_str), if.fuzzy.codon = "solve"))
  aa_obs_str <- as.character(Biostrings::translate(Biostrings::DNAString(obs_gene_str),  if.fuzzy.codon = "solve"))

  del_aa_start <- ceiling(del_start / 3)
  del_aa_end   <- ceiling(del_end   / 3)
  del_aa_str   <- substr(aa_ref_str, del_aa_start, del_aa_end)

  aa_change <- if (aa_ref_str == aa_obs_str) "Synonymous" else
    paste0(gene_name, ":", del_aa_str, del_aa_start, "del")

  data.frame(
    SNP                   = as.character(genome_snp),
    snp_position_genome   = paste(position_vec, collapse = "|"),
    snp_position_gene     = paste(snp_in_gene,  collapse = "|"),
    Theoretical_Reference = paste(theo_refs,    collapse = "|"),
    Gene                  = as.character(gene_name),
    Product               = as.character(gene_product),
    AA                    = aa_change,
    SNP_Gene              = paste0(theo_refs[1], del_start, "del", n_comp, "bp")
  )
}

# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

snps.to.amino <- function(snp_db, ref_seq, cores = NULL) {
  reference    <- suppressWarnings(genbankr::readGenBank(ref_seq))
  Reference_DF <- data.frame(reference@cds)

  Reference_DF <- as.data.frame(lapply(Reference_DF, function(x) if (is.list(x)) sapply(x, paste, collapse = ";") else x))
  Reference_DF$experiment <- NULL
  Reference_DF$inference  <- NULL
  Reference_DF <- Reference_DF %>%
    mutate(across(where(~inherits(.x, "CharacterList")), ~sapply(.x, paste, collapse = "; ")))

  Reference_DF$sequence <- "No Seq"
  for (i in 1:nrow(Reference_DF)) {
    Reference_DF$sequence[i] <- substr(as.character(reference@sequence),
                                       Reference_DF[i, 2], Reference_DF[i, 3])
  }

  # SNP strings are 1 or 2 "<ref><pos><mut>" tokens, pipe-joined for
  # codon-merge combo rows (e.g. "T5118A|T5119A"); parse each token
  # separately with the existing reference/position/mutation regex, in
  # token order (ascending genome position, per expand_codon_merges()).
  SNP_Tokens <- snp_db %>%
    mutate(.snp_row = row_number()) %>%
    select(.snp_row, SNP) %>%
    mutate(token = str_split(SNP, "\\|")) %>%
    unnest(token) %>%
    separate(token, into = c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = FALSE) %>%
    separate(snp_position, into = c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)") %>%
    mutate(snp_position = as.numeric(snp_position)) %>%
    select(-token)

  Tokens_By_Row <- split(SNP_Tokens, SNP_Tokens$.snp_row)

  Amino_Acid_List <- snp_db %>%
    mutate(.snp_row = row_number()) %>%
    left_join(
      SNP_Tokens %>%
        group_by(.snp_row) %>%
        summarise(snp_mutation = paste(snp_mutation, collapse = "|"), .groups = "drop"),
      by = ".snp_row"
    ) %>%
    select(-.snp_row)

  if (!is.null(cores) && cores > 1) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  } else {
    registerDoSEQ()
  }

  Temp <- foreach(SNP = 1:nrow(Amino_Acid_List), .combine = rbind,
                  .export  = c(".translate_snp", ".translate_insertion", ".translate_deletion"),
                  .packages = c("dplyr", "Biostrings")) %dopar% {

    GENOME_SNP   <- Amino_Acid_List$SNP[SNP]
    Components   <- Tokens_By_Row[[as.character(SNP)]]
    POSITION_vec <- Components$snp_position
    MUTATION_vec <- Components$snp_mutation
    n_comp       <- length(POSITION_vec)

    Out <- data.frame()

    is_snp       <- all(MUTATION_vec != "_" & nchar(MUTATION_vec) == 1)
    is_insertion <- length(MUTATION_vec) == 1 && nchar(MUTATION_vec[1]) %% 3 == 1 && nchar(MUTATION_vec[1]) > 1
    is_deletion  <- all(MUTATION_vec == "_") && length(MUTATION_vec) %% 3 == 0

    if (is_snp) {
      for (GENE in 1:nrow(Reference_DF)) {
        if (all(POSITION_vec >= Reference_DF$start[GENE] & POSITION_vec <= Reference_DF$end[GENE])) {
          Out <- rbind(Out, .translate_snp(
            POSITION_vec, MUTATION_vec, GENOME_SNP,
            Reference_DF$sequence[GENE], Reference_DF$start[GENE],
            Reference_DF$end[GENE],      Reference_DF$strand[GENE],
            Reference_DF$gene[GENE],     Reference_DF$product[GENE]
          ))
        }
      }
    } else if (is_insertion) {
      for (GENE in 1:nrow(Reference_DF)) {
        if (POSITION_vec[1] >= Reference_DF$start[GENE] & POSITION_vec[1] <= Reference_DF$end[GENE]) {
          Out <- rbind(Out, .translate_insertion(
            POSITION_vec[1], MUTATION_vec[1], GENOME_SNP,
            Reference_DF$sequence[GENE], Reference_DF$start[GENE],
            Reference_DF$end[GENE],      Reference_DF$strand[GENE],
            Reference_DF$gene[GENE],     Reference_DF$product[GENE]
          ))
        }
      }
    } else if (is_deletion) {
      for (GENE in 1:nrow(Reference_DF)) {
        if (all(POSITION_vec >= Reference_DF$start[GENE] & POSITION_vec <= Reference_DF$end[GENE])) {
          Out <- rbind(Out, .translate_deletion(
            POSITION_vec, GENOME_SNP,
            Reference_DF$sequence[GENE], Reference_DF$start[GENE],
            Reference_DF$end[GENE],      Reference_DF$strand[GENE],
            Reference_DF$gene[GENE],     Reference_DF$product[GENE]
          ))
        }
      }
    }
    Out
  }

  Out <- full_join(Amino_Acid_List, Temp)

  Out$snp_mutation <- as.character(unlist(Out$snp_mutation))
  mut_components <- str_split(Out$snp_mutation, "\\|")
  has_insertion  <- sapply(mut_components, function(m) any(nchar(m) > 1))
  has_deletion   <- sapply(mut_components, function(m) any(m == "_"))

  is_inframe_ins <- sapply(mut_components, function(m) length(m) == 1 && nchar(m[1]) > 1 && nchar(m[1]) %% 3 == 1)
  is_inframe_del <- sapply(mut_components, function(m) all(m == "_") && length(m) %% 3 == 0)

  Out$AA[is.na(Out$AA) & has_insertion & !is_inframe_ins] <- "Insertion Not In-frame"
  Out$AA[is.na(Out$AA) & has_deletion  & !is_inframe_del] <- "Deletion Not In-frame"
  Out$AA[is.na(Out$AA)]                                   <- "Non-coding SNP"

  Out$SNP_Gene[is.na(Out$SNP_Gene) & has_insertion & !is_inframe_ins] <- "Insertion Not In-frame"
  Out$SNP_Gene[is.na(Out$SNP_Gene) & has_deletion  & !is_inframe_del] <- "Deletion Not In-frame"
  Out$SNP_Gene[is.na(Out$SNP_Gene)]                                    <- "Non-coding SNP"

  select(Out, SNP, SNP_Gene, AA, Gene, Product, Theoretical_Reference)
}
