snps.to.amino <- function(snp_db, ref_seq, cores = parallelly::availableCores()) {
  library(tidyverse)
  library(genbankr)
  library(Biostrings)
  library(foreach)
  library(doParallel)
  library(parallelly)

  reference    <- suppressWarnings(genbankr::readGenBank(ref_seq))
  Reference_DF <- data.frame(reference@cds)

  # Convert Bioconductor CharacterLists to standard character vectors
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

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Temp <- foreach(SNP = 1:nrow(Amino_Acid_List), .combine = rbind) %dopar% {
    library(dplyr)
    library(Biostrings)

    GENOME_SNP   <- Amino_Acid_List$SNP[SNP]
    Components   <- Tokens_By_Row[[as.character(SNP)]]
    POSITION_vec <- Components$snp_position
    MUTATION_vec <- Components$snp_mutation
    n_comp       <- length(POSITION_vec)

    Out <- data.frame()

    if (all(MUTATION_vec != "_" & nchar(MUTATION_vec) == 1)) {
      for (GENE in 1:nrow(Reference_DF)) {
        if (all(POSITION_vec >= Reference_DF$start[GENE] & POSITION_vec <= Reference_DF$end[GENE])) {

          SNP_in_gene_vec <- (POSITION_vec - Reference_DF$start[GENE]) + 1
          Reference_Seq   <- Biostrings::DNAString(Reference_DF$sequence[GENE])
          Observed_Seq    <- Biostrings::DNAString(Reference_DF$sequence[GENE])

          Theoretical_Ref_vec <- character(n_comp)
          for (i in seq_len(n_comp)) {
            Theoretical_Ref_vec[i] <- as.character(Reference_Seq[SNP_in_gene_vec[i]])
            Observed_Seq[SNP_in_gene_vec[i]] <- MUTATION_vec[i]
          }

          MUTATION_out_vec <- MUTATION_vec

          if (Reference_DF$strand[GENE] == "-") {
            Reference_Seq <- Biostrings::reverseComplement(Reference_Seq)
            Observed_Seq  <- Biostrings::reverseComplement(Observed_Seq)
            for (i in seq_len(n_comp)) {
              Theoretical_Ref_vec[i] <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(Theoretical_Ref_vec[i])))
              MUTATION_out_vec[i]    <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(MUTATION_vec[i])))
            }
            SNP_in_gene_vec <- (Reference_DF$end[GENE] - POSITION_vec) + 1
          }

          AA_Seq      <- Biostrings::translate(Reference_Seq, if.fuzzy.codon = "solve")
          AA_Observed <- Biostrings::translate(Observed_Seq,  if.fuzzy.codon = "solve")

          Mutations <- Biostrings::pairwiseAlignment(AA_Seq, AA_Observed) %>%
            Biostrings::mismatchTable() %>%
            mutate(AA_Change = paste(Reference_DF$gene[GENE], ":", PatternSubstring,
                                     PatternStart, SubjectSubstring, sep = ""))

          SNP_Gene_vec <- paste0(Theoretical_Ref_vec, SNP_in_gene_vec, MUTATION_out_vec)

          Out <- rbind(Out, data.frame(
            SNP                   = as.character(GENOME_SNP),
            snp_position_genome   = paste(POSITION_vec, collapse = "|"),
            snp_position_gene     = paste(SNP_in_gene_vec, collapse = "|"),
            Theoretical_Reference = paste(Theoretical_Ref_vec, collapse = "|"),
            Gene                  = as.character(Reference_DF$gene[GENE]),
            Product               = as.character(Reference_DF$product[GENE]),
            AA                    = as.character(ifelse(
              length(as.character(unique(Mutations$AA_Change))) > 0,
              as.character(unique(Mutations$AA_Change)),
              "Synonymous"
            )),
            SNP_Gene = paste(SNP_Gene_vec, collapse = "|")
          ))
        }
      }
    }
    Out
  }

  stopCluster(cl)

  Out <- full_join(Amino_Acid_List, Temp)

  Out$snp_mutation <- as.character(unlist(Out$snp_mutation))
  mut_components <- str_split(Out$snp_mutation, "\\|")
  has_insertion  <- sapply(mut_components, function(m) any(nchar(m) > 1))
  has_deletion   <- sapply(mut_components, function(m) any(m == "_"))

  Out$AA[is.na(Out$AA) & has_insertion] <- "Insertions Not Supported"
  Out$AA[is.na(Out$AA) & has_deletion]  <- "Deletions Not Supported"
  Out$AA[is.na(Out$AA)]                 <- "Non-coding SNP"

  Out$SNP_Gene[is.na(Out$SNP_Gene) & has_insertion] <- "Insertions Not Supported"
  Out$SNP_Gene[is.na(Out$SNP_Gene) & has_deletion]  <- "Deletions Not Supported"
  Out$SNP_Gene[is.na(Out$SNP_Gene)]                  <- "Non-coding SNP"

  select(Out, SNP, SNP_Gene, AA, Gene, Product, Theoretical_Reference)
}
