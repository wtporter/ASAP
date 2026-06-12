genome.snp.to.gene.snp <- function(snp_db, ref_seq, cores = parallelly::availableCores()) {
  library(tidyverse)
  library(genbankr)
  library(Biostrings)
  library(foreach)
  library(doParallel)
  library(parallelly)

  reference    <- suppressWarnings(genbankr::readGenBank(ref_seq))
  Reference_DF <- left_join(
    data.frame(reference@genes),
    data.frame(reference@cds) %>% select(locus_tag, product, translation)
  )

  Reference_DF$sequence <- "No Seq"
  for (i in 1:nrow(Reference_DF)) {
    Reference_DF$sequence[i] <- substr(as.character(reference@sequence),
                                       Reference_DF[i, 2], Reference_DF[i, 3])
  }

  Reference_DF <- Reference_DF %>%
    mutate(gene = ifelse(is.na(gene), locus_tag, gene))

  SNP_List <- snp_db

  # SNP strings are 1 or 2 "<ref><pos><mut>" tokens, pipe-joined for
  # codon-merge combo rows (e.g. "T5118A|T5119A"); parse each token
  # separately with the existing reference/position/mutation regex, in
  # token order (ascending genome position, per expand_codon_merges()).
  SNP_Tokens <- SNP_List %>%
    mutate(.snp_row = row_number()) %>%
    select(.snp_row, SNP) %>%
    mutate(token = str_split(SNP, "\\|")) %>%
    unnest(token) %>%
    separate(token, into = c("reference", "snp_position"), sep = "(?<=\\D)(?=\\d)", remove = FALSE) %>%
    separate(snp_position, into = c("snp_position", "snp_mutation"), sep = "(?<=\\d)(?=\\D)") %>%
    mutate(snp_position = as.numeric(snp_position)) %>%
    select(-token)

  Tokens_By_Row <- split(SNP_Tokens, SNP_Tokens$.snp_row)

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  Temp <- foreach(SNP = 1:nrow(SNP_List), .combine = rbind) %dopar% {
    library(dplyr)
    library(Biostrings)

    GENOME_SNP <- SNP_List$SNP[SNP]
    Components <- Tokens_By_Row[[as.character(SNP)]]

    POSITION_vec <- Components$snp_position
    MUTATION_vec <- Components$snp_mutation
    n_comp <- length(POSITION_vec)

    Out <- data.frame()

    for (GENE in 1:nrow(Reference_DF)) {
      if (all(POSITION_vec >= Reference_DF$start[GENE] & POSITION_vec <= Reference_DF$end[GENE])) {

        Reference_Seq   <- Biostrings::DNAString(Reference_DF$sequence[GENE])
        SNP_in_gene_vec <- (POSITION_vec - Reference_DF$start[GENE]) + 1

        Theoretical_Ref_vec <- character(n_comp)
        for (i in seq_len(n_comp)) {
          Theoretical_Ref_vec[i] <- as.character(Reference_Seq[SNP_in_gene_vec[i]])
        }

        MUTATION_out_vec <- MUTATION_vec

        if (Reference_DF$strand[GENE] == "-") {
          for (i in seq_len(n_comp)) {
            Theoretical_Ref_vec[i] <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(Theoretical_Ref_vec[i])))
            MUTATION_out_vec[i] <- if (MUTATION_vec[i] == "_") "_" else
              as.character(Biostrings::reverseComplement(Biostrings::DNAString(MUTATION_vec[i])))
          }
          SNP_in_gene_vec <- (Reference_DF$end[GENE] - POSITION_vec) + 1
        }

        SNP_Gene_vec <- paste0(Theoretical_Ref_vec, SNP_in_gene_vec, MUTATION_out_vec)

        Out <- rbind(Out, data.frame(
          SNP                   = GENOME_SNP,
          snp_position_genome   = paste(POSITION_vec, collapse = "|"),
          snp_position_gene     = paste(SNP_in_gene_vec, collapse = "|"),
          snp_mutation          = paste(MUTATION_out_vec, collapse = "|"),
          Theoretical_Reference = paste(Theoretical_Ref_vec, collapse = "|"),
          Gene                  = Reference_DF$gene[GENE],
          SNP_Gene              = paste(SNP_Gene_vec, collapse = "|")
        ))
      }
    }
    Out
  }

  stopCluster(cl)

  Out <- full_join(select(SNP_List, SNP), Temp)
  Out$Gene[is.na(Out$Gene)]         <- "Non-gene region"
  Out$SNP_Gene[is.na(Out$SNP_Gene)] <- "Non-gene region"

  select(Out, SNP, Gene, SNP_Gene)
}
