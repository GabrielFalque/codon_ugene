#!/usr/bin/env Rscript

################################
### ---- Some functions ---- ###
################################

# Function to compute CAI from DNAStringSet
calculate_cai <- function(dna_string_set, ref_codon_weights) {
  sapply(dna_string_set, function(seq) {
    seq <- toupper(as.character(seq))
    seq_vector <- str_extract_all(seq, boundary("character"))[[1]]
    cai(seq_vector, w = ref_codon_weights)
  })
}

# Function to compute nucleotides frequency from DNAStringSet
calculate_nuc_freq <- function(dna_string_set) {
  nuc_freq <- alphabetFrequency(dna_string_set, baseOnly = TRUE, 
                                as.prob = TRUE)
  rownames(nuc_freq) <- names(dna_string_set)
  transposed_nuc_freq <- t(nuc_freq)
  return(transposed_nuc_freq)
}

# Function to extract codons and their frequency from DNAStringSet
count_codons_local <- function(gene_sequences) {
  gene_codons_freq_list <- list()
  for (i in seq_along(gene_sequences)) {
    seq_name <- names(gene_sequences)[i]
    codons_freq <- oligonucleotideFrequency(gene_sequences[i], width = 3, step = 3)  # Divise la séquence en codons
    gene_codons_freq_list[[seq_name]] <- codons_freq
  }
  return(gene_codons_freq_list)
}

# Function to compute relative dinucleotides abundancy from DNAStringSet
calculate_dinuc_freq <- function(dna_string_set) {
  dinuc_freq <- oligonucleotideFrequency(dna_string_set, width = 2, as.prob = TRUE)
  rownames(dinuc_freq) <- names(dna_string_set)
  transposed_dinuc_freq <- t(dinuc_freq)
  return(transposed_dinuc_freq)
}

# Function to compute RSCU from one sequence
calculate_rscu <- function(sequence, CODON_DICT, sequence_codons_count) {
  
  # Compute RSCU for each codon
  rscu_values <- numeric()
  CODON_DICT[c('Stop',"Trp","Met")] <- NULL
  for (aa in names(CODON_DICT)) {
    codon_list <- CODON_DICT[[aa]]
    usage <- sequence_codons_count[,codon_list]
    total_usage <- sum(usage, na.rm = TRUE)
    n_codons <- length(codon_list)
    
    # RSCU computation
    rscu <- usage / (total_usage / n_codons)

    # Detect Nan and 0 and replace them by 0.0001
    rscu[is.nan(rscu)] <- 0.0001
    rscu[rscu<=0.0001] <- 0.0001
    rscu_values[codon_list] <- rscu
  }
  return(rscu_values)
}

# Function to compute RSCU for each sequence of a DNAStringSet and return
# a matrix.
calculate_rscu_for_set <- function(dna_set, CODON_DICT, codons_count) {
  
  # Create list to keep results
  rscu_matrix_list <- list()
  
  # Loop over all sequences
  for (i in seq_along(dna_set)) {
    seq_name <- names(dna_set)[i]
    sequence_codons_count <- codons_count[[seq_name]]
    rscu_values <- calculate_rscu(dna_set[[i]], CODON_DICT, sequence_codons_count)
    
    if (!is.null(rscu_values)) {
      rscu_matrix_list[[seq_name]] <- rscu_values
    }
  }
  
  # Transform list into matrix
  rscu_matrix <- do.call(cbind, rscu_matrix_list)
  rownames(rscu_matrix) <- names(rscu_values)
  
  return(rscu_matrix)
}

# Function to compute RCDI from a sequence
calculate_rcdi <- function(codon_counts, codon_weights) {
  
  # Verify codons name
  if (!all(names(codon_counts) %in% names(codon_weights))) {
    stop(cat("All codons present in ‘codon_counts’ must ",
             "exist in ‘codon_weights’."))
  }
  
  # Keep only codons with non null count
  codon_counts_filtered <- codon_counts[,codon_counts > 0]
  
  # Normalize weights to consider only present codons
  relevant_weights <- codon_weights[names(codon_counts_filtered)]
  
  # Compute relative codons frequency
  total_codons <- sum(codon_counts_filtered)
  codon_frequencies <- codon_counts_filtered / total_codons
  
  # Compute RCDI
  rcdi <- prod(relevant_weights ^ codon_frequencies)
  return(rcdi)
}

# Function to compute RCDI from DNAStringSet
calculate_dna_set_rcdi <- function(gene_codons_freq_list, weights){
  rcdi_matrix_list <- list()
  
  # Loop over sequences from DNAStringSet
  for (i in seq_along(gene_codons_freq_list)){
    seq_name <- names(gene_codons_freq_list)[i]
    rcdi_matrix_list[[seq_name]] <- calculate_rcdi(gene_codons_freq_list[[i]], weights)
  }
  
  # Transform list to matrix
  rcdi_matrix <- do.call(cbind, rcdi_matrix_list)
  rownames(rcdi_matrix) <- names("rcdi")
  return(rcdi_matrix)
}

# Function to compute values of GC1, GC2, GC3 and GC12 for each sequence of
# DNAStringSet
calc_GC12_GC3 <- function(dna_set) {
  results <- list()
  gc1_results <- list()
  gc2_results <- list()
  gc3_results <- list()
  gc12_results <- list()
  for (i in seq_along(dna_set)) {
    seq_name <- names(dna_set)[i]
    codon_counts <- oligonucleotideFrequency(dna_set[i], width = 3, step = 3)
    gc12_total = 0
    gc1_total = 0
    gc2_total = 0
    gc3_total = 0
    total_codons = sum(codon_counts)
    
    for (codon in colnames(codon_counts)){
      pos1 <- unlist(strsplit(substr(codon, 1, 1), ""))
      pos2 <- unlist(strsplit(substr(codon, 2, 2), ""))
      pos3 <- unlist(strsplit(substr(codon, 3, 3), ""))
      
      if (pos1 %in% c("C", "G")){
        gc12_total = gc12_total + as.numeric(codon_counts[, codon])
        gc1_total = gc1_total + as.numeric(codon_counts[, codon])
      }
      if (pos2 %in% c("C", "G")){
        gc12_total = gc12_total + as.numeric(codon_counts[, codon])
        gc2_total = gc2_total + as.numeric(codon_counts[, codon])
      }
      if (pos3 %in% c("G", "C")){
        gc3_total = gc3_total + as.numeric(codon_counts[, codon])
      }
    }
    gc12 = gc12_total /(total_codons*2)
    gc1 = gc1_total /total_codons
    gc2 = gc2_total /total_codons
    gc3 = gc3_total /total_codons
    gc1_results[[seq_name]] <- gc1
    gc2_results[[seq_name]] <- gc2
    gc3_results[[seq_name]] <- gc3
    gc12_results[[seq_name]] <- gc12
  }
  
  # Transform lists to matrix
  gc1_matrix <- do.call(cbind, gc1_results)
  rownames(gc1_matrix) <- c("gc1")
  
  gc2_matrix <- do.call(cbind, gc2_results)
  rownames(gc2_matrix) <- c("gc2")
  
  gc3_matrix <- do.call(cbind, gc3_results)
  rownames(gc3_matrix) <- c("gc3")
  
  gc12_matrix <- do.call(cbind, gc12_results)
  rownames(gc12_matrix) <- c("gc12")
  
  results[["gc1_matrix"]] <- gc1_matrix
  results[["gc2_matrix"]] <- gc2_matrix
  results[["gc3_matrix"]] <- gc3_matrix
  results[["gc12_matrix"]] <- gc12_matrix
  
  return(results)
}

# Function to filter Amino acid sequences (AAStringSet)
# 1. Remove sequences with more than 5% of ambiguous amino acids
# 2. For remaining sequences, don't consider ambiguous amino acids ("X")
clean_aa_sequences <- function(aa_set, threshold = 0.05){
  
  # Tolerance threshold for ambiguous amino acids ("X")
  threshold <- 0.05  # 5%
  
  # Function to compute "X" percentage in an AA sequence
  calc_ambiguious_aa <- function(aa_seq) {
    aa_length <- length(aa_seq)
    
    # Count ambiguous AAs ("X")
    ambiguous_count <- sum(letterFrequency(aa_seq, "X"))
    ambiguous_ratio <- ambiguous_count / aa_length
    return(ambiguous_ratio)
  }
  
  # Keep only sequences with less than 5% ambiguous AAs ("X")
  filtered_sequences <- aa_set[
    sapply(aa_set, calc_ambiguious_aa) <= threshold
  ]
  
  # For remaining sequences, remove ambiguous AAS (in order not to disturb 
  # GRAVY and amino acids frequency computations)
  filtered_sequences_no_X <- AAStringSet(sapply(
    filtered_sequences, function(x) gsub("X", "", as.character(x))
  ))
  
  return(filtered_sequences_no_X)
}

# Function to compute the following features for a DNAStringSet :
# - RSCU
# - GC (%)
# - GC1 (%)
# - GC2 (%)
# - GC3 (%)
# - GC12 (%)
# - RCDI
# - ENC
# - CAI
# - Dinucleotides frequency
# - Nucleotides frequency
# - Aromatic amino acids frequency
# - GRAVY score
# Return a matrix with those features value
gene_sequences_features_computation <- function(
    gene_sequences,
    CODON_DICT,
    gene,
    ref_codon_weights) {

  # get codon frequency
  gene_sequences_cf <- cubar::count_codons(gene_sequences, as.prob = TRUE)
  #print(gene)
  #print(gene_sequences)
  # get ENC
  #print("A")
  codon_data <- coRdon::codonTable(gene_sequences)
  #print("B")
  enc<- coRdon::ENC(codon_data)
  enc_matrix <- matrix(unlist(enc), ncol = length(enc), byrow = TRUE)
  rownames(enc_matrix) <- c("enc")
  colnames(enc_matrix) <- names(gene_sequences)
  enc_df <- as.data.frame(enc_matrix)

  # get GC content
  gc <- cubar::get_gc(gene_sequences_cf)
  gc_matrix <- matrix(unlist(gc), ncol = length(gc), byrow = TRUE)
  rownames(gc_matrix) <- c("gc")
  colnames(gc_matrix) <- names(gene_sequences)
  gc_df <- as.data.frame(gc_matrix)

  gc3s <- cubar::get_gc3s(gene_sequences_cf)

  # get CAI
  cai <- calculate_cai(gene_sequences, ref_codon_weights)
  cai_matrix <- matrix(unlist(cai), ncol = length(cai), byrow = TRUE)
  rownames(cai_matrix) <- c("cai")
  colnames(cai_matrix) <- names(gene_sequences)
  cai_df <- as.data.frame(cai_matrix)

  # compute nucleotides frequency
  gene_nuc_freq <- as.data.frame(calculate_nuc_freq(gene_sequences))

  # compute nucleotides frequency
  gene_dinuc_freq <- calculate_dinuc_freq(gene_sequences)

  #count codons
  codons_count <- count_codons_local(gene_sequences)

  # compute RSCU
  rscu_results <- calculate_rscu_for_set(gene_sequences, CODON_DICT, codons_count)

  # Compute RCDI
  rcdi_results <- calculate_dna_set_rcdi(codons_count, ref_codon_weights)
  rownames(rcdi_results) <- c("rcdi")

  # Compute mean gene RSCU
  mean_rscu <- rowMeans(rscu_results)

  # Create Dataframe with all RSCU results
  mean_rscu_df <- data.frame(
    mean_rscu
  )
  colnames(mean_rscu_df) <- gene

  # Compute GC1, GC2, GC3 and GC12 values
  GC12_GC3_results <- calc_GC12_GC3(gene_sequences)
  gc1_matrix <- GC12_GC3_results$gc1_matrix
  gc2_matrix <- GC12_GC3_results$gc2_matrix
  gc3_matrix <- GC12_GC3_results$gc3_matrix
  gc12_matrix <- GC12_GC3_results$gc12_matrix

  # replace gaps with "N"
  gene_sequences_with_N <- gsub("-", "N", as.character(gene_sequences))
  gene_sequences_with_N <- DNAStringSet(gene_sequences_with_N)

  # Traduction from nucleotides to Amino acids
  gene_aa_sequences <- Biostrings::translate(
    gene_sequences_with_N, 
    no.init.codon=FALSE,
    if.fuzzy.codon="X")
  
  # Clean AAs sequences
  # 1. Remove sequences with more than 5% of ambiguous amino acids
  # 2. For remaining sequences, don't consider ambiguous amino acids ("X")
  gene_aa_sequences_cleaned <- clean_aa_sequences(gene_aa_sequences)
  
  # Initialization of a GRAVY scores vector with NAs
  gravy <- rep(NA, length(gene_aa_sequences))
  names(gravy) <- names(gene_aa_sequences)  # Assign sequences name to scores
  
  # Compute GRAVY scores of DNAStringSet cleaned sequences
  gravy_cleaned <- sapply(
    as.character(gene_aa_sequences_cleaned), Peptides::hydrophobicity
  )
  
  # Fill GRAVY vector with computed scores
  gravy[names(gravy_cleaned)] <- gravy_cleaned
  
  # Conversion into a matrix
  gravy_matrix <- matrix(gravy, nrow = 1, byrow = TRUE)
  rownames(gravy_matrix) <- c("gravy")
  colnames(gravy_matrix) <- names(gene_aa_sequences)
  
  # Conversion into a data.frame
  gravy_df <- as.data.frame(gravy_matrix)

  # Aromatic residues percentage computation for every translated sequence
  # Initialization of a Aromatic residues percentage vector with NAs
  aromatic_prcts <- rep(NA, length(gene_aa_sequences))
  # Assign sequences name to scores
  names(aromatic_prcts) <- names(gene_aa_sequences)  

  # Compute aromatic percentage residues for cleaned sequences
  for (i in seq_along(gene_aa_sequences_cleaned)) {
    seq_name <- names(gene_aa_sequences_cleaned)[i]
    arom_tab <- aaComp(as.character(gene_aa_sequences_cleaned[i]))
    aromatic_prcts[seq_name] <- arom_tab[[seq_name]]["Aromatic", "Mole%"]
  }
  
  # Transform vector into a matrix
  aromatic_matrix <- matrix(aromatic_prcts, nrow = 1, byrow = TRUE)
  rownames(aromatic_matrix) <- c("aromatic_prct")
  colnames(aromatic_matrix) <- names(gene_aa_sequences)
  
  # Transform matrix into a Dataframe
  aromatic_df <- as.data.frame(aromatic_matrix)
  
  # Annotate sequences with gene name
  gene_list <- (rep(gene, length(gene_sequences)))
  gene_matrix <- matrix(unlist(gene_list), ncol = length(gene_sequences), byrow = TRUE)
  rownames(gene_matrix) <- c("gene")
  colnames(gene_matrix) <- names(gene_sequences)
  gene_df <- as.data.frame(gene_matrix)

  # All gene sequences features matrix creation
  total_gene_df <- as.data.frame(rbind(
    rscu_results,
    gc1_matrix,
    gc2_matrix,
    gc3_matrix,
    gc12_matrix,
    rcdi_results,
    enc_df,
    cai_df,
    gene_dinuc_freq,
    gene_nuc_freq,
    gc_df,
    aromatic_matrix,
    gravy_df
  ))
  total_gene_df <- rbind(total_gene_df, gene_df)
  
  # Return gene sequences features value and RSCU mean values matrix.
  gene_results_list <- list()
  gene_results_list[["total_gene_df"]] <- total_gene_df
  gene_results_list[["mean_rscu_df"]] <- mean_rscu_df
  
  return(gene_results_list)
}

# Function to compute the following features for all sequences of all genes :
# - RSCU
# - GC (%)
# - GC1 (%)
# - GC2 (%)
# - GC3 (%)
# - GC12 (%)
# - RCDI
# - ENC
# - CAI
# - Dinucleotides frequency
# - Nucleotides frequency
# - Aromatic amino acids frequency
# - GRAVY score
# Return a matrix with those features value
all_genes_matrix_computation <- function(
    all_gene_sequences,
    CODON_DICT,
    ref_codon_weights,
    n_threads,
    scripts_directory) {

  all_genes_list <- list()
  all_genes_rscu_mean_list <- list()
  
  # Configure parallelization plan
  plan(multisession, workers = n_threads)
  
  # Parallelization with future_lapply
  results <- future_lapply(names(all_gene_sequences), function(gene) {
    gene_sequences <- all_gene_sequences[[gene]]
    
    # Call gene_sequences_features_computation function for each gene
    gene_results_list <- gene_sequences_features_computation(
      gene_sequences,
      CODON_DICT,
      gene,
      ref_codon_weights
    )
    
    # Return list with gene name and the two features matrix
    return(list(
      gene = gene,  # Include gene name for traceability
      total_gene_df = gene_results_list$total_gene_df,
      mean_rscu_df = gene_results_list$mean_rscu_df
    ))
  })
  
  # Name lists by gene names
  all_genes_list <- setNames(lapply(
    results, `[[`, "total_gene_df"), sapply(results, `[[`, "gene")
    )
  all_genes_rscu_mean_list <- setNames(lapply(
    results, `[[`, "mean_rscu_df"), sapply(results, `[[`, "gene")
    )
  
  # Verify lines name from each dataframes are the same
  row_names <- lapply(all_genes_list, rownames)
  if (!all(sapply(row_names, function(x) identical(x, row_names[[1]])))) {
    stop("Row names do not match between dataframes.")
  }
  
  # Concatenate by columns
  all_genes_df <- do.call(cbind, all_genes_list)
  all_genes_rscu_mean_df <- do.call(cbind, all_genes_rscu_mean_list)

  # Transpose Dataframe
  translated_matrix_1 <- as.data.frame(t(all_genes_df))
  translated_matrix <- subset(translated_matrix_1, select = -c(gene, other))
  
  # Transform variables from character to numeric
  translated_matrix <- cbind(
    mutate_all(
      translated_matrix, function(x) as.numeric(as.character(x))
    ), translated_matrix_1[,"gene", drop=FALSE]
  )
  
  # Return a dataframe with features value for each sequence of each gene
  # and a dataframe with RSCU mean values for each gene
  all_genes_results_list <- list()
  all_genes_results_list[["total_matrix"]] <- translated_matrix
  all_genes_results_list[["all_genes_rscu_mean_df"]] <- all_genes_rscu_mean_df
  
  return(all_genes_results_list)
}

# Function to get RSCU comparison with reference (human by default)
# Add RSCU values from human genome to our mean RSCU values dataframe
create_RSCU_comparison_file <- function(
    all_genes_rscu_mean_df,
    ref_codon_usage,
    save=FALSE,
    file_path=NULL){
  
  comparison_rscu_df <- all_genes_rscu_mean_df
  comparison_rscu_df$reference <- ref_codon_usage$RSCU[
    match(rownames(comparison_rscu_df), ref_codon_usage$CODON)
  ]
  
  if(save){
    if (is.null(file_path)){
      stop("Indicate file path as file_path=/path/to/file/my_file.txt")
    } else{
      write.csv(comparison_rscu_df, file_path)
    }
  }
  return(comparison_rscu_df)
}

######################
### ---- Main ---- ###
######################

# Compute the following features for all sequences of all genes :
# - RSCU
# - GC (%)
# - GC1 (%)
# - GC2 (%)
# - GC3 (%)
# - GC12 (%)
# - RCDI
# - ENC
# - CAI
# - Dinucleotides frequency
# - Nucleotides frequency
# - Aromatic amino acids frequency
# - GRAVY score
# Compute a Dataframe with those features value as well as a dataframe with RSCU
# mean values for each gene compared with RSCU values from human genome
features_computation <- function(
    all_gene_sequences,
    CODON_DICT,
    ref_codon_weights,
    ref_codon_usage_df,
    n_threads,
    scripts_directory,
    outdir
    ){
  
  # ---- RSCU and features value computation ----

  all_genes_results_list <- all_genes_matrix_computation(
    all_gene_sequences = all_gene_sequences,
    CODON_DICT = CODON_DICT,
    ref_codon_weights = ref_codon_weights,
    n_threads = n_threads,
    scripts_directory = scripts_directory)

  total_matrix <- all_genes_results_list$total_matrix
  all_genes_rscu_mean_df <- all_genes_results_list$all_genes_rscu_mean_df

  # save raw data
  write.csv(total_matrix, file.path(outdir, "total_matrix.csv"))

  # ---- RSCU comparison with ref ----

  comparison_rscu_df <- create_RSCU_comparison_file(all_genes_rscu_mean_df, 
                                                    ref_codon_usage_df)
  
  write.csv(comparison_rscu_df, file.path(outdir, "rscu_with_ref.csv"))
  all_genes_results_list[["comparison_rscu_df"]] <- comparison_rscu_df
  
  return(all_genes_results_list)
}