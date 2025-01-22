#!/usr/bin/env Rscript

################################
### ---- Some functions ---- ###
################################

# Function to compute observed/expected ratio
calculate_ratio_dinucleotides_freq <- function(observed, expected) {
  ratio <- observed / expected
  return(ratio)
}

# Function to compute expected dinucleotide frequencies
calculate_expected_dinucleotides_freq <- function(data, nucleotides, dinucleotides) {
  expected <- data[dinucleotides]  # Copy for identical structure
  for (dinuc in dinucleotides) {
    nuc1 <- substr(dinuc, 1, 1)  # First nucleotide
    nuc2 <- substr(dinuc, 2, 2)  # Second nucleotide
    expected[[dinuc]] <- data[[nuc1]] * data[[nuc2]]  # Produce frequencies
  }
  return(expected)
}

# Function to calculate expected dinuc freq dataframe and ratio
# observed/expected dinucleotides frequency as dataframes
calculate_exp_ratio_dinuc_freq <- function(total_data, NUCLEOTIDES, DINUCLEOTIDES){
  exp_ratio_dinuc_freq_list <- list()
  data <- total_data[,c(NUCLEOTIDES, DINUCLEOTIDES)]
  expected_dinuc_df<- calculate_expected_dinucleotides_freq(
    data, NUCLEOTIDES, DINUCLEOTIDES
  )
  exp_ratio_dinuc_freq_list[["expected_dinuc_df"]] <- expected_dinuc_df
  exp_ratio_dinuc_freq_list[["ratio_dinuc_df"]] <- calculate_ratio_dinucleotides_freq(
    data[DINUCLEOTIDES], expected_dinuc_df
  )
  return(exp_ratio_dinuc_freq_list)
}

######################
### ---- Main ---- ###
######################

# Produces "ratio_dinuc_freq_df.csv" file in outdir with Dinucleotides relative 
# abundance per sequence.
dinucleotide_relative_abundance <- function(
    total_matrix,
    NUCLEOTIDES,
    DINUCLEOTIDES,
    outdir
){
  # compute expected dinuc freq and ratio obs/exp dinuc freq
  exp_ratio_dinuc_list <- calculate_exp_ratio_dinuc_freq(
    total_matrix, NUCLEOTIDES, DINUCLEOTIDES
  )
  
  exp_dinuc_freq_df <- exp_ratio_dinuc_list$expected_dinuc_df
  ratio_dinuc_freq_df <- exp_ratio_dinuc_list$ratio_dinuc_df
  write.csv(ratio_dinuc_freq_df, file.path(outdir, "ratio_dinuc_freq_df.csv"))
  
  return(ratio_dinuc_freq_df)
}

