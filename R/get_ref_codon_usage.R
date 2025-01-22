#!/usr/bin/env Rscript

################################
### ---- Some functions ---- ###
################################

# Function to validate CSV reference file 
validate_csv <- function(file_path, verbose=FALSE) {
  
  # read file
  if (!file.exists(file_path)) {
    stop("Error : Specified reference file does not exist.")
  }
  data <- read.csv(file_path, stringsAsFactors = FALSE, sep="\t")
  
  # Verify if required columns are present
  required_columns <- c("CODON", "Amino.acid", "Fraction")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Error : Following columns are missing :", 
               paste(missing_columns, collapse = ", ")))
  }
  
  # Vérifier que tous les codons sont répertoriés (64 codons)
  # Verify that all 64 codons are indexed
  unique_codons <- unique(data$CODON)
  if (length(unique_codons) != 64) {
    stop(paste("Error : File does not contain 64 unique codons. Found codons :", 
               length(unique_codons)))
  }
  
  # Vérifier le type de données pour chaque colonne
  # Verify data type for each required column
  issues <- list()
  if (!is.character(data$CODON)) {
    issues <- c(issues, "'CODON' column must be string.")
  }
  if (!is.character(data$Amino.acid)) {
    issues <- c(issues, "'Amino acid' column must be string.")
  }
  if (!is.numeric(data$Fraction)) {
    issues <- c(issues, "'FRACTION' column must be of numeric type")
  }
  
  # Stop execution if data type issues are detected
  if (length(issues) > 0) {
    stop(paste("Error : Detected issues in data types :\n", 
               paste(issues, collapse = "\n"),
               "\n You can download RSCU table from the desired species",
               "directly from http://codonstatsdb.unr.edu/index.html."))
  }
  
  # If everything is correct
  if (verbose){
    message("CSV reference file is correct.")
  }
  
  return(data)
}

# get codon weights df from reference file
get_ref_codon_weights <- function(ref_codon_usage_df){
  
  # Initialize codon weights vector
  ref_codon_weights <- numeric(length = nrow(ref_codon_usage_df))
  names(ref_codon_weights) <- ref_codon_usage_df$CODON
  
  # Compute codon weights for each amino acid
  for (aa in unique(ref_codon_usage_df$Amino.acid)) {
    
    # Get codons associated with the amino acid
    aa_codons <- ref_codon_usage_df[ref_codon_usage_df$Amino.acid == aa, ]
    
    # Find maximal frequency (preferred codon)
    max_freq <- max(aa_codons$Fraction)
    
    # Compute weight for each codon (= frequency / (preferred codon frequency))
    ref_codon_weights[aa_codons$CODON] <- aa_codons$Fraction / max_freq
  }
  
  return(ref_codon_weights)
}

######################
### ---- Main ---- ###
######################

# Verify provided reference codon usage csv file then
# computes codons weights in reference species genome
get_ref_codon_usage <- function(
    codon_usage_ref_file,
    verbose){
  
  # ---- Reference codon weights computation ----
  
  ref_codon_weights <- NULL
  
  ref_codon_usage_df <- validate_csv(codon_usage_ref_file, verbose=FALSE)
  ref_codon_weights <- get_ref_codon_weights(ref_codon_usage_df)
  
  ref_codons_infos_list <- list()
  ref_codons_infos_list[["ref_codon_usage_df"]] <- ref_codon_usage_df
  ref_codons_infos_list[["ref_codon_weights"]] <- ref_codon_weights
  
  return(ref_codons_infos_list)
}