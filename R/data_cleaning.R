#!/usr/bin/env Rscript

################################
### ---- Some functions ---- ###
################################

# Function to remove "bad sequences" i.e. sequences with more than 15% of 
# ambiguious nucleotides ("n")
remove_bad_sequences <- function(gene_sequences){
  # Tolerance threshold for "n"s
  threshold <- 0.15  # 15%
  
  # Function to compute sequence "n" percentage
  calc_n_percentage <- function(seq) {
    n_frequency <-  1 -letterFrequency(seq, "ACGT", as.prob=TRUE)
    return(n_frequency)
  }
  
  # Filter sequences with less than 15% "n"
  filtered_sequences <- gene_sequences[sapply(gene_sequences, calc_n_percentage) <= threshold]
  
  return(filtered_sequences)
}

######################
### ---- Main ---- ###
######################

# Bad sequences (>15% ambiguious nucleotides) removal and complete genomes 
# selection.
data_cleaning <- function(genes_sequences_list,
                          no_complete){
  # ---- Bad sequences removal ----

  genes_good_sequences_list <- lapply(genes_sequences_list, remove_bad_sequences)
  
  if (!no_complete){
    # ---- Complete genomes selection ----
    # Identify genomes with one sequence for each gene
    # 1. Extract sequences name for each gene
    # 2. Find common genomes (sequence name) for each gene
    all_genomes <- lapply(genes_good_sequences_list, names)  
    common_genomes <- Reduce(intersect, all_genomes)  
    
    # Filter sequences to keep only common genomes
    filtered_gene_sets <- lapply(genes_good_sequences_list, function(gene_set) {
      gene_set[common_genomes]
    })
    
    # Concatenate sequences of each gene of selected genomes to get full CDS
    complete_genomes <- DNAStringSet(sapply(common_genomes, function(genome) {
      paste0(sapply(filtered_gene_sets, 
                    function(gene_set) as.character(gene_set[genome])), 
             collapse = "")
    }))
    

    # Assign genome names to concatenated sequences
    names(complete_genomes) <- common_genomes
    
    genes_good_sequences_list[["complete"]] <- complete_genomes
  }
  
  return(genes_good_sequences_list)
}