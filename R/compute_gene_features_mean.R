#!/usr/bin/env Rscript

################################
### ---- Some functions ---- ###
################################

# compute features mean per gene as dataframe and can save it
compute_features_mean_per_gene <- function(
    total_matrix,
    save=FALSE,
    file_path=NULL
){
  # Group by sum of multiple columns
  df_mean_features <- total_matrix %>% group_by(gene) %>% 
    summarise(gc1=mean(gc1, na.rm = TRUE),
              gc2= mean(gc2, na.rm = TRUE),
              gc3=mean(gc3, na.rm = TRUE),
              gc12=mean(gc12, na.rm = TRUE),
              rcdi=mean(rcdi, na.rm = TRUE),
              enc=mean(enc, na.rm = TRUE),
              cai=mean(cai, na.rm = TRUE),
              gc=mean(gc, na.rm = TRUE),
              aromatic_prct=mean(aromatic_prct, na.rm = TRUE),
              gravy=mean(gravy, na.rm = TRUE),
              .groups = 'drop') %>%
    as.data.frame() %>% data.frame(row.names = 1) 
  
  #transpose
  df_mean_features <- t(df_mean_features)
  if (save){
    if (is.null(file_path)){
      stop("Indicate file path as file_path=/path/to/file/my_file.txt")
    } else{
      write.csv(df_mean_features, file_path)
    }
  }
  return(df_mean_features)
}

######################
### ---- Main ---- ###
######################

# ---- Mean features value per gene ----

# Generate and save a dataframe with mean values per gene of the 
# following features :
# - gc1
# - gc2
# - gc3
# - gc12
# - rcdi
# - enc
# - cai
# - gc
# - aromatic_prct
# - gravy

compute_gene_features_mean <- function(
    total_matrix,
    outdir
){
  # get mean features value per gene
  gene_features_mean <- compute_features_mean_per_gene(total_matrix)
  write.csv(gene_features_mean, file.path(outdir, "gene_features_mean.csv"))
  
  return(gene_features_mean)
}

