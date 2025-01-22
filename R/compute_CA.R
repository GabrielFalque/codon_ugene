#!/usr/bin/env Rscript

# Computes Correspondance Analysis (CA) on RSCU values from all sequences.
# Creates chart of CA with points colored by gene so it is possible to see
# different clusters. It also plots variables contribution to the first 2 
# dimensions. And plots how variance is explained by the 10 first dimensions.
# Produces two .csv files :
# - "correlation_features_dimensions.csv" : correlation between every feature
#   and the two first CA dimensions along with p-value.
# - "significant_features_correlations.csv" : correlation between every pair of
#   features for which p-value < 0.05 (significant)

################################
### ---- Some functions ---- ###
################################

compute_genes_CA <- function(total_matrix, gene_color_map, ca_output_repo){
  results <- list()
  
  all_features <- c(FEATURES_TO_USE, QUANTI_SUPP_FEATURES)
  
  # Select columns of interest
  wanted_data <- total_matrix[, FEATURES_TO_USE]
  
  # Remove lines with infinite values or NaN
  wanted_data <- wanted_data[complete.cases(wanted_data), ]
  
  # Correspondence Analysis (CA)
  res_cca <- cca(wanted_data)

  # Get points score (sites) and variables (species)
  points_scores <- as.data.frame(scores(res_cca, display = "sites"))
  variables_scores <- as.data.frame(scores(res_cca, display = "species"))
  
  # Clean lines name
  rownames(total_matrix) <- trimws(as.character(rownames(total_matrix)))
  rownames(points_scores) <- trimws(as.character(rownames(points_scores)))
  
  # Add "gene" column
  points_scores$gene <- total_matrix[rownames(points_scores), "gene"]
  
  # Clean genes name
  points_scores$gene <- trimws(as.character(points_scores$gene))  
  unique_genes <- unique(points_scores$gene)

  # Verify if gene names are matching with color palette
  missing_genes <- setdiff(unique_genes, names(gene_color_map))
  if (length(missing_genes) > 0) {
    stop("Some gene values in points_scores are not present in gene_color_map : ", 
         paste(missing_genes, collapse = ", "))
  }

  # Define output SVG file
  svg(paste0(ca_output_repo, "ca_clustered.svg"), 
      width = 7, height = 5)  # Adjust dimensions
  
  # Adjust the margins to leave more room for the caption
  par(mar = c(2.5, 2.5, 5, 8))  # (bottom, left, top, right)
  
  # Create the Correspondence Analysis plot
  plot(res_cca, type = "n", main = "Correspondance Analysis")

  # Add points for each gene using gene_color_map
  for (gene_name in names(gene_color_map)) {
    points_subset <- subset(points_scores, gene == gene_name)
    points(points_subset[, 1], points_subset[, 2], 
           col = gene_color_map[[gene_name]], pch = 19, cex = 0.8)
  }
  
  # Add legend on the right
  legend("topright", inset = c(-0.3, 0),  # Adjust legend position
         legend = names(gene_color_map), col = gene_color_map, 
         pch = 19, cex = 1,               # Adjust legend size
         xpd = TRUE)                      # Allows the legend to extend beyond 
                                          # the graphics area

  # Close graphic tool
  dev.off()
 
  # Extract variables score
  var_scores <- as.data.frame(scores(res_cca, display = "species"))
  
  # Compute variables contribution to the first two dimensions
  contrib_dim1 <- abs(var_scores[, 1]) / sum(abs(var_scores[, 1])) * 100
  contrib_dim2 <- abs(var_scores[, 2]) / sum(abs(var_scores[, 2])) * 100
  
  # Create dataframe with contributions
  contrib_df <- data.frame(
    Variable = rownames(var_scores),
    Contribution_Dim1 = contrib_dim1,
    Contribution_Dim2 = contrib_dim2
  )
  
  # Select 10 most contributing variables for each dimension
  top10_dim1 <- contrib_df[
    order(contrib_df$Contribution_Dim1, decreasing = TRUE), ][1:10, ]
  top10_dim2 <- contrib_df[
    order(contrib_df$Contribution_Dim2, decreasing = TRUE), ][1:10, ]
  
  # Barplot for dimension 1
  plot_dim1 <- ggplot(top10_dim1, aes(x = reorder(Variable, Contribution_Dim1), y = Contribution_Dim1)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +
    labs(title = "Top 10 Variables Contribuant à la Dimension 1", y = "Contribution (%)", x = "Variable")
  
  # Save plot as SVG
  ggsave(paste0(ca_output_repo, "top10_dim1.svg"), plot = plot_dim1, width = 10, height = 7)
  results[["plot_dim1"]] <- plot_dim1
  
  # Barplot for dimension 2
  plot_dim2 <- ggplot(top10_dim2, aes(x = reorder(Variable, Contribution_Dim2), y = Contribution_Dim2)) +
    geom_bar(stat = "identity", fill = "lightgreen") +
    coord_flip() +
    labs(title = "Top 10 Variables Contribuant à la Dimension 2", y = "Contribution (%)", x = "Variable")
  
  # Save plot as SVG
  ggsave(paste0(ca_output_repo, "top10_dim2.svg"), plot = plot_dim2, width = 10, height = 7)
  results[["plot_dim2"]] <- plot_dim2
  
  # Extract eigenvalues and compute explained variance
  eigenvalues <- res_cca$CA$eig
  var_explained <- eigenvalues / sum(eigenvalues) * 100
  
  # Create dataframe for 10 first dimensions
  var_df <- data.frame(
    Dimension = c("1","2","3","4","5","6","7","8","9","10"),
    Variance_Explained = var_explained[1:10]
  )
  
  # add percentage to display
  var_df <- var_df |> 
    dplyr::mutate(perc = paste0(round(Variance_Explained, 1), "%"))
  
  # Braplot of variance explained
  variance_explained_plot <- ggplot(var_df, aes(x = Dimension, y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "orange") +
    scale_x_discrete(limits=var_df$Dimension)+
    ## add percentage labels
    geom_text(aes(label = perc), 
              hjust = 0.5, nudge_y = +.7) + 
    labs(title = "Variance explained by the 10 first dimensions", 
         y = "Variance explained (%)", 
         x = "Dimension") 
  
  # Save plot as SVG
  ggsave(paste0(ca_output_repo, "variance_explained.svg"), 
         plot = variance_explained_plot, width = 10, height = 7)
  results[["variance_explained_plot"]] <- variance_explained_plot
  
  ### Correlation between features
  
  # Extract observation scores for first 2 dimensions
  obs_scores <- as.data.frame(scores(res_cca, display = "sites"))
  
  # Extract numerical variables of total_matrix (excluding 'gene')
  numeric_vars <- total_matrix[, sapply(total_matrix, is.numeric)]
  
  # Initialize list to save results
  cor_results <- list()
  
  # Function for computing Pearson correlation and p-value
  get_cor_and_pval <- function(var, axis_scores) {
    result <- cor.test(var, axis_scores, method = "pearson")
    list(correlation = result$estimate, p_value = result$p.value)
  }
  
  # Loop to compute correlation between each variable and each dimension
  for (var_name in colnames(numeric_vars)) {
    var_data <- numeric_vars[[var_name]]
    
    # Correlation with Dim1
    cor_dim1 <- get_cor_and_pval(var_data, obs_scores[, 1])
    
    # Correlation with Dim2
    cor_dim2 <- get_cor_and_pval(var_data, obs_scores[, 2])
    
    # Add results to the list
    cor_results[[var_name]] <- c(
      Correlation_Dim1 = cor_dim1$correlation,
      P_Value_Dim1 = cor_dim1$p_value,
      Correlation_Dim2 = cor_dim2$correlation,
      P_Value_Dim2 = cor_dim2$p_value
    )
  }
  
  # Convert results to dataframe
  cor_df <- do.call(rbind, cor_results)
  cor_df <- as.data.frame(cor_df)
  
  # Function to compute correlation matrix and p-values between features 
  calculate_correlation_between_features <- function(data) {
    vars <- colnames(data)
    results <- data.frame(
      var1 = character(),
      var2 = character(),
      correlation = numeric(),
      p_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:(ncol(data) - 1)) {
      for (j in (i + 1):ncol(data)) {
        x <- data[[i]]
        y <- data[[j]]
        test <- cor.test(x, y, method = "pearson")
        results <- rbind(results, data.frame(
          var1 = vars[i],
          var2 = vars[j],
          correlation = test$estimate,
          p_value = test$p.value
        ))
      }
    }
    return(results)
  }
  
  # Compute correlation and p-value between features
  correlation_results <- calculate_correlation_between_features(numeric_vars)
  
  # Keep only pair with p-value < 0.05
  significant_features_corr <- correlation_results %>%
    filter(p_value < 0.05)
  
  # Save significative results in a CSV file
  write.csv(significant_features_corr, 
            paste0(ca_output_repo,"significant_features_correlations.csv"), 
            row.names = FALSE)
  
  # Save correlation between features and dimensions in CSV file
  write.csv(cor_df,paste0(ca_output_repo,"correlation_features_dimensions.csv"))
  
  # Return results
  results[["correlation_features_dimensions"]] <- cor_df
  results[["numeric_vars"]] <- numeric_vars # ADDED FOR TEST
  
  return(results)
}

######################
### ---- Main ---- ###
######################

# Computes Correspondance Analysis (CA) on RSCU values from all sequences.
# Creates chart of CA with points colored by gene so it is possible to see
# different clusters. It also plots variables contribution to the first 2 
# dimensions. And plots how variance is explained by the 10 first dimensions.
# Produces two .csv files :
# - "correlation_features_dimensions.csv" : correlation between every feature
#   and the two first CA dimensions along with p-value.
# - "significant_features_correlations.csv" : correlation between every pair of
#   features for which p-value < 0.05 (significant)

compute_CA <- function(
    total_matrix,
    gene_color_map,
    outdir
){
  # ---- CA computation ----
  
  ca_results <- compute_genes_CA(
    total_matrix,
    gene_color_map,
    outdir
    )
  
  return(ca_results)
}