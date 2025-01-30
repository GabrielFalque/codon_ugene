#!/usr/bin/env Rscript

# ---- Plots generation ----
# Generates plots for :
# - RSCU distribution : one plot per gene; each plot shows RSCU distribution
#   for each codon, grouped by amino acid.
# - GC plot : shows GC, GC1, GC12, GC2 and GC3 content percentage for each gene
# - ENC plot : shows ENC as a function of GC3. Points colored by gene.
# - Neutrality plot : shows GC12 content as a function of GC3 content.
# - Dinucleotides relative abundance plot : shows relative abundance (computed
#   as (observed frequency of dinucleotide)/(expected frequency)) of each
#   dinucleotide. Points are colored by gene.

################################
### ---- Some functions ---- ###
################################

# save RSCU distribution  per gene in asked repo
display_rscu_distribution <- function(total_matrix, rscu_output_repo) {
  # Filter the necessary columns
  data_filtered <- total_matrix %>% 
    select(all_of(FEATURES_TO_USE), gene)
  
  # Restructure in long format
  data_long <- data_filtered %>%
    pivot_longer(cols = all_of(FEATURES_TO_USE), names_to = "Codon", values_to = "RSCU")
  
  codons_dict_reduced <- CODON_DICT
  #Met, Trp, Stop
  codons_dict_reduced[["Met"]] <- NULL
  codons_dict_reduced[["Trp"]] <- NULL
  codons_dict_reduced[["Stop"]] <- NULL
  
  # Determine codons order and their annotations
  codon_order <- unlist(codons_dict_reduced)
  codon_labels <- rep(names(codons_dict_reduced), 
                      times = sapply(codons_dict_reduced, length))
  
  # Add a column to annotate amino acids
  data_long <- data_long %>%
    mutate(
      AminoAcid = factor(codon_labels[match(Codon, codon_order)], 
                         levels = names(CODON_DICT)),
      Codon = factor(Codon, levels = codon_order) # Codons order
    )

  # Determine amino acid groups mean position on plot
  amino_acid_positions <- data_long %>%
    group_by(AminoAcid) %>%
    summarise(
      Position = mean(as.numeric(Codon)),
      Start = min(as.numeric(Codon)) - 0.5,
      End = max(as.numeric(Codon)) + 0.5,
      .groups = "drop"
    )
  
  data_long <- data_long %>%
    # Compute by gene and AminoAcid
    group_by(gene, AminoAcid) %>%
    # Compute MaxValue and filter abnormal values
    mutate(
      MaxValue = max(RSCU[RSCU > 0.0001], na.rm = TRUE),  
    ) %>%
    ungroup()%>%
    group_by(gene, AminoAcid, Codon) %>%
    mutate(
      MaxRSCU = max(RSCU[RSCU > 0.0001], na.rm = TRUE), 
      MaxCodon = near(MaxRSCU, MaxValue) & MaxValue > 0.0001,
      StyledLabel = ifelse(
        MaxCodon,  # if MaxCodon is TRUE
        sprintf("<b><span style='color:red;'>%s</span></b>", Codon),  # Red and bold
        sprintf("<span style='color:black;'>%s</span>", Codon)  # Black and standard
      )
      ) %>%
    ungroup()
  
  # Create a plot for each gene
  plots <- data_long %>%
    group_by(gene) %>%
    group_split() %>%
    lapply(function(df) {
      
      # Generate a list of styled labels
      styled_labels <- setNames(df$StyledLabel, df$Codon)
      
      ggplot(df, aes(x = Codon, y = RSCU, fill = AminoAcid)) +
        geom_boxplot(outlier.size = 1.5) +
        
        # Add amino acids annotations
        annotate("text", 
                 x = amino_acid_positions$Position, 
                 y = max(df$RSCU, na.rm = TRUE) + 0.5,
                 label = amino_acid_positions$AminoAcid, 
                 size = 4, fontface = "bold") +
        
        # Add vertical lines betwee amino acid groups
        geom_vline(xintercept = amino_acid_positions$End[-nrow(amino_acid_positions)], 
                   color = "black", linetype = "dashed") +
        
        # Modify codons name appearence on x axis 
        scale_x_discrete(labels = styled_labels) +
        labs(title = paste("RSCU distribution per gene", unique(df$gene)),
             x = "Codons (per amino acid)", y = "RSCU") +
        theme_bw() +
        theme(
          legend.position = "none", # Remove legend
          axis.text.x = ggtext::element_markdown(
            angle = 90, hjust = 1, # Labels rotation
            vjust = 0.5, size = 10 # and styles
            ), 
          plot.title = element_text(hjust = 0.5)
        )
    })

  results <- list()
  # Save plots for each gene
  for (i in seq_along(plots)) {
    g <- plots[[i]]
    g_name <- unique(data_long$gene)[i]
    ggsave(paste0(rscu_output_repo,"/RSCU_distribution_", g_name, ".png"), 
           g, width = 12, height = 6)
    results[[g_name]] <- g
  }
  
  return(results)
}

# Generate GC plot : shows GC, GC1, GC12, GC2 and GC3 content percentage 
# for each gene
get_gc_plot <- function(gene_features_mean, output_repo){
  gc_features <- c("gc1", "gc2", "gc3", "gc12", "gc")
  reduced_mean_df <- gene_features_mean[gc_features, ]

  # Convert matrix to long dataframe
  data_long <- data.frame(
    type_gc = rep(rownames(reduced_mean_df), times = ncol(reduced_mean_df)),
    gene = rep(colnames(reduced_mean_df), each = nrow(reduced_mean_df)),
    valeur_gc = as.vector(reduced_mean_df)
  )

  # Create grouped barplot
  gc_plot <- 
    ggplot(data_long, aes(x = gene, y = valeur_gc, fill = type_gc)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_brewer(palette = "Set1", name = "GC feature") +
    labs(
      title = "GC content percentage per gene",
      x = "Gene",
      y = "GC percentage (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),  # Fond blanc
      plot.background = element_rect(fill = "white", color = NA),   # Fond global blanc
      legend.position = "right"
    )

  # Save plot as SVG file
  ggsave(paste0(output_repo, "gc_plot.svg"), plot = gc_plot, width = 7, height = 5)
}


# Produces and save ENC plot : shows ENC as a function of GC3. 
# Points colored by gene.
get_enc_plot <- function(total_matrix, gene_color_map, enc_plot_repo){

  # Draw ENC (Effective Number of Codons) vs GC3
  enc_curv <- function(x) { 2 + x + 29 / (x^2 + (1 - x)^2) } #max = 60.5
  straigth_line <- function(x) {x}

  # get min and max of ENC values
  enc_max <- max(total_matrix$enc)
  enc_min <- min(total_matrix$enc)
  
  enc_plot <- ggplot(total_matrix, aes(x = gc3, y = enc, color = gene)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(stat='function', fun=enc_curv, color='green', 
              linetype='dashed', linewidth=1) +
    xlim(0, 1) +
    ylim(enc_min - 10, min(c(enc_max, 70)) + 10)+
    labs(title = "ENC Plot",
         x = "GC3 Content",
         y = "Effective Number of Codons (ENC)",
         color = "Gene") +
    scale_color_manual(values = gene_color_map) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),   
      legend.position = "right"
    )

  # Enregistrer le graphique au format SVG
  ggsave(paste0(enc_plot_repo, "enc_plot.svg"), plot = enc_plot, 
         width = 7, height = 5)

  return(enc_plot)
}


# Generate Neutrality plot : shows GC12 content as a function of GC3 content.
get_neutrality_plot <- function(total_matrix, 
                                gene_color_map, 
                                neutrality_plot_repo){

  # Draw GC12 vs GC3
  neutrality_plot <- ggplot(total_matrix, aes(x = gc3, y = gc12, color = gene)) +
    geom_point(size = 3, alpha = 0.8)+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "green")+
    labs(title = "Neutrality Plot",
         x = "GC3 Content",
         y = "GC12 Content",
         color = "Gene") +
    scale_color_manual(values = gene_color_map) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    )

  # Save plot as SVG file
  ggsave(paste0(neutrality_plot_repo, "neutrality_plot.svg"), 
         plot = neutrality_plot, width = 7, height = 5)
  
  return(neutrality_plot)
}

# Dinucleotides relative abundance plot : shows relative abundance (computed
# as (observed frequency of dinucleotide)/(expected frequency)) of each
# dinucleotide. Points are colored by gene.
get_dinuc_ratio_plot <- function(ratio_dinuc_freq_df, gene_color_map, outdir){
  # Add column with gene name (extraction before ".")
  ratio_dinuc_freq_df <- ratio_dinuc_freq_df %>%
    rownames_to_column(var = "id") %>%
    mutate(gene = sub("\\..*", "", id)) %>%
    select(-id)
  
  # Compute mean relative abundance for each dinucleotide per gene.
  data_summary <- ratio_dinuc_freq_df %>%
    group_by(gene) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  
  # Convert to long format (for ggpot)
  data_long <- data_summary %>%
    pivot_longer(-gene, names_to = "dinucleotide", values_to = "oe_dinuc_ratio")
  
  # Create points plot
  dinuc_plot <- ggplot(data_long, aes(x = dinucleotide, y = oe_dinuc_ratio, color = gene, group = gene)) +
    geom_point(size = 3) +
    geom_line() +
    labs(
      title = "Observed/Expected dinucleotides frequency ratio per gene",
      x = "Dinucleotide",
      y = "O/E ratio",
      color = "gene"
    ) +
    scale_color_manual(values = gene_color_map) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),  # White background
      plot.background = element_rect(fill = "white", color = NA),   
      legend.position = "right"
    )
  
  # Save plot as SVG file
  ggsave(paste0(outdir, "oe_dinuc_freq_plot.svg"), plot = dinuc_plot, 
         width = 7, height = 5)
}

######################
### ---- Main ---- ###
######################

# ---- Plots generation ----
# Generates plots for :
# - RSCU distribution : one plot per gene; each plot shows RSCU distribution
#   for each codon, grouped by amino acid.
# - GC plot : shows GC, GC1, GC12, GC2 and GC3 content percentage for each gene
# - ENC plot : shows ENC as a function of GC3. Points colored by gene.
# - Neutrality plot : shows GC12 content as a function of GC3 content.
# - Dinucleotides relative abundance plot : shows relative abundance (computed
#   as (observed frequency of dinucleotide)/(expected frequency)) of each
#   dinucleotide. Points are colored by gene.

plots_generation <- function(
    total_matrix,
    gene_features_mean,
    ratio_dinuc_freq_df,
    gene_color_map,
    outdir
    ){

  # ---- RSCU distribution ----
  options(warn = -1) 
  # generate RSCU distribution charts per gene
  rscu_distribution_plots <- display_rscu_distribution(total_matrix, outdir) 
  options(warn = 1)

  # ---- GC Plot -----

  # Generate GC plot : shows GC, GC1, GC12, GC2 and GC3 content percentage 
  # for each gene
  get_gc_plot(gene_features_mean, outdir)

  # ---- ENC Plot ----

  # Produces and save ENC plot : shows ENC as a function of GC3. 
  # Points colored by gene.
  options(warn = -1) 
  enc_plot <- get_enc_plot(total_matrix, gene_color_map, outdir)
  options(warn = 1)
  # ---- Neutrality Plot ----

  # generate Neutrality plot : shows GC12 content as a function of GC3 content.
  neutrality_plot <- get_neutrality_plot(total_matrix, gene_color_map, outdir)

    # ---- Dinucleotides relative abundance Plot ----

  # Dinucleotides relative abundance plot : shows relative abundance (computed
  # as (observed frequency of dinucleotide)/(expected frequency)) of each
  # dinucleotide. Points are colored by gene.
  get_dinuc_ratio_plot(ratio_dinuc_freq_df, gene_color_map, outdir)
}