#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

# Liste des librairies nécessaires
required_packages <- c(
  "cubar", "Biostrings", "ggplot2", "coRdon", "seqinr", "stringr", "tidyverse",
  "ca", "Peptides", "FactoMineR", "factoextra", "vegan", "reshape2",
  "AnaCoDa", "dplyr", "svglite", "optparse", "here", "this.path"
)

cat("Checking libraries...\n")
# Fonction pour installer et charger les packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      invisible(install.packages(pkg, dependencies = TRUE)) #added invisible()
    }
    suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
  }
}

# Appeler la fonction avec la liste des packages
install_and_load(required_packages)
cat("Libraries checked !\n")

cat("Script directory : ", this.dir(), "\n")
data_directory = paste0(this.dir(), "/data/")

################################
###### ---- Parsing ---- ######
################################

# Définir les options en ligne de commande
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = paste0("Directory with gene sequences (mandatory).\n", 
              "Files must be named like gene_{gene_name}.fasta with gene_name", 
              " being Orf1 for example. (Ex : gene_Orf1.fasta"), 
              metavar = "DIR"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Directory where to store output files (mandatory).", 
              metavar = "OUTDIR"),
  make_option(c("-r", "--ref"), type = "character", 
              default = paste0(data_directory, "human_codon_usage.csv"),
              help = paste0("Species codon usage CSV (by default :",
                            "'human_codon_usage.csv'). \nYou can download RSCU",
                            " table from the desired species directly from",
                            " http://codonstatsdb.unr.edu/index.html."), 
              metavar = "FILE"),
  make_option(c("--cov"), action = "store_true", default = FALSE,
              help = "If the sequences are from SARS-CoV-2."),
  make_option(c("--verbose"), action = "store_true", default = FALSE,
              help = "Verbose mode.")
)

# Parser les options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
outdir <- opt$outdir

# Afficher l'aide si aucun argument ou si -h/--help est utilisé
if (is.null(opt$dir) || is.null(outdir) || opt$help) {
  print_help(opt_parser)
  quit(status = 0)
}

# Vérifier que le répertoire d'entrée existe
if (!dir.exists(opt$dir)) {
  stop("Input directory does not exist : ", opt$dir, "\n")
}

# Créer le répertoire de sortie s'il n'existe pas
if (!dir.exists(outdir)) {
  if (opt$verbose) message("Output directory creation : ", outdir,"\n")
  dir.create(outdir, recursive = TRUE)
}

# Récupérer tous les fichiers fasta correspondant au motif
fasta_files <- list.files(path = opt$dir, pattern = "^gene_.*\\.fasta$", full.names = TRUE)
if (length(fasta_files) == 0) {
  stop("No file corresponding to 'gene_*.fasta' pattern in the directory : ", opt$dir)
}

# Indiquer si l'option --cov est activée
if (opt$cov) {
  message("Option --cov activated. Sequences are from SARS-CoV-2.\n")
} else {
  message("Option --cov non activated. Sequences are not from SARS-CoV-2.\n")
}

# Charger le fichier de référence
codon_usage_ref_file <- opt$ref
if (!file.exists(codon_usage_ref_file)) {
  stop("Specified species codon usage file does not exist : ", codon_usage_ref_file)
}
if (opt$verbose) message("Loading species codon usage file : ", codon_usage_ref_file)

# Charger les séquences dans une liste
genes_sequences_list <- list()
for (file in fasta_files) {
  gene_name <- sub("^gene_(.*)\\.fasta$", "\\1", basename(file))  # Extract genes name
  if (opt$verbose) message("Reading gene : ", gene_name)
  
  if (opt$cov && gene_name=="ORF1ab"){
    ORF1ab_sequences <- readDNAStringSet(file)
    raw_ORF1a_sequences <- subseq(ORF1ab_sequences, start = 1, end = 13203)
    raw_ORF1b_sequences <- subseq(ORF1ab_sequences, start = 13203) 
    nogap_ORF1a_sequences<- DNAStringSet(sapply(
      raw_ORF1a_sequences, function(x) gsub("-", "", as.character(x))
    ))
    genes_sequences_list[["ORF1a"]] <- check_cds(
      nogap_ORF1a_sequences,
      check_len = TRUE,
      check_start = TRUE,
      check_stop = FALSE,
      check_istop = TRUE,
      rm_start = TRUE,
      rm_stop = FALSE,
    )
    nogap_ORF1b_sequences <- DNAStringSet(sapply(
      raw_ORF1b_sequences, function(x) gsub("-", "", as.character(x))
    ))
    genes_sequences_list[["ORF1b"]] <- check_cds(
      nogap_ORF1b_sequences,
      check_len = TRUE,
      check_start = FALSE,
      check_stop = TRUE,
      check_istop = TRUE,
      rm_start = FALSE,
      rm_stop = TRUE,
    )
  } else {
    raw_gene_sequences <- readDNAStringSet(file)
    nogap_gene_sequences <- DNAStringSet(sapply(
      raw_gene_sequences, function(x) gsub("-", "", as.character(x))
    ))
    genes_sequences_list[[gene_name]] <- check_cds(nogap_gene_sequences)
  }
}




################################
### ---- Some variables ---- ###
################################

# separation between features to use in Correspondance analysis
features_to_use <- c(
  "GCT",
  "GCC",
  "GCA",
  "GCG",
  "CGT",
  "CGC",
  "CGA",
  "CGG",
  "AGA",
  "AGG",
  "AAT",
  "AAC",
  "GAT",
  "GAC",
  "TGT",
  "TGC",
  "CAA",
  "CAG",
  "GAA",
  "GAG",
  "GGT",
  "GGC",
  "GGA",
  "GGG",
  "CAT",
  "CAC",
  "ATT",
  "ATC",
  "ATA",
  "TTA",
  "TTG",
  "CTT",
  "CTC",
  "CTA",
  "CTG",
  "AAA",
  "AAG",
  #"ATG",
  "TTT",
  "TTC",
  "CCT",
  "CCC",
  "CCA",
  "CCG",
  "TCT",
  "TCC",
  "TCA",
  "TCG",
  "AGT",
  "AGC",
  "ACT",
  "ACC",
  "ACA",
  "ACG",
  #"TGG",
  "TAT",
  "TAC",
  "GTT",
  "GTC",
  "GTA",
  "GTG"#,
  #"TAA",
  #"TAG",
  #"TGA"
)

quanti_supp_features <- c(
  "gc3",
  "gc12",
  "enc",
  "cai_values",
  "AA",
  "AC",
  "AG",
  "AT",
  "CA",
  "CC",
  "CG",
  "CT",
  "GA",
  "GC",
  "GG",
  "GT",
  "TA",
  "TC",
  "TG",
  "TT",
  "A",
  "C",
  "G",
  "T", #"other"
  "gc",
  "aromatic_prct",
  "gravy_scores"
)

# Liste des nucléotides et dinucléotides
nucleotides <- c("A", "C", "G", "T")
dinucleotides <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
                   "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

quali_supp_features <- c("gene")

# Dictionnaire des acides aminés et leurs codons
codon_dict <- list(
  Ala = c("GCT", "GCC", "GCA", "GCG"),
  Arg = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
  Asn = c("AAT", "AAC"),
  Asp = c("GAT", "GAC"),
  Cys = c("TGT", "TGC"),
  Gln = c("CAA", "CAG"),
  Glu = c("GAA", "GAG"),
  Gly = c("GGT", "GGC", "GGA", "GGG"),
  His = c("CAT", "CAC"),
  Ile = c("ATT", "ATC", "ATA"),
  Leu = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
  Lys = c("AAA", "AAG"),
  Met = c("ATG"),
  Phe = c("TTT", "TTC"),
  Pro = c("CCT", "CCC", "CCA", "CCG"),
  Ser = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
  Thr = c("ACT", "ACC", "ACA", "ACG"),
  Trp = c("TGG"),
  Tyr = c("TAT", "TAC"),
  Val = c("GTT", "GTC", "GTA", "GTG"),
  Stop = c("TAA", "TAG", "TGA")
)

################################
### ---- Some functions ---- ###
################################

# get codon weights df from reference file
get_ref_codon_weights <- function(ref_codon_usage_df){
  # Initialiser le vecteur de poids des codons
  ref_codon_weights <- numeric(length = nrow(ref_codon_usage_df))
  names(ref_codon_weights) <- ref_codon_usage_df$CODON
  
  # Calculer les poids des codons
  for (aa in unique(ref_codon_usage_df$Amino.acid)) {
    # Filtrer les codons pour un acide aminé donné
    aa_codons <- ref_codon_usage_df[ref_codon_usage_df$Amino.acid == aa, ]
    
    # Trouver la fréquence maximale (codon préféré)
    max_freq <- max(aa_codons$Fraction)
    
    # Calculer le poids pour chaque codon (fréquence / fréquence du codon préféré)
    ref_codon_weights[aa_codons$CODON] <- aa_codons$Fraction / max_freq
  }
  return(ref_codon_weights)
}

# Fonction pour valider le fichier CSV
validate_csv <- function(file_path, verbose=FALSE) {
  # Lire le fichier
  if (!file.exists(file_path)) {
    stop("Error : Specified reference file does not exist.")
  }
  data <- read.csv(file_path, stringsAsFactors = FALSE, sep="\t")
  
  # Vérifier la présence des colonnes requises
  required_columns <- c("CODON", "Amino.acid", "Fraction")
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0) {
    stop(paste("Error : Following columns are missing :", 
               paste(missing_columns, collapse = ", ")))
  }
  
  # Vérifier que tous les codons sont répertoriés (64 codons)
  unique_codons <- unique(data$CODON)
  if (length(unique_codons) != 64) {
    stop(paste("Error : File does not contain 64 unique codons. Found codons :", length(unique_codons)))
  }
  
  # Vérifier le type de données pour chaque colonne
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
  
  # Arrêter l'exécution si des problèmes de type de données sont détectés
  if (length(issues) > 0) {
    stop(paste("Error : Detected issues in data types :\n", 
               paste(issues, collapse = "\n"),
               "\n You can download RSCU table from the desired species",
               "directly from http://codonstatsdb.unr.edu/index.html."))
  }
  
  # Si tout est correct
  if (verbose){
    message("CSV reference file is correct.")
  }
  
  return(data)
}

# Fonction pour calculer le CAI à partir d'un DNAStringSet
calculate_cai <- function(dna_string_set, ref_codon_weights) {
  sapply(dna_string_set, function(seq) {
    seq <- toupper(as.character(seq))
    seq_vector <- str_extract_all(seq, boundary("character"))[[1]]
    cai(seq_vector, w = ref_codon_weights)
  })
}

# Calculer la fréquence des nucléotides
calculate_nuc_freq <- function(dna_string_set) {
  nuc_freq <- alphabetFrequency(dna_string_set, baseOnly = TRUE, 
                                as.prob = TRUE)
  rownames(nuc_freq) <- names(dna_string_set)
  transposed_nuc_freq <- t(nuc_freq)
  return(transposed_nuc_freq)
}

# Calculer l'abondance relative des dinucléotides
calculate_dinuc_freq <- function(dna_string_set) {
  dinuc_freq <- oligonucleotideFrequency(dna_string_set, width = 2, as.prob = TRUE)
  rownames(dinuc_freq) <- names(dna_string_set)
  transposed_dinuc_freq <- t(dinuc_freq)
  return(transposed_dinuc_freq)
}

# Fonction pour extraire les codons et leur fréquence
count_codons_local <- function(gene_sequences) {
  gene_codons_freq_list <- list()
  for (i in seq_along(gene_sequences)) {
    seq_name <- names(gene_sequences)[i]
    codons_freq <- oligonucleotideFrequency(gene_sequences[i], width = 3, step = 3)  # Divise la séquence en codons
    gene_codons_freq_list[[seq_name]] <- codons_freq
  }
  return(gene_codons_freq_list)
}


# Fonction pour calculer le RSCU pour une séquence
calculate_rscu <- function(sequence, codon_dict, sequence_codons_count) {
  # Compte la fréquence des codons
  # codon_counts <- oligonucleotideFrequency(sequence, width = 3, step = 3)
  
  # Calculer le RSCU pour chaque codon
  rscu_values <- numeric()
  codon_dict[c('Stop',"Trp","Met")] <- NULL
  for (aa in names(codon_dict)) {
    codon_list <- codon_dict[[aa]]
    usage <- sequence_codons_count[,codon_list]
    total_usage <- sum(usage, na.rm = TRUE)
    n_codons <- length(codon_list)
    
    # Calcul du RSCU
    rscu <- usage / (total_usage / n_codons)
    # Détecter les NaN et les remplacer par 0.0001
    rscu[is.nan(rscu)] <- 0.0001
    rscu[rscu==0] <- 0.0001
    rscu_values[codon_list] <- rscu
  }
  return(rscu_values)
}

# Fonction principale pour appliquer le calcul à un DNAStringSet
calculate_rscu_for_set <- function(dna_set, codon_dict, codons_count) {
  # Créer une liste pour stocker les résultats
  rscu_matrix_list <- list()
  
  for (i in seq_along(dna_set)) {
    seq_name <- names(dna_set)[i]
    sequence_codons_count <- codons_count[[seq_name]]
    rscu_values <- calculate_rscu(dna_set[[i]], codon_dict, sequence_codons_count)
    
    if (!is.null(rscu_values)) {
      rscu_matrix_list[[seq_name]] <- rscu_values
    }
  }
  # Transformer la liste en matrice
  rscu_matrix <- do.call(cbind, rscu_matrix_list)
  rownames(rscu_matrix) <- names(rscu_values)
  
  # Calcul du RSCU moyen pour l'ensemble des séquences par codon
  return(rscu_matrix)
}

# RSCU Distribution

# save RSCU distribution  per gene in asked repo
display_rscu_distribution <- function(total_matrix, rscu_output_repo) {
  # Filtrer les colonnes nécessaires
  data_filtered <- total_matrix %>% 
    select(all_of(features_to_use), gene)
  
  # Restructurer en format long
  data_long <- data_filtered %>%
    pivot_longer(cols = all_of(features_to_use), names_to = "Codon", values_to = "RSCU")
  
  codons_dict_reduced <- codon_dict
  #Met, Trp, Stop
  codons_dict_reduced[["Met"]] <- NULL
  codons_dict_reduced[["Trp"]] <- NULL
  codons_dict_reduced[["Stop"]] <- NULL
  
  # Déterminer l'ordre des codons et leurs annotations (acides aminés)
  codon_order <- unlist(codons_dict_reduced)
  codon_labels <- rep(names(codons_dict_reduced), 
                      times = sapply(codons_dict_reduced, length))
  
  # Ajouter une colonne pour les annotations des acides aminés
  data_long <- data_long %>%
    mutate(
      AminoAcid = factor(codon_labels[match(Codon, codon_order)], 
                         levels = names(codon_dict)),
      Codon = factor(Codon, levels = codon_order) # Ordre des codons
    )
  
  # Déterminer les positions moyennes pour chaque groupe d'acides aminés
  amino_acid_positions <- data_long %>%
    group_by(AminoAcid) %>%
    summarise(
      Position = mean(as.numeric(Codon)),
      Start = min(as.numeric(Codon)) - 0.5,
      End = max(as.numeric(Codon)) + 0.5,
      .groups = "drop"
    )
  
  # Créer un graphe pour chaque gène
  plots <- data_long %>%
    group_by(gene) %>%
    group_split() %>%
    lapply(function(df) {
      ggplot(df, aes(x = Codon, y = RSCU, fill = AminoAcid)) +
        geom_boxplot(outlier.size = 1.5) +
        # Ajouter des annotations des acides aminés
        annotate("text", 
                 x = amino_acid_positions$Position, 
                 y = max(df$RSCU, na.rm = TRUE) + 0.5,
                 label = amino_acid_positions$AminoAcid, 
                 size = 4, fontface = "bold") +
        # Ajouter des lignes verticales entre les groupes
        geom_vline(xintercept = amino_acid_positions$End[-nrow(amino_acid_positions)], 
                   color = "black", linetype = "dashed") +
        labs(title = paste("Distribution des RSCU pour le gène", unique(df$gene)),
             x = "Codons (par acide aminé)", y = "RSCU") +
        theme_bw() +
        theme(
          legend.position = "none",            # Supprimer la légende
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotation des labels
          plot.title = element_text(hjust = 0.5)
        )
    })
  
  results <- list()
  # Sauvegarder les graphes pour chaque gène
  for (i in seq_along(plots)) {
    g <- plots[[i]]
    g_name <- unique(data_long$gene)[i]
    ggsave(paste0(rscu_output_repo,"/RSCU_distribution_", g_name, ".png"), 
           g, width = 12, height = 6)
    results[[g_name]] <- g
  }
  
  return(results)
}

# Fonction pour calculer le RCDI
calculate_rcdi <- function(codon_counts, codon_weights) {
  
  # Vérification des noms des codons
  if (!all(names(codon_counts) %in% names(codon_weights))) {
    stop("Tous les codons présents dans 'codon_counts' doivent exister dans 'codon_weights'.")
  }
  
  # Extraire uniquement les codons avec un compte non nul
  codon_counts_filtered <- codon_counts[,codon_counts > 0]
  
  # Normaliser les poids pour ne considérer que les codons présents
  relevant_weights <- codon_weights[names(codon_counts_filtered)]
  
  # Calculer la fréquence relative des codons
  total_codons <- sum(codon_counts_filtered)
  codon_frequencies <- codon_counts_filtered / total_codons
  
  # Calculer le RCDI
  rcdi <- prod(relevant_weights ^ codon_frequencies)
  return(rcdi)
}

calculate_dna_set_rcdi <- function(gene_codons_freq_list, weights){
  rcdi_matrix_list <- list()
  
  for (i in seq_along(gene_codons_freq_list)){
    seq_name <- names(gene_codons_freq_list)[i]
    rcdi_matrix_list[[seq_name]] <- calculate_rcdi(gene_codons_freq_list[[i]], weights)
  }
  # Transformer la liste en matrice
  rcdi_matrix <- do.call(cbind, rcdi_matrix_list)
  rownames(rcdi_matrix) <- names("rcdi")
  return(rcdi_matrix)
}

# Calcul pour chacune des séquences d'un DNAStringSet de gc1, gc2, gc3 et gc12
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
  # Transformer la liste en matrice
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

# Calcul des fréquences attendues des dinucléotides
calculate_expected_dinucleotides_freq <- function(data, nucleotides, dinucleotides) {
  expected <- data[dinucleotides]  # Copie pour structure identique
  for (dinuc in dinucleotides) {
    nuc1 <- substr(dinuc, 1, 1)  # Premier nucléotide
    nuc2 <- substr(dinuc, 2, 2)  # Second nucléotide
    expected[[dinuc]] <- data[[nuc1]] * data[[nuc2]]  # Produit des fréquences
  }
  return(expected)
}

# Calcul du ratio observé/attendu
calculate_ratio_dinucleotides_freq <- function(observed, expected) {
  ratio <- observed / expected
  return(ratio)
}

# Calculate expected dinuc freq dataframe and ratio observed/expected
# dinuc freq as dataframes
calculate_exp_ratio_dinuc_freq <- function(total_data, nucleotides, dinucleotides){
  exp_ratio_dinuc_freq_list <- list()
  data <- total_data[,c(nucleotides, dinucleotides)]
  expected_dinuc_df<- calculate_expected_dinucleotides_freq(
    data, nucleotides, dinucleotides
  )
  exp_ratio_dinuc_freq_list[["expected_dinuc_df"]] <- expected_dinuc_df
  exp_ratio_dinuc_freq_list[["ratio_dinuc_df"]] <- calculate_ratio_dinucleotides_freq(
    data[dinucleotides], expected_dinuc_df
  )
  return(exp_ratio_dinuc_freq_list)
}

### Get RSCU compaison with reference (human here)
# Ajout de la colonne RSCU de "human_codon_usage_df" à notre df
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

# compute features mean per gene as dataframe and can save it
compute_features_mean_per_gene <- function(
    total_matrix,
    save=FALSE,
    file_path=NULL
){
  # Group by sum of multiple columns
  df_mean_features <- total_matrix %>% group_by(gene) %>% 
    summarise(gc1=mean(gc1),
              gc2= mean(gc2),
              gc3=mean(gc3),
              gc12=mean(gc12),
              rcdi=mean(rcdi),
              enc=mean(enc),
              cai=mean(cai),
              gc=mean(gc),
              aromatic_prct=mean(aromatic_prct),
              gravy=mean(gravy),
              .groups = 'drop') %>%
    as.data.frame() %>% data.frame(row.names = 1) 
  
  # Calculer les moyennes pour chaque colonne
  #all_genes <- colMeans(df_mean_features, na.rm = TRUE)
  
  # Ajouter la nouvelle ligne au dataframe
  #df_mean_features["all_genes", ] <- all_genes
  
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

gene_sequences_features_computation <- function(
    gene_sequences,
    codon_dict,
    gene,
    ref_codon_weights) {
  
  # get codon frequency
  gene_sequences_cf <- cubar::count_codons(gene_sequences, as.prob = TRUE)
  
  # get codon table for the standard genetic code
  ctab <- get_codon_table(gcid = '1')
  
  # get ENC
  codon_data <- codonTable(gene_sequences)
  enc<- ENC(codon_data)
  enc_matrix <- matrix(unlist(enc), ncol = length(enc), byrow = TRUE)
  rownames(enc_matrix) <- c("enc")
  colnames(enc_matrix) <- names(gene_sequences)
  enc_df <- as.data.frame(enc_matrix)
  
  # get GC content
  gc <- get_gc(gene_sequences_cf)
  gc_matrix <- matrix(unlist(gc), ncol = length(gc), byrow = TRUE)
  rownames(gc_matrix) <- c("gc")
  colnames(gc_matrix) <- names(gene_sequences)
  gc_df <- as.data.frame(gc_matrix)
  
  gc3s <- get_gc3s(gene_sequences_cf)
  
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
  rscu_results <- calculate_rscu_for_set(gene_sequences, codon_dict, codons_count)
  
  #compute RCDI
  rcdi_results <- calculate_dna_set_rcdi(codons_count, ref_codon_weights)
  rownames(rcdi_results) <- c("rcdi")
  
  #compute mean gene RSCU
  # Calculer la moyenne des lignes
  mean_rscu <- rowMeans(rscu_results)
  # Créer un dataframe avec les résultats
  mean_rscu_df <- data.frame(
    mean_rscu
  )
  colnames(mean_rscu_df) <- gene
  
  #compute GC3 and GC12
  GC12_GC3_results <- calc_GC12_GC3(gene_sequences)
  gc1_matrix <- GC12_GC3_results$gc1_matrix
  gc2_matrix <- GC12_GC3_results$gc2_matrix
  gc3_matrix <- GC12_GC3_results$gc3_matrix
  gc12_matrix <- GC12_GC3_results$gc12_matrix
  
  # replace gaps with "N"
  gene_sequences_with_N <- gsub("-", "N", as.character(gene_sequences))
  gene_sequences_with_N <- DNAStringSet(gene_sequences_with_N)
  
  # Traduction
  gene_aa_sequences <- Biostrings::translate(
    gene_sequences_with_N, 
    no.init.codon=FALSE,
    if.fuzzy.codon="X")
  
  # GRAVY score computation for every translated sequence
  gravy <- sapply(
    as.character(gene_aa_sequences), Peptides::hydrophobicity
  )
  gravy_matrix <- matrix(unlist(gravy), ncol = length(gravy), byrow = TRUE)
  rownames(gravy_matrix) <- c("gravy")
  colnames(gravy_matrix) <- names(gene_sequences)
  gravy_df <- as.data.frame(gravy_matrix)  
  
  # Calcul du pourcentage de résidus aromatiques pour chaque séquence
  # aromatic residues percentage computation for every translated sequence
  aromatic_prcts <- list()
  for (i in seq_along(gene_aa_sequences)) {
    seq_name <- names(gene_aa_sequences)[i]
    arom_tab <- aaComp(as.character(gene_aa_sequences[i]))
    aromatic_prct <- arom_tab[[seq_name]]["Aromatic","Mole%"]
    aromatic_prcts[[seq_name]] <- aromatic_prct
  }
  # Transformer la liste en matrice
  aromatic_matrix <- do.call(cbind, aromatic_prcts)
  rownames(aromatic_matrix) <- c("aromatic_prct")
  aromatic_matrix <- aromatic_matrix
  
  # annotate seq with gene name
  gene_list <- (rep(gene, length(gene_sequences)))
  gene_matrix <- matrix(unlist(gene_list), ncol = length(gene_sequences), byrow = TRUE)
  rownames(gene_matrix) <- c("gene")
  colnames(gene_matrix) <- names(gene_sequences)
  gene_df <- as.data.frame(gene_matrix)
  
  # all gene sequences features matrix creation
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
  
  gene_results_list <- list()
  gene_results_list[["total_gene_df"]] <- total_gene_df
  gene_results_list[["mean_rscu_df"]] <- mean_rscu_df
  
  return(gene_results_list)
}


all_genes_matrix_computation <- function(
    all_gene_sequences,
    codon_dict,
    ref_codon_weights) {
  all_genes_list <- list()
  all_genes_rscu_mean_list <- list()
  
  # computation of features matrix for each gene sequences
  for (gene in names(all_gene_sequences)) {
    gene_sequences <- all_gene_sequences[[gene]]
    
    gene_results_list <- gene_sequences_features_computation(
      gene_sequences,
      codon_dict,
      gene,
      ref_codon_weights)
    
    
    all_genes_list[[gene]] <- gene_results_list$total_gene_df
    all_genes_rscu_mean_list[[gene]] <- gene_results_list$mean_rscu_df
    
  }
  
  # Vérifier que tous les noms de lignes sont identiques
  row_names <- lapply(all_genes_list, rownames)
  if (!all(sapply(row_names, function(x) identical(x, row_names[[1]])))) {
    stop("Les noms des lignes ne correspondent pas entre les dataframes.")
  }
  
  # Concaténer par colonne
  all_genes_df <- do.call(cbind, all_genes_list)
  all_genes_rscu_mean_df <- do.call(cbind, all_genes_rscu_mean_list)
  
  # transpose df
  translated_matrix_1 <- as.data.frame(t(all_genes_df))
  translated_matrix <- subset(translated_matrix_1, select = -c(gene, other))
  
  #transform variables from character to numeric
  translated_matrix <- cbind(
    mutate_all(
      translated_matrix, function(x) as.numeric(as.character(x))
    ), translated_matrix_1[,"gene", drop=FALSE]
  )
  
  all_genes_results_list <- list()
  all_genes_results_list[["total_matrix"]] <- translated_matrix
  all_genes_results_list[["all_genes_rscu_mean_df"]] <- all_genes_rscu_mean_df
  return(all_genes_results_list)
}

# computes Correspondance Analysis, produces plots and matrix, save
# and return them.
compute_CA <- function(total_matrix, ca_output_repo){
  results <- list()
  
  all_features <- c(features_to_use, quanti_supp_features)
  
  # Sélection des colonnes d'intérêt
  wanted_data <- total_matrix[, features_to_use]
  head(wanted_data)
  
  # Suppression des lignes avec des valeurs infinies ou NaN (juste pour vérifier)
  wanted_data <- wanted_data[complete.cases(wanted_data), ]
  
  # Analyse de Correspondance (CA)
  res_cca <- cca(wanted_data)
  
  # Récupérer les scores des points (sites) et des variables (species)
  points_scores <- as.data.frame(scores(res_cca, display = "sites"))
  variables_scores <- as.data.frame(scores(res_cca, display = "species"))
  
  # Nettoyer les noms de lignes
  rownames(total_matrix) <- trimws(as.character(rownames(total_matrix)))
  rownames(points_scores) <- trimws(as.character(rownames(points_scores)))
  
  # Ajout de la colonne "gene"
  points_scores$gene <- total_matrix[rownames(points_scores), "gene"]
  
  # Définir le fichier de sortie SVG
  svg(paste0(ca_output_repo, "ca_clustered.svg"), width = 7, height = 5)  # Ajuster les dimensions
  
  # Ajuster les marges pour laisser plus de place à la légende
  par(mar = c(2.5, 2.5, 5, 8))  # (bottom, left, top, right)
  
  # Tracer le plot de l'Analyse de Correspondance
  plot(res_cca, type = "n", main = "Correspondance Analysis")
  
  # Définir les couleurs par gène
  unique_genes <- unique(points_scores$gene)
  colors <- rainbow(length(unique_genes))
  
  # Ajouter les points pour chaque gène
  for (i in seq_along(unique_genes)) {
    gene <- unique_genes[i]
    points_subset <- subset(points_scores, gene == unique_genes[i])
    points(points_subset[, 1], points_subset[, 2], 
           col = colors[i], pch = 19, cex = 0.8)
  }
  
  # Ajouter une légende à droite
  legend("topright", inset = c(-0.3, 0),  # Ajuster la position de la légende
         legend = unique_genes, col = colors, 
         pch = 19, cex = 1, #cex = 0.7             # Réduire la taille de la légende
         xpd = TRUE)                      # Permet à la légende de dépasser la zone graphique
  
  # Fermer l'appareil graphique
  dev.off()
  
  # Extraire les scores des variables
  var_scores <- as.data.frame(scores(res_cca, display = "species"))
  
  # Calculer la contribution des variables pour les deux premières dimensions
  contrib_dim1 <- abs(var_scores[, 1]) / sum(abs(var_scores[, 1])) * 100
  contrib_dim2 <- abs(var_scores[, 2]) / sum(abs(var_scores[, 2])) * 100
  
  # Créer un dataframe avec les contributions
  contrib_df <- data.frame(
    Variable = rownames(var_scores),
    Contribution_Dim1 = contrib_dim1,
    Contribution_Dim2 = contrib_dim2
  )
  
  # Sélectionner les 10 variables les plus contributrices pour chaque dimension
  top10_dim1 <- contrib_df[order(contrib_df$Contribution_Dim1, decreasing = TRUE), ][1:10, ]
  top10_dim2 <- contrib_df[order(contrib_df$Contribution_Dim2, decreasing = TRUE), ][1:10, ]
  
  # Diagramme en barres pour la dimension 2
  plot_dim1 <- ggplot(top10_dim1, aes(x = reorder(Variable, Contribution_Dim1), y = Contribution_Dim1)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +
    labs(title = "Top 10 Variables Contribuant à la Dimension 1", y = "Contribution (%)", x = "Variable")
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(ca_output_repo, "top10_dim1.svg"), plot = plot_dim1, width = 10, height = 7)
  results[["plot_dim1"]] <- plot_dim1
  
  # Diagramme en barres pour la dimension 2
  plot_dim2 <- ggplot(top10_dim2, aes(x = reorder(Variable, Contribution_Dim2), y = Contribution_Dim2)) +
    geom_bar(stat = "identity", fill = "lightgreen") +
    coord_flip() +
    labs(title = "Top 10 Variables Contribuant à la Dimension 2", y = "Contribution (%)", x = "Variable")
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(ca_output_repo, "top10_dim2.svg"), plot = plot_dim2, width = 10, height = 7)
  results[["plot_dim2"]] <- plot_dim2
  
  # Extraire les valeurs propres et calculer la variance expliquée
  eigenvalues <- res_cca$CA$eig
  var_explained <- eigenvalues / sum(eigenvalues) * 100
  
  # Créer un dataframe pour les 10 premières dimensions
  var_df <- data.frame(
    Dimension = c("1","2","3","4","5","6","7","8","9","10"),
    Variance_Explained = var_explained[1:10]
  )
  
  # add percentage to display
  var_df <- var_df |> 
    dplyr::mutate(perc = paste0(round(Variance_Explained, 1), "%"))
  
  # Diagramme en barres de la variance expliquée
  variance_explained_plot <- ggplot(var_df, aes(x = Dimension, y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "orange") +
    scale_x_discrete(limits=var_df$Dimension)+
    ## add percentage labels
    geom_text(aes(label = perc), 
              hjust = 0.5, nudge_y = +.7) + 
    labs(title = "Variance Expliquée par les 10 Premières Dimensions", y = "Variance Expliquée (%)", x = "Dimension") 
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(ca_output_repo, "variance_explained.svg"), plot = variance_explained_plot, width = 10, height = 7)
  results[["variance_explained_plot"]] <- variance_explained_plot
  
  ### Corrélations entre features
  
  # Extraire les scores des observations pour les deux premières dimensions
  obs_scores <- as.data.frame(scores(res_cca, display = "sites"))
  
  # Extraire les variables numériques de total_matrix (en excluant 'gene')
  numeric_vars <- total_matrix[, sapply(total_matrix, is.numeric)]
  
  # Initialiser une liste pour stocker les résultats
  cor_results <- list()
  
  # Fonction pour calculer la corrélation de Pearson et la p-value
  get_cor_and_pval <- function(var, axis_scores) {
    result <- cor.test(var, axis_scores, method = "pearson")
    list(correlation = result$estimate, p_value = result$p.value)
  }
  
  # Boucle pour calculer la corrélation pour chaque variable et chaque axe
  for (var_name in colnames(numeric_vars)) {
    var_data <- numeric_vars[[var_name]]
    
    # Corrélation avec Dim1
    cor_dim1 <- get_cor_and_pval(var_data, obs_scores[, 1])
    
    # Corrélation avec Dim2
    cor_dim2 <- get_cor_and_pval(var_data, obs_scores[, 2])
    
    # Ajouter les résultats à la liste
    cor_results[[var_name]] <- c(
      Correlation_Dim1 = cor_dim1$correlation,
      P_Value_Dim1 = cor_dim1$p_value,
      Correlation_Dim2 = cor_dim2$correlation,
      P_Value_Dim2 = cor_dim2$p_value
    )
  }
  
  # Convertir les résultats en dataframe
  cor_df <- do.call(rbind, cor_results)
  cor_df <- as.data.frame(cor_df)
  
  # Fonction pour calculer la matrice de corrélation et les p-values entre features
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
  
  # Calcul de la corrélation et p-value entre features
  correlation_results <- calculate_correlation_between_features(numeric_vars)
  
  # Filtrer les paires avec p-value < 0.05
  significant_features_corr <- correlation_results %>%
    filter(p_value < 0.05)
  
  # Enregistrer les résultats significatifs dans un fichier CSV
  write.csv(significant_features_corr, 
            paste0(ca_output_repo,"significant_features_correlations.csv"), 
            row.names = FALSE)
  
  # Afficher le tableau final
  write.csv(cor_df,paste0(ca_output_repo,"correlation_features_dimensions.csv"))
  results[["correlation_features_dimensions"]] <- cor_df
  results[["numeric_vars"]] <- numeric_vars # ADDED FOR TEST
  
  return(results)
}

# produces and save ENC plot
get_enc_plot <- function(total_matrix, enc_plot_repo){
  # Tracer ENC (Effective Number of Codons) vs GC3
  enc_curv <- function(x) { 2 + x + 29 / (x^2 + (1 - x)^2) }
  straigth_line <- function(x) {x}
  
  enc_plot <- ggplot(total_matrix, aes(x = gc3, y = enc, color = gene)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(stat='function', fun=enc_curv, color='green', linetype='dashed', linewidth=1) +
    xlim(0, 1) +
    ylim(30, 70)+
    labs(title = "ENC Plot",
         x = "GC3 Content",
         y = "Effective Number of Codons (ENC)",
         color = "Gene") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # Fond blanc
      plot.background = element_rect(fill = "white", color = NA),   # Fond global blanc
      legend.position = "right"
    )
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(enc_plot_repo, "enc_plot.svg"), plot = enc_plot, width = 7, height = 5)
  
  return(enc_plot)
}

# produces and save neutrality plot
get_neutrality_plot <- function(total_matrix, neutrality_plot_repo){
  # Tracer GC12 vs GC3
  neutrality_plot <- ggplot(total_matrix, aes(x = gc3, y = gc12, color = gene)) +
    geom_point(size = 3, alpha = 0.8)+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "green")+
    #geom_line(stat='function', fun=straigth_line, color='green', linetype='dashed', linewidth=1) + #
    labs(title = "Neutrality Plot",
         x = "GC3 Content",
         y = "GC12 Content",
         color = "Gene") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # Fond blanc
      plot.background = element_rect(fill = "white", color = NA),   # Fond global blanc
      legend.position = "right"
    )
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(neutrality_plot_repo, "neutrality_plot.svg"), plot = neutrality_plot, width = 7, height = 5)
  
  return(neutrality_plot)
}


get_dinuc_ratio_plot <- function(ratio_dinuc_freq_df, outdir){
  # Ajouter une colonne pour les gènes (extraction avant le ".")
  ratio_dinuc_freq_df <- ratio_dinuc_freq_df %>%
    rownames_to_column(var = "id") %>%
    mutate(gene = sub("\\..*", "", id)) %>%
    select(-id)
  
  # Calculer la moyenne pour chaque dinucléotide par gène
  data_summary <- ratio_dinuc_freq_df %>%
    group_by(gene) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  
  # Convertir en format long pour ggplot
  data_long <- data_summary %>%
    pivot_longer(-gene, names_to = "dinucleotide", values_to = "oe_dinuc_ratio")
  
  
  # Créer le graphe à points
  dinuc_plot <- ggplot(data_long, aes(x = dinucleotide, y = oe_dinuc_ratio, color = gene, group = gene)) +
    geom_point(size = 3) +
    geom_line() +
    labs(
      title = "Observed/Expected dinucleotides frequency ratio per gene",
      x = "Dinucleotide",
      y = "O/E ratio",
      color = "gene"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),  # Fond blanc
      plot.background = element_rect(fill = "white", color = NA),   # Fond global blanc
      legend.position = "right"
    )
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(outdir, "oe_dinuc_freq_plot.svg"), plot = dinuc_plot, width = 7, height = 5)
}


get_gc_plot <- function(mean_df, output_repo){
  gc_features <- c("gc1", "gc2", "gc3", "gc12", "gc")
  reduced_mean_df <- mean_df[gc_features, ]
  # Convertir la matrice en un dataframe long
  data_long <- data.frame(
    type_gc = rep(rownames(reduced_mean_df), times = ncol(reduced_mean_df)),
    gène = rep(colnames(reduced_mean_df), each = nrow(reduced_mean_df)),
    valeur_gc = as.vector(reduced_mean_df)
  )
  
  # Créer un diagramme en barres groupées
  gc_plot <- # Créer le diagramme en barres groupé
    ggplot(data_long, aes(x = gène, y = valeur_gc, fill = type_gc)) +
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
  
  # Enregistrer le graphique au format SVG
  ggsave(paste0(output_repo, "gc_plot.svg"), plot = gc_plot, width = 7, height = 5)
}

remove_bad_sequences <- function(gene_sequences){
  # Seuil de tolérance pour les "n"
  threshold <- 0.15  # 15%
  
  # Fonction pour calculer le pourcentage de "n" dans une séquence
  calc_n_percentage <- function(seq) {
    # n_count <- sum(as.character(seq) == "n" | as.character(seq) == "N")  # Nombre de "n"
    # n_percentage <- n_count / length(seq)  # Pourcentage 
    n_frequency <-  1 -letterFrequency(seq, "ACGT", as.prob=TRUE)
    return(n_frequency)
  }
  
  # Filtrer les séquences avec moins de 15% de "n"
  filtered_sequences <- gene_sequences[sapply(gene_sequences, calc_n_percentage) <= threshold]
  
  return(filtered_sequences)
}

################################
######## ---- Main ---- ########
################################

# ---- Reference codon weights computation ----

ref_codon_weights <- NULL

ref_codon_usage_df <- validate_csv(codon_usage_ref_file, verbose=FALSE)
ref_codon_weights <- get_ref_codon_weights(ref_codon_usage_df)

# ---- Bad sequences removal ----

genes_good_sequences_list <- lapply(genes_sequences_list, remove_bad_sequences)

# ---- Complete genomes selection ----

# Identifier les génomes ayant une séquence pour chaque gène
all_genomes <- lapply(genes_good_sequences_list, names)  # Extraire les noms des séquences pour chaque gène
common_genomes <- Reduce(intersect, all_genomes)  # Trouver les génomes communs à tous les gènes

# Filtrer les séquences pour ne conserver que les génomes communs
filtered_gene_sets <- lapply(genes_good_sequences_list, function(gene_set) {
  gene_set[common_genomes]
})

# Concaténer les séquences de chaque gène pour les génomes sélectionnés
complete_genomes <- DNAStringSet(sapply(common_genomes, function(genome) {
  paste0(sapply(filtered_gene_sets, 
                function(gene_set) as.character(gene_set[genome])), 
         collapse = "")
}))

# Attribuer les noms des génomes aux séquences concaténées
names(complete_genomes) <- common_genomes

genes_good_sequences_list[["complete"]] <- complete_genomes

# ---- RSCU and features value computation ----

all_genes_results_list <- all_genes_matrix_computation(
  genes_good_sequences_list,
  codon_dict,
  ref_codon_weights)

total_matrix <- all_genes_results_list$total_matrix
all_genes_rscu_mean_df <- all_genes_results_list$all_genes_rscu_mean_df

# save raw data
write.csv(total_matrix, file.path(outdir, "total_matrix.csv"))

# ---- Dinucleotide ratio ----

# compute expected dinuc freq and ratio obs/exp dinuc freq
exp_ratio_dinuc_list <- calculate_exp_ratio_dinuc_freq(
  total_matrix, nucleotides, dinucleotides
)

exp_dinuc_freq_df <- exp_ratio_dinuc_list$expected_dinuc_df
ratio_dinuc_freq_df <- exp_ratio_dinuc_list$ratio_dinuc_df
write.csv(ratio_dinuc_freq_df, file.path(outdir, "ratio_dinuc_freq_df.csv"))

get_dinuc_ratio_plot(ratio_dinuc_freq_df, outdir)

# ---- Mean features value per gene ----

# get mean features value per gene
mean_df <- compute_features_mean_per_gene(total_matrix)
write.csv(mean_df, file.path(outdir, "mean_df.csv"))

# ---- RSCU comparison with ref ----

comparison_rscu_df <- create_RSCU_comparison_file(all_genes_rscu_mean_df, 
                                                  ref_codon_usage_df)
write.csv(comparison_rscu_df, file.path(outdir, "rscu_df.csv"))

# ---- RSCU distribution ----

# generate RSCU distribution charts
rscu_distribution_plots <- display_rscu_distribution(total_matrix, outdir) 

# ---- CA computation ----

# generate Comparison Analysis matrix and charts
ca_results <- compute_CA(total_matrix, outdir)

# ---- GC Plot -----

get_gc_plot(mean_df, outdir)

# ---- ENC Plot ----

enc_plot <- get_enc_plot(total_matrix, outdir)

# ---- Neutrality Plot ----

neutrality_plot <- get_neutrality_plot(total_matrix, outdir)