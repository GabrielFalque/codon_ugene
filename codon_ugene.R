#!/usr/bin/env Rscript
options(warn=1)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

if (!requireNamespace("this.path", quietly = TRUE)) {
  invisible(install.packages("this.path", dependencies = TRUE)) #added invisible()
}
suppressWarnings(suppressMessages(library("this.path", character.only = TRUE)))

cat("Script directory : ", this.dir(), "\n")
data_directory = paste0(this.dir(), "/data/")
scripts_directory <- paste0(this.dir(), "/R")
  
source(paste0(this.dir(), "/dependencies.R"), local = FALSE, echo = FALSE)

source_without_messages <- function(
    script_path
){
  source(script_path, , local = FALSE, echo = FALSE)
}

# Dynamically load all scripts from R/ directory
script_files <- list.files(scripts_directory, pattern = "\\.R$", full.names = TRUE)
invisible(sapply(script_files, source_without_messages))

# compute elapsed time
tic()

################################
###### ---- Parsing ---- ######
################################

# Define command line options
option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help = paste0("Directory with gene sequences (mandatory).\n", 
              "Files must be named like gene_{gene_name}.fasta with gene_name", 
              " being Orf1 for example. (Ex : gene_Orf1.fasta)"), 
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
  
  make_option(c("--no_complete"), action = "store_true", default = FALSE,
              help = "If you want to remove complete genome CDS computation."),
  
  make_option(c("--cov"), action = "store_true", default = FALSE,
              help = "If the sequences are from SARS-CoV-2."),
  
  make_option(c("-t", "--threads"), type= "integer",
              default = 1,
              help = "Number of threads.",
              metavar = "INTEGER"),
  
  make_option(c("--verbose"), action = "store_true", default = FALSE,
              help = "Verbose mode.")
)

# Options parser
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
outdir <- opt$outdir

# Display help if no argument or if -h/--help is used
if (is.null(opt$dir) || is.null(outdir) || opt$help) {
  print_help(opt_parser)
  quit(status = 0)
}

# Verify if input directoy exists
if (!dir.exists(opt$dir)) {
  stop("Input directory does not exist : ", opt$dir, "\n")
}

# Create output directory if it does not exist
if (!dir.exists(outdir)) {
  if (opt$verbose) message("Output directory creation : ", outdir,"\n")
  dir.create(outdir, recursive = TRUE)
}

# Get all FASTA files from genes (with gene_*.fasta pattern)
fasta_files <- list.files(path = opt$dir, pattern = "^gene_.*\\.fasta$", full.names = TRUE)
if (length(fasta_files) == 0) {
  stop("No file corresponding to 'gene_*.fasta' pattern in the directory : ", opt$dir)
}

n_threads <- opt$threads
# Verify if number of threads is an integer
if (!n_threads %% 1 == 0){
  stop("-t/--threads (number of threads) option must be an integer")
} else {
  if (opt$verbose) message("Number of threads : ", n_threads,"\n")
}
n_threads <- min(n_threads, detectCores() - 1)



# Indicate if --no_complete option is activated
if (opt$no_complete) {
  message("Option --no_complete activated. No complete genome CDS codon analysis will be performed.\n")
} else {
  if (opt$verbose){
    message("Option --no_complete non activated. Complete CDS codon analysis will be performed (by default).\n")
  }
}

# Indicate if --cov option is activated
if (opt$cov) {
  message("Option --cov activated. Sequences are from SARS-CoV-2.\n")
} else {
  message("Option --cov non activated. Sequences are not from SARS-CoV-2.\n")
}

# Load reference file
codon_usage_ref_file <- opt$ref
if (!file.exists(codon_usage_ref_file)) {
  stop("Specified species codon usage file does not exist : ", codon_usage_ref_file)
}
if (opt$verbose) message("Loading species codon usage file : ", codon_usage_ref_file)

# Charger les séquences dans une liste
# Load sequences in a list
genes_sequences_list <- list()
for (file in fasta_files) {
  
  # Extract genes name
  gene_name <- sub("^gene_(.*)\\.fasta$", "\\1", basename(file))
  
  if (opt$verbose) message("Reading gene : ", gene_name)
    
    if (opt$cov && grepl("ORF1a", gene_name, fixed = TRUE)){ 
      
      raw_ORF1a_sequences <- readDNAStringSet(file)
      ORF1a_sequences <- raw_ORF1a_sequences[
        names(raw_ORF1a_sequences) != "NC_045512.2"]
      
      # Remove gaps if input file is an Multiple Alignement Sequences (MSA)
      nogap_ORF1a_sequences<- DNAStringSet(sapply(
        ORF1a_sequences, function(x) gsub("-", "", as.character(x))
      ))
      
      # Does not check for STOP codons 
      ORF1a_gene_name <- gene_name
      sequence_temp <- check_cds(
        nogap_ORF1a_sequences,
        check_len = TRUE,
        check_start = TRUE,
        check_stop = FALSE,
        check_istop = TRUE,
        rm_start = TRUE,
        rm_stop = FALSE,
      )
      
      genes_sequences_list[[ORF1a_gene_name]] <- sequence_temp
      
    } else if (opt$cov && grepl("ORF1b", gene_name, fixed = TRUE)) {
      
      # separate ORF1ab gene into ORF1a and ORF1b
      raw_ORF1b_sequences <- readDNAStringSet(file)
      ORF1b_sequences <- raw_ORF1b_sequences[
        names(raw_ORF1b_sequences) != "NC_045512.2"]
      
      # Remove gaps if input file is an MSA
      nogap_ORF1b_sequences <- DNAStringSet(sapply(
        raw_ORF1b_sequences, function(x) gsub("-", "", as.character(x))
      ))
      
      ORF1b_gene_name <- gene_name
      # Does not check for START codons 
      genes_sequences_list[[ORF1b_gene_name]] <- check_cds(
        nogap_ORF1b_sequences,
        check_len = TRUE,
        check_start = FALSE,
        check_stop = TRUE,
        check_istop = TRUE,
        rm_start = FALSE,
        rm_stop = TRUE,
      )
    
  } else {
    # Get gene sequences from file as DNAStringSet 
    raw_gene_sequences <- readDNAStringSet(file)
    raw_gene_sequences <- raw_gene_sequences[names(raw_gene_sequences) != "NC_045512.2"]
    nogap_gene_sequences <- DNAStringSet(sapply(
      raw_gene_sequences, function(x) gsub("-", "", as.character(x))
    ))
    genes_sequences_list[[gene_name]] <- check_cds(nogap_gene_sequences)
  }
}

################################
######## ---- Main ---- ########
################################

# ---- Reference codon weights computation ----
# Verify reference codon usage csv file then
# computes codons weights in reference species genome

ref_codons_infos_list <- get_ref_codon_usage(
  codon_usage_ref_file,
  opt$verbose
)

ref_codon_usage_df <- ref_codons_infos_list$ref_codon_usage_df
ref_codon_weights <- ref_codons_infos_list$ref_codon_weights

# ---- Data cleaning ----
# Bad sequences (>15% ambiguious nucleotides) removal and complete genomes 
# selection.
cat("Data cleaning...")
all_gene_sequences <- data_cleaning(genes_sequences_list,
                                    opt$no_complete)
cat("Done !\n")

# ---- RSCU and features value computation ----
# Produces "total_matrix.csv" file in outdir with all features computed per
# sequence.
# Also produces "rscu_with_ref.csv" with RSCU values for each gene and for
# reference species genome.
# Following features for all sequences of all genes are present in total_matrix :
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

cat("Features computation...")

all_genes_results_list <- features_computation(
  all_gene_sequences = all_gene_sequences,
    CODON_DICT = CODON_DICT,
    ref_codon_weights = ref_codon_weights,
    ref_codon_usage_df = ref_codon_usage_df,
    n_threads = n_threads,
    scripts_directory = scripts_directory,
    outdir = outdir
)

total_matrix <- all_genes_results_list$total_matrix
all_genes_rscu_mean_df <- all_genes_results_list$all_genes_rscu_mean_df
comparison_rscu_df <- all_genes_results_list$comparison_rscu_df
cat("Done !\n")

# ---- Set color palette for plots ----
# Set unique color palette with a color assigned to each gene in order to keep
# an homogeneity through the plots

# Get unique gene names
unique_genes <- unique(total_matrix$gene)

if (!opt$no_complete){
  # Reorder to put ‘complete’ first
  gene_list <- unique_genes[order(unique_genes != "complete", unique_genes)]
} else {
  gene_list <- unique_genes
}

# Associate unique color to each gene
# Use color palette from ggplot2
palette <- scales::hue_pal()(length(gene_list))
gene_color_map <- setNames(palette, gene_list)

# ---- Dinucleotides relative abundance ----
# Produces "ratio_dinuc_freq_df.csv" file in outdir with Dinucleotides relative 
# abundance per sequence.

ratio_dinuc_freq_df <- dinucleotide_relative_abundance(
    total_matrix,
    NUCLEOTIDES,
    DINUCLEOTIDES,
    outdir
)

# ---- Features mean value per gene ----
# Produces "gene_features_mean.csv" file in outdir with mean values per gene for 
# the following features :
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

gene_features_mean <- compute_gene_features_mean(
  total_matrix,
  outdir
)

# ---- CA computation ----
# Computes Correspondence Analysis (CA) on RSCU values from all sequences.
# Creates chart of CA with points colored by gene so it is possible to see
# different clusters. It also plots variables contribution to the first 2 
# dimensions. And plots how variance is explained by the 10 first dimensions.
# Produces two .csv files :
# - "correlation_features_dimensions.csv" : correlation between every feature
#   and the two first CA dimensions along with p-value.
# - "significant_features_correlations.csv" : correlation between every pair of
#   features for which p-value < 0.05 (significant)

cat("Correspondence Analysis computation...")
ca_results <- compute_CA(total_matrix, gene_color_map, outdir)
cat("Done!\n")

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

plots_generation(
    total_matrix,
    gene_features_mean,
    ratio_dinuc_freq_df,
    gene_color_map,
    outdir
)

cat(paste0("Plots generated at ", normalizePath(dirname(outdir)), "/", outdir, "\n"))
cat("Job finished !\n")
toc()