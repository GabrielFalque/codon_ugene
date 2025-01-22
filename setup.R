#!/usr/bin/env Rscript

# Setup script to prepare the environment for codon_ugene.R

cat("Checking required R packages...\n")

# List of required libraries
required_packages <- c(
  "cubar", "Biostrings", "ggplot2", "coRdon", "tidyverse", "seqinr",
  "Peptides", "vegan", "dplyr", "optparse", "this.path", "tictoc", "ggtext", 
  "parallel", "future", "future.apply"
)

# Install missing libraries
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    cat("Installing missing packages:", new_packages, "\n")
    install.packages(new_packages)
  }
}

install_if_missing(required_packages)
cat("Setup complete! All required packages are installed.\n")