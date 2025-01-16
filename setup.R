#!/usr/bin/env Rscript

# Setup script to prepare the environment for codon_ugene.R

cat("Checking required R packages...\n")

# Liste des librairies nécessaires
required_packages <- c(
  "cubar", "Biostrings", "ggplot2", "coRdon", "seqinr", "stringr", "tidyverse",
  "ca", "Peptides", "FactoMineR", "factoextra", "vegan", "reshape2",
  "AnaCoDa", "dplyr", "svglite", "optparse", "here", "this.path"
)

# Installer les packages manquants
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) {
    cat("Installing missing packages:", new_packages, "\n")
    install.packages(new_packages)
  }
}

install_if_missing(required_packages)
cat("Setup complete! All required packages are installed.\n")
