if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required libraries
required_packages <- c(
  "cubar", "Biostrings", "ggplot2", "coRdon", "tidyverse", "seqinr",
  "Peptides", "vegan", "dplyr", "optparse", "this.path", "tictoc", "ggtext", 
  "parallel", "future", "future.apply"
)

cat("Checking libraries...\n")

# Function to install and load libraries
install_and_load <- function(packages) {
  
  # Get missing libraries
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  # Install missing libraries
  if (length(missing_packages) > 0) {
    cat("Installing missing packages...\n")
    invisible(install.packages(missing_packages, dependencies = TRUE))
  }
  
  # Load all required libraries
  suppressWarnings(suppressMessages(lapply(packages, library, character.only = TRUE)))
}

bioc_packages <- c("Biostrings", "cubar", "coRdon")
cran_packages <- setdiff(required_packages, bioc_packages)
options(warn = -1) 
# Install and load CRAN libraries
install_and_load(cran_packages)
options(warn = 1) 
# Install and load Bioconductor libraries
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("Installing missing Bioconductor packages...\n")
  BiocManager::install(missing_bioc)
}
install_and_load(bioc_packages)

cat("Libraries checked!\n")