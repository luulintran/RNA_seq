# setup_environment.R

# REFERENCE FILE STRUCTURE: ----------------------------------------------------
#project_dir/
#  ├── data/
#  │   ├── meta_data/
#  │   ├── processed_data/
#  ├── logs/
#  ├── notebooks/
#  │   ├── exploratory_analysis/
#  │   ├── final_reports/
#  ├── results/
#  │   ├── figures/
#  │   ├── r_objects/
#  │   ├── reports/
#  │   ├── tables/
#  ├── scripts/
#  │   ├── analysis/
#  │   ├── preprocessing/
#  │   ├── setup/

# ------------------------------------------------------------------------------
# Load renv for package management
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# Restore the environment from the renv.lock file
if (file.exists("renv.lock")) {
  cat("\nRestoring environment from renv.lock...\n")
  renv::restore()
} else {
  stop("renv.lock file not found. Check that it is in the project directory.")
}

# Load libraries
library(clusterProfiler)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(readr)
library(org.Mm.eg.db)
library(ggplot2)
library(extrafont)
library(EnhancedVolcano)

# import Arial font
font_import(pattern = "Arial", prompt = FALSE)

# Set up directory structure
required_directories <- c(
  "data",
  "data/meta_data",
  "data/processed_data",
  "data/processed_data/raw_counts",
  "notebooks",
  "notebooks/exploratory_analysis",
  "notebooks/final_reports",
  "results",
  "results/figures",
  "results/tables",
  "results/r_objects",
  "results/reports",
  "logs",
  "scripts",
  "scripts/analysis",
  "scripts/preprocessing",
  "scripts/setup"
)

for (dir in required_directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(paste0("Created directory: ", dir, "\n"))
  }
}

cat("\nEnvironment setup complete.\n")