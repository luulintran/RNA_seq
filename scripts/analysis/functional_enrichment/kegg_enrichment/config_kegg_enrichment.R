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
  
# 1. DEFINE FILES AND PATHS: ---------------------------------------------------
project <- "rnaseq_e16"
analysis <- "kegg"

input_coldata <- "data/meta_data/coldata.csv"
input_cts <- "data/processed_data/normalized_counts/rnaseq_e16_normalized_counts.tsv"

output_dir_robj <- file.path("results", "r_objects", analysis)
output_dir_tables <- file.path("results", "tables", analysis)
output_dir_figures <- file.path("results", "figures", analysis)

# Make sure output directory exists, and if it doesn't, create one
if (!dir.exists(output_dir_robj)) {
  dir.create(output_dir_robj, recursive = TRUE)
}

if (!dir.exists(output_dir_tables)) {
  dir.create(output_dir_tables, recursive = TRUE)
}

if (!dir.exists(output_dir_figures)) {
  dir.create(output_dir_figures, recursive = TRUE)
}



# 2. TREATMENT GROUPS: --------------------------------------

# Define control and mutant groups
control_group <- "control"     
mutant_group <- "NICD"  
