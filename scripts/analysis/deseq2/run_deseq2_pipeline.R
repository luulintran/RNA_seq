# RUN ALL SCRIPTS IN ORDER: ----------------------------------------------------
# These will run all lines of code in each script listed below:

# MAKE SURE YOU EDIT THIS CONFIG SCRIPT WITH THE CORRECT DATA AND VARIABLES FOR YOUR EXPERIMENT
source("scripts/analysis/deseq2/config_deseq2.R")
print("Ran config script and defined variables.")
# ******************************************************************************

# 1. Run DESeq2
source("scripts/analysis/deseq2/01_run_deseq2.R")
print("Ran DESeq2 successfully. Now retrieving results.")

# 2. Retrieve DESeq2 results
source("scripts/analysis/deseq2/02_deseq2_results.R")
print("Retrieved DESeq2 results successfully. Now annotating results with gene names and Entrez ID's.")

# 3. Annotate DESeq2 results
source("scripts/analysis/deseq2/03_deseq2_results_annotated.R")
print("Annotated DESeq2 results successfully. Now filtering for genes based on gene list in script 04.")

# 4. Filtered DESeq2 results based on gene list
source("scripts/analysis/deseq2/04_filter_genes.R")
print("Filtered DESeq2 results successfully for genes based on gene list. Now making PCA plot.")

# 5. Make PCA plot
source("scripts/analysis/deseq2/05_pcaplot.R")
print("Made PCA plot successfully. Now making volcano plots.")

# 6. Make Volcano plots
source("scripts/analysis/deseq2/06_selected_degs_volcanoplot.R")
print("Made Volcano plots successfully. Now making heatmaps.")

# 7. Make Heatmaps
source("scripts/analysis/deseq2/07_selected_degs_heatmap.R")
print("Made heatmaps successfully.")

print("DESeq2 downstream analysis pipeline done.")
