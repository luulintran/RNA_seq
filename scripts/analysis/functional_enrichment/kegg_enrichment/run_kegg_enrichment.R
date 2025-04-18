library(clusterProfiler)

# FUNCTIONAL ENRICHMENT ANALYSIS: ----------------------------------------------
# READ IN DATA:-----------------------------------------------------------------
up_list <- read.csv(file.path("results/tables/deseq", 
                              paste0(project, "_deseq2_results_upgenes.csv")))

down_list <- read.csv(file.path("results/tables/deseq", 
                                paste0(project, "_deseq2_results_downgenes.csv")))

# FUNCTIONAL ENRICHMENT ANALYSIS: ----------------------------------------------
# With clusterProfiler

# Make a list of genes for upregulated and downregulated genes
#Extract the gene entrez id's into a list
up_genes <- up_list$entrez
down_genes <- down_list$entrez

# REACTOME PATHWAY ANALYSIS: ---------------------------------------------------

# run Reactome Pathway Analysis
# UP genes
pathway_up <- enrichKEGG(up_genes,
                            organism = "mouse")

# DOWN genes
pathway_down <- enrichKEGG(down_genes,
                              organism = "mouse")

# SAVE RESULTS: ----------------------------------------------------------------
# Write to csv files
filename <- "kegg_results_upgenes.csv"

write.csv(
  pathway_up,   # specify dataframe, not string
  file = file.path(output_dir_tables, filename)
)

filename <- "kegg_results_downgenes.csv"

write.csv(
  pathway_down,   # specify dataframe, not string
  file = file.path(output_dir_tables, filename)
)


print(paste0("KEGG enrichment results saved in ", output_dir_tables))

