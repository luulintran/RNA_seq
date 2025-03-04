library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(dplyr)
library(readr)


# READ IN DATA:------------------------------------------------------------------------------------------
up_list <- read.csv(file.path("results/tables/deseq", paste0(project, "_deseq2_results_upgenes.csv")))
down_list <- read.csv(file.path("results/tables/deseq", paste0(project, "_deseq2_results_downgenes.csv")))

# FUNCTIONAL ENRICHMENT ANALYSIS: -----------------------------------------------------------------------------------------
# With clusterProfiler

# Make a list of genes for upregulated and downregulated genes
#Extract the gene symbols into a list
up_genes <- up_list$symbol
down_genes <- down_list$symbol

# GO term enrichment analysis
GO_results_up <- enrichGO(gene = up_genes,
                            OrgDb = "org.Mm.eg.db", #annotation database
                            keyType = "SYMBOL", #gene id type
                            ont = "BP") #ontology: BP (biological process), MP (molecular function), CC (cellular component)

GO_results_down <- enrichGO(gene = down_genes,
                            OrgDb = "org.Mm.eg.db", #annotation database
                            keyType = "SYMBOL", #gene id type
                            ont = "BP") #ontology: BP (biological process), MF (molecular function), CC (cellular component)

# Save results as CSV files
filename <- file.path(output_dir_tables, "GO_results_upgenes.csv")
write.csv(GO_results_up, 
          file = filename)

filename <- file.path(output_dir_tables, "GO_results_downgenes.csv")
write.csv(GO_results_down,
          file = filename)

print(paste0("GO Results saved in ", output_dir_tables))
