# SET UP
library(clusterProfiler)
library(org.Mm.eg.db)
library(DESeq2)
library(dplyr)
library(ggplot2)

# LOAD RDS FILE OF DDS OBJECT: -------------------------------------------------
dds <- readRDS(
  file = file.path("results/r_objects/deseq", paste0(project, "_deseq2_dds.rds"))
)

# STORE DESEQ2 RESULTS: --------------------------------------------------------
res <- results(dds)

# ANNOTATE RESULTS WITH GENE SYMBOLS AND ENTREZ IDS: ---------------------------
ensembl_ids <- rownames(res)

# annotate with gene symbols using org.Mm.eg.db package
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ensembl_ids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
# annotate with entrez ids
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ensembl_ids,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# Re-order res based on padj
resOrdered <- res[order(res$padj),]

# Move symbol and entrez columns to the front
resOrdered <- resOrdered[, c("symbol", 
                             "entrez", 
                             setdiff(colnames(resOrdered), 
                                     c("symbol", "entrez")))]

deseq2_results <- as.data.frame(resOrdered)

# FILTER DESEQ2 RESULTS FOR SIGNIFICANT UPREGULATED GENES: ---------------------
genes_up_list <- deseq2_results %>%
  dplyr:: filter(log2FoldChange > 0 & padj < 0.05)

# This is a df of all the significant upregulated genes in order by padj
genes_up_list <- genes_up_list[order(genes_up_list$padj), ]

# Display number of upregulated and downregulated genes
cat("Number of upregulated genes:", nrow(genes_up_list))

# FILTER DESEQ2 RESULTS FOR SIGNIFICANT DOWNREGULATED GENES: -------------------
genes_down_list <- deseq2_results %>%
  dplyr:: filter(log2FoldChange < 0 & padj < 0.05)

# This is a df of all the significant upregulated genes in order by padj
genes_down_list <- genes_down_list[order(genes_down_list$padj), ]

# Display number of upregulated and downregulated genes
cat("Number of downregulated genes:", nrow(genes_down_list))

# GO ANALYSIS: -----------------------------------------------------------------
# Extract the up gene symbols into a list
up_genes <- genes_up_list$symbol

GO_results_up <- enrichGO(gene = up_genes,
                          OrgDb = "org.Mm.eg.db", #annotation database
                          keyType = "SYMBOL", #gene id type
                          ont = "BP") #ontology: BP (biological process)

# Save GO results for upregulated genes as a dataframe
GO_up_genes_df <- as.data.frame(GO_results_up)

# Extract the down gene symbols into a list
down_genes <- genes_down_list$symbol

GO_results_down <- enrichGO(gene = down_genes,
                            OrgDb = "org.Mm.eg.db", #annotation database
                            keyType = "SYMBOL", #gene id type
                            ont = "BP") #ontology: BP (biological process)

# Save GO results for upregulated genes as a dataframe
GO_down_genes_df <- as.data.frame(GO_results_down)

# PLOT GO RESULTS FOR GLIOGENIC TERMS: -----------------------------------------

# Define the GO terms of interest
GO_gliogenesis <- c('gliogenesis', 
                    'glial cell differentiation', 
                    'myelination', 
                    'glial cell development', 
                    'oligodendrocyte development')

# Filter the dataframe for the relevant GO terms
GO_up_genes_df_glia <- GO_up_genes_df %>% 
  dplyr::filter(Description %in% GO_gliogenesis)

# Calculate -log10(p.adjust) values
GO_up_genes_df_glia$log_p.adjust <- -log10(GO_up_genes_df_glia$p.adjust)

# Make the bar plot
GO_glia_barplot <- ggplot(GO_up_genes_df_glia, 
                           aes(x = reorder(Description, log_p.adjust),
                               y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = "#f48c67") +
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(title = "Upregulated GO Terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
filename = file.path(output_dir_figures, paste0("selected_upregulated_go_barplot.pdf"))
pdf(filename, width = 5, height = 3)
print(GO_glia_barplot)

dev.off()

print(paste0("Upregulated GO barplot saved in ", output_dir_figures))




# GO ANALYSIS: -----------------------------------------------------------------
# Define the GO terms of interest
GO_neurogenesis <- c('regulation of neurogenesis', 
                     'dendrite development', 
                     'positive regulation of neurogenesis', 
                     'regulation of synapse structure or activity', 
                     'regulation of neuron differentiation')

# Filter the dataframe for the relevant GO terms
GO_down_genes_df_neuron <- GO_down_genes_df %>% 
  dplyr::filter(Description %in% GO_neurogenesis)

# Calculate -log10(p.adjust) values
GO_down_genes_df_neuron$log_p.adjust <- -log10(GO_down_genes_df_neuron$p.adjust)

# Create the bar plot
GO_neuron_barplot <- ggplot(GO_down_genes_df_neuron, 
                            aes(x = reorder(Description, log_p.adjust), 
                                y = log_p.adjust)) + 
  geom_bar(stat = "identity", fill = "#ca8bc9") + 
  coord_flip() +  # make horizontal barplot
  labs(title = "Downregulated GO terms", 
       x = "GO Term", 
       y = "-log10(p.adjust)") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# SAVE PLOT: -------------------------------------------------------------------
filename = file.path(output_dir_figures, paste0("selected_downregulated_go_barplot.pdf"))
pdf(filename, width = 5, height = 3)
print(GO_neuron_barplot)

dev.off()

print(paste0("Downregulated GO barplot saved in ", output_dir_figures))