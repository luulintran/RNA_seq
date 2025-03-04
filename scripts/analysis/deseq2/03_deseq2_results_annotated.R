# Run this after running 'analysis/01_deseq2_e16.R'

# SET UP
library(DESeq2)
library(org.Mm.eg.db)
library(tidyverse)
library(readr)
library(dplyr)



# LOAD RDS FILE OF DDS OBJECT: -------------------------------------------------
dds <- readRDS(
  file = file.path(output_dir_robj, paste0(project, "_deseq2_dds.rds"))
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

# Rename resOrdered to a simpler name as dataframe
out <- as.data.frame(resOrdered)

# Export CSV file with gene symbols and entrez ids
filename <- paste0(project, "_deseq2_results_annotated.csv")
write.csv(out, file = file.path(output_dir_tables, filename))


# UP AND DOWN REGULATED GENE LISTS: --------------------------------------------
# Separate up and downregulated significant genes
out_up <- out %>%
  dplyr:: filter(log2FoldChange > 0, padj < 0.05)

# Separate up and downregulated significant genes
out_down <- out %>%
  dplyr:: filter(log2FoldChange < 0, padj < 0.05)

# Export CSV file with gene symbols and entrez ids
filename <- paste0(project, "_deseq2_results_upgenes.csv")
write.csv(out_up, file = file.path(output_dir_tables, filename))

# Export CSV file with gene symbols and entrez ids
filename <- paste0(project, "_deseq2_results_downgenes.csv")
write.csv(out_down, file = file.path(output_dir_tables, filename))


print(paste0("DESeq2 results were annotated and saved in ", output_dir_tables))
