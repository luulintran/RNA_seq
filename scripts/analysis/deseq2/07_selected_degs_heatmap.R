# READ DDS OBJECT FOLLOWING DESEQ2 ANALYSIS FROM RDS FILE: --------------------
dds <- readRDS(
  file = file.path(output_dir_robj, paste0(project, "_deseq2_dds.rds"))
)

# Store results in res
res <- results(dds)

# MAKE GENE LISTS: -------------------------------------------------------------
# Gene lists are based on cell type markers to determine 
# what cell identity progenitors are based on gene expression
progenitor_genes <- c('Fabp7', 'Nes', 'Pax6', 'Slc1a3', 'Sox2', 
                      'Vim', 'Nr2e1', 'Hes1', 'Hes5', 'Ednrb', 
                      'Eomes', 'Ccne2', 'Clspn', 'Gins2', 'Pcna', 
                      'Atad2', 'Mcm7', 'Mcm3', 'Slbp', 'Gmnn', 
                      'Kiaa0101', 'Mcm10', 'Rad51', 'Cdc45', 'Exo1', 
                      'Hist1h4c', 'Cdk1', 'Hist1h1b', 'Hist1h1c', 
                      'Hist1h1e', 'Ube2c', 'Rrm2', 'Zwint', 'Hmgb2', 
                      'Ccna', 'Cdca5', 'Esco2', 'Aurkb', 
                      'Kif18b', 'Ckap2l', 'Hjurp', 'Cdca8', 'Ccnb1', 
                      'Cenpf', 'Cks2', 'Pttg1', 'Cdc20', 'Top2a', 
                      'Nusap1', 'Cenpa', 'Psrc1', 'Gas2l3', 'Plk1', 
                      'Kif20a')

RGC_genes <- c('Fabp7', 'Nes', 'Pax6', 'Slc1a3', 'Sox2', 'Vim', 
               'Nr2e1', 'Hes1', 'Hes5', 'Ednrb') 

IPC_genes <- c('Eomes', 'Sema3c', 'Neurod1', 'Neurog2', 'Sstr2', 
               'Gadd45g', 'Neurog1')

proliferative_genes <- c('Fabp7', 'Nes', 'Pax6', 'Slc1a3', 'Sox2', 
                         'Vim', 'Nr2e1')

neurogenic_genes <- c('Eomes', 'Neurog2', 'Neurog', 'Tuba1a', 'Btg2')

neuronal_genes <- c('Map2', 'Mapt', 'Rbfox3', 'Tbr1', 'Tubb3', 'Neurod6', 
                    'Neurod2','Satb2', 'Gria2', 'Nrp1', 'Dab1', 
                    'Nrxn3', 'Neurod4')

newborn_neurons <- c('Foxg1', 'Neurod1', 'Unc5d', 'Rnd2', 'Rnd3', 'Dcx', 
                     'Pafah1b1', 'Cdk5')

OPC_genes <- c('Sox10', 'Pdgfra', 'Olig1', 'Olig2', 'Ascl1', 'Gng12', 
               'Cnp', 'Cspg4', 'Matn4', 'Brinp3', 
               'Lhfpl3', 'Cntn1')

preOPC_genes <- c('Ascl1', 'Ccnd1', 'Dleu7', 'Egfr', 'Egr1', 'Qk', 'Gas1', 
                  'Sall3', 'Gng12', 'Ncald', 'Gsx2', 
                  'Fam181b', 'Rfx4', 'Bcan')

astrocyte_genes <- c('Aldh1l1', 'Fabp7', 'Aldoc', 'Hes5', 'Aqp4', 'Sox9')


# Combine some of the lists to be more general
neurogenic_all <- c(IPC_genes, neurogenic_genes, neuronal_genes, 
                    newborn_neurons)

gliogenic_all <- c(OPC_genes, preOPC_genes, astrocyte_genes)

progenitor_all <- c(progenitor_genes, RGC_genes, proliferative_genes)

opc_all <- c(preOPC_genes, OPC_genes)

# ANNOTATE DESEQ2 RESULTS WITH GENE SYMBOLS AND ENTREZ IDS: --------------------
ensembl_ids <- rownames(res)

# annotate with gene symobols using org.Mm.eg.db package
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
                                     c("symbol", "entrez")))
]


# PREPARE DATA FOR HEATMAP SHOWING SPECIFIC GENE SET: --------------------------

# Combine progenitor, neurogenic, and gliogenic gene lists for heatmap
combined_gene_list <- c(progenitor_all, gliogenic_all, neurogenic_all)

# Remove NAs from padj before filtering for significant genes with padj < 0.05
sig_res <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05, ]

# Extract rownames/ensembl id's of significant DEGs
sig_genes <- rownames(sig_res)

# Map ensembl id's to symbols
sig_symbols <- mapIds(org.Mm.eg.db,
                      keys = sig_genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Filter sig_genes list for genes in combined_gene list by symbol
selected_genes <- sig_genes[sig_symbols %in% combined_gene_list]

# EXTRACT TRANSFORMED VALUES: **************************************************
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% combined_gene_list]

# Center the data
mat <- mat - rowMeans(mat)


# MAKE HEATMAP: ---------------------------------------------
default_clustered_heatmap <- pheatmap(
  mat,
  cluster_rows = TRUE,  # Default clustering behavior
  cluster_cols = TRUE, 
  show_rownames = TRUE,  
  annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c("#76BAE0", "white", "#B8396B"))(50)
)

# Save heatmap
filename <- file.path(output_dir_figures, paste0(project, "_degs_heatmap_default.png"))
png(filename, width = 5, height = 8, units = "in", res = 300)
print(default_clustered_heatmap)
dev.off()

print(paste0("Heatmap saved in ", output_dir_figures))

# MAKE HEATMAP FOR NOTCH GENES: ------------------------------------------------

genes_list <- c("Notch1", "Notch2", "Hes1", "Hes5", "Hey1", "Hes6", "Notch3", 
                "Notch4", "Rbpj", "Maml", "Dll1", "Dll3", "Jag1","Jag2")

# Remove NAs from padj before filtering for significant genes with padj < 0.05
sig_res <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05, ]

# Extract rownames/ensembl id's of significant DEGs
sig_genes <- rownames(sig_res)

# Map ensembl id's to symbols
sig_symbols <- mapIds(org.Mm.eg.db,
                      keys = sig_genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Filter sig_genes list for genes in combined_gene list by symbol
selected_genes <- sig_genes[sig_symbols %in% genes_list]

# EXTRACT TRANSFORMED VALUES: **************************************************
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% genes_list]

# Center the data
mat <- mat - rowMeans(mat)


# MAKE HEATMAP: ---------------------------------------------
default_clustered_heatmap <- pheatmap(
  mat,
  cluster_rows = TRUE,  # Default clustering behavior
  cluster_cols = TRUE, 
  show_rownames = TRUE,  
  annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c("#76BAE0", "white", "#B8396B"))(50)
)

# Save heatmap
filename <- file.path(output_dir_figures, paste0(project, "_degs_heatmap_notch.png"))
png(filename, width = 5, height = 8, units = "in", res = 300)
print(default_clustered_heatmap)
dev.off()

print(paste0("Heatmap saved in ", output_dir_figures))

# MAKE HEATMAP FOR SHH GENES: ------------------------------------------------

genes_list <- c("Smo", "Boc", "Cdon", "Gas1", "Gli1", "Gli2", "Gli3", "Sufu", 
                "Disp1", "Iqce", "Efcab7", "Ptch1", "Ptch2", "Hhip", "Sufu", 
                "Gsk3b", "Ck1", "Pcaf", "Cul3", "Kif7")

# Remove NAs from padj before filtering for significant genes with padj < 0.05
sig_res <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05, ]

# Extract rownames/ensembl id's of significant DEGs
sig_genes <- rownames(sig_res)

# Map ensembl id's to symbols
sig_symbols <- mapIds(org.Mm.eg.db,
                      keys = sig_genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Filter sig_genes list for genes in combined_gene list by symbol
selected_genes <- sig_genes[sig_symbols %in% genes_list]

# EXTRACT TRANSFORMED VALUES: **************************************************
vsd <- vst(dds, blind=FALSE)

# Subset assay matrix for the selected genes
mat <- assay(vsd)[selected_genes, ]

# Update row names with gene symbols for the selected genes
rownames(mat) <- sig_symbols[sig_symbols %in% genes_list]

# Center the data
mat <- mat - rowMeans(mat)


# MAKE HEATMAP: ---------------------------------------------
default_clustered_heatmap <- pheatmap(
  mat,
  cluster_rows = TRUE,  # Default clustering behavior
  cluster_cols = TRUE, 
  show_rownames = TRUE,  
  annotation_col = as.data.frame(colData(vsd)[, "condition", drop=FALSE]),
  color = colorRampPalette(c("#76BAE0", "white", "#B8396B"))(50)
)

# Save heatmap
filename <- file.path(output_dir_figures, paste0(project, "_degs_heatmap_shh.png"))
png(filename, width = 5, height = 8, units = "in", res = 300)
print(default_clustered_heatmap)
dev.off()

print(paste0("Heatmap saved in ", output_dir_figures))