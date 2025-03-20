# MAKE LIST OF GENES RELATED TO FILTER: -----------------------------------
genes_list <- c("Smo", "Boc", "Cdon", "Gas1", "Gli1", "Gli2", "Gli3", "Sufu", 
                "Disp1", "Iqce", "Efcab7", "Ptch1", "Ptch2", "Hhip", "Sufu", 
                "Gsk3b", "Ck1", "Pcaf", "Cul3", "Kif7")

# READ IN DESeq2 RESULTS CSV FILE: ---------------------------------------------
assign(
  "resOrdered",  
  read.csv(file.path(output_dir_tables, paste0(project, "_deseq2_results_annotated.csv"))) 
)

# Save as dataframe
out <- as.data.frame(resOrdered)

# FILTER DESEQ2 RESULTS FOR SHH GENES: -----------------------------------------

genes_list_deseq2_results <- out %>%
  filter(symbol %in% genes_list)


# Make rownames a column and name ensembl because it contains ensembl id's
genes_list_deseq2_results <- rownames_to_column(genes_list_deseq2_results, var = "ensembl")

# SAVE TO FILE: ----------------------------------------------------------------
filename <- paste0(project, "_deseq2_filtered_genes.csv")

write.csv(
  genes_list_deseq2_results, 
  file = file.path(output_dir_tables, filename),
  row.names = F)

