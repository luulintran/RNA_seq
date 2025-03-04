# SET UP
library(DESeq2)
library(tidyverse)
library(dplyr)
library(readr)


# READ IN DATA:-----------------------------------------------------------------
# raw gene counts data from star_salmon; remove gene_name column
cts <- read.csv(input_cts, sep= "\t", row.names = 1) # makes the rownames the first col
cts <- as.data.frame(cts[, -c(1)]) # remove the gene_name col

# sample metadata
coldata <- read.csv(input_coldata, row.names = 1)

# CREATE DESEQ2 DATASET:--------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(cts), 
                              colData = coldata, 
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "control")

# RUN DESEQ2 ANALYSIS: ----------------------------------------------------------
dds <- DESeq(dds)

# SAVE DESEQ2 OBJECT TO RDS FOR LATER: -----------------------------------------
filename <- file.path(output_dir_robj, paste0(project, "_deseq2_dds.rds"))

saveRDS(dds, file = filename)

print(paste0("DESeq2 analysis done. RDS file saved in ", output_dir_robj))


