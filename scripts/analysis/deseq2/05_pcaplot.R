# SET UP
library(ggplot2)
library(DESeq2)
library(extrafont)

# import Arial font
font_import(pattern = "Arial", prompt = FALSE)

# READ DDS OBJECT FROM RDS FILE: -----------------------------------------------
dds <- readRDS(
  file = file.path(output_dir_robj, paste0(project, "_deseq2_dds.rds"))
)

# EXTRACT TRANSFORMED VALUES: --------------------------------------------------
vsd <- vst(dds, blind=FALSE)

# EXTRACT PCA DATA FROM VSD: ---------------------------------------------------
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# get percentage of variance for first two principal components
percentVar <- round(100 * attr(pcaData, "percentVar"))

# DEFINE TREATMENT GROUPS: -----------------------------------------------------
# Define the colors
control_color <- "#76BAE0"
mutant_color <- "#B8396B"

# PLOT CUSTOM PCA: -------------------------------------------------------------
pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA",
       x = paste0("PC1[", percentVar[1], "%]"),
       y = paste0("PC2[", percentVar[2], "%]")) +
  scale_color_manual(values = c("control" = control_color, 
                                "NICD" = mutant_color)) + 
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# SAVE PLOT: -------------------------------------------------------------------
filename <- file.path(output_dir_figures, "pcaplot.png")

png(filename, width = 3, height = 3, units = "in", res = 300)
print(pcaplot)

dev.off()

print(paste0("PCA plot saved in ", output_dir_figures))