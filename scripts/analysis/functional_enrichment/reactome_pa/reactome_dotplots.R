library(ReactomePA)

# MAKE DOTPLOTS: ---------------------------------------------------------------

# Upregulated genes:
filename <- file.path(output_dir_figures, 
                      paste0("kegg_dotplot_upgenes.png"))

png(filename, width = 10, height = 10, units = "in", res = 300)

up_dotplot <- dotplot(pathway_up)
print(up_dotplot)

dev.off()


# Downregulated genes: 
filename <- file.path(output_dir_figures, 
                      paste0("kegg_dotplot_downgenes.png"))

png(filename, width = 10, height = 10, units = "in", res = 300)

down_dotplot <- dotplot(pathway_down)
print(down_dotplot)

dev.off()

print(paste0("KEGG results saved in ", output_dir_figures))

