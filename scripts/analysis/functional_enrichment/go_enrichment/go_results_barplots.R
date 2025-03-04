
library(clusterProfiler)


# MAKE BARPLOTS: ---------------------------------------------------------------
# UP GENES: --------------------------------------------------------------------
filename <- file.path(output_dir_figures, 
                      "GO_results_barplot_upgenes.png")

png(filename, width = 10, height = 10, units = "in", res = 300)

upgenes_barplot <- barplot(GO_results_up, showCategory = 20)
print(upgenes_barplot)

dev.off()


# DOWN GENES: ------------------------------------------------------------------
filename <- file.path(output_dir_figures, 
                      "GO_results_barplot_downgenes.png")

png(filename, width = 10, height = 10, units = "in", res = 300)

downgenes_barplot <- barplot(GO_results_down, showCategory = 20)
print(downgenes_barplot)

dev.off()

print(paste0("GO results saved in ", output_dir_figures))