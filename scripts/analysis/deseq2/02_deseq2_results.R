# STORE DESEQ2 RESULTS: --------------------------------------------------------
res <- results(dds)

# Make it into a dataframe and save for later
out <- as.data.frame(res)

# Write to file
filename <- paste0(project, "_deseq2_results.csv")

write.csv(out, 
            file=file.path(output_dir_tables, filename)
          )

print(paste0(filename, " saved in ", output_dir_tables))