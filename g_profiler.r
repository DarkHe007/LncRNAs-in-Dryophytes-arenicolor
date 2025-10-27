#!/usr/bin/env Rscript

library(gprofiler2)

# Simple g:Profiler analysis
stages <- c("Pre", "Climax", "Adult")

for (stage in stages) {
  # Read significant hits
  hits_file <- paste0("filtered_hits/", stage, "_significant_ids.txt")
  if (file.exists(hits_file)) {
    genes <- readLines(hits_file)
    genes <- unique(genes[genes != ""])
    
    if (length(genes) > 0) {
      # Run g:Profiler
      result <- gost(genes, organism = "hsapiens")  
      
      # Save results
      if (!is.null(result$result)) {
        write.csv(result$result, 
                  paste0("gprofiler_results/", stage, "_simple.csv"),
                  row.names = FALSE)
        message(paste(stage, ":", nrow(result$result), "enriched terms"))
      }
    }
  }
}