# This script is to calculate cluster stats for each k
suppressPackageStartupMessages({
  library(CancerSubtypes)
  library(tidyverse)
  library(cluster)
})

run_clusterstats <- function(dat, wt = NULL, output_dir, maxK, clusterNum, K_range, alpha_range) {
  # Check if wt is NULL and assign default values
  if (is.null(wt)) {
    wt <- if (is.list(dat)) rep(1, length(dat)) else 1
  }
  
  # Define the file paths
  snf_fname <- file.path(output_dir, "snfcc_fit_all.rds")
  output_fname <- file.path(output_dir, "snfcc_clusterstats.tsv")
  
  # Initialize empty data frame
  output_df <- data.frame()
  
  # Check if snf_fname exists
  if (!file.exists(snf_fname)) {
    snf_output <- list()
    
    for (K in K_range) {
      for (alpha in alpha_range) {
        print(paste("Running with clusterNum =", "K =", K, "alpha =", alpha))
        var <- paste0("k", "_K", K, "_alpha", gsub("\\.", "_", as.character(alpha)))
        
        # ExecuteSNF.CC function call goes here
        
        # Add error handling for ExecuteSNF.CC and file operations
        
        snf_output[[var]] <- ExecuteSNF.CC(dat, clusterNum = clusterNum, K = K, alpha = alpha, t = 50, 
                                           maxK = maxK, pItem = 0.8, reps = 500, 
                                           title = "MB-Noclus", 
                                           plot = "pdf", finalLinkage = "average")
      }
    }
    
    # Save snf_output as an RDS file
    saveRDS(snf_output, file = snf_fname)
  } else {
    # Load snf_output from an existing RDS file
    snf_output <- readRDS(snf_fname)
  }
  
  # Iterate over K and alpha ranges
  for (K in K_range) {
    for (alpha in alpha_range) {
      var <- paste0("k", "_K", K, "_alpha", gsub("\\.", "_", as.character(alpha)))
      for (m in 2:20) {
        # Calculate silhouette width
        sil <- silhouette_SimilarityMatrix(snf_output[[var]][["originalResult"]][[m]][["consensusClass"]], snf_output[[var]][["originalResult"]][[m]][["consensusMatrix"]])
        sil <- summary(sil)
        sil_width <- sil[["avg.width"]]
        
        # Add the results to output_df
        output_df <- rbind(output_df, data.frame(k = var, m = m, avg_sil_width = sil_width))
      }
    }
  }
  
  # Save output_df as a TSV file
  write.table(output_df, file = output_fname, sep = "\t", row.names = FALSE)
  
  # Find the maximum avg_sil_width and corresponding k and m
  max_k <- output_df[which.max(output_df$avg_sil_width), ]
  
  # Check if max_k is not NULL
  if (!is.null(max_k)) {
    # Save the best fit to a file
    saveRDS(object = snf_output[[max_k$k]][["originalResult"]][[max_k$m]], file = file.path(output_dir, "snfcc_best_fit.rds"))
    
    # Return the snf output corresponding to the most optimal k for downstream analyses
    return(snf_output[[max_k$k]][["originalResult"]][[max_k$m]])
  } else {
    # Return a message indicating that no optimal k was found
    return("No optimal k found.")
  }
}
