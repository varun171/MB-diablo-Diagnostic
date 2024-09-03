# This script is to calculate cluster stats for each k and alpha
suppressPackageStartupMessages({
  library(CancerSubtypes)
  library(tidyverse)
})
run_clusterstats <- function(dat, wt, output_dir, clusterNum) {
  
  # if weight is not assigned, use default
  if(is.null(wt)){
    wt = if(is.list(dat)) rep(1,length(dat)) else 1
  }
  
  snf_fname <- file.path(output_dir, "snfcc_fit_all.rds")
  
  if (!file.exists(snf_fname)) {
    snf_output <- list()
    
    for (alpha in seq(0.3, 0.8, by = 0.1)) {
      alpha_var <- as.character(alpha)
      
        for (i in 4:clusterNum) {
        print(paste("alpha =", alpha, "k =", i))
        var <- paste0("alpha_", alpha_var, "_k", i)
        snf_output[[alpha_var]][[var]] <- ExecuteSNF.CC(dat, clusterNum = i, K = 20,
                                                              alpha = alpha, t = 50,
                                                              maxK = 20, pItem = 0.8, reps = 500,
                                                              title = "clusstat",
                                                              plot = "png", finalLinkage = "average")
      }
    }
    saveRDS(object = snf_output, file = snf_fname)
  } else {
    snf_output <- readRDS(snf_fname)
  }
  
  output_df <- data.frame()
  
  for (alpha in seq(0.3, 0.8, by = 0.1)) {
    for (j in 1:length(dat)) {
      input_dat <- dat[[j]]
      count_matrix_dist <- as.matrix(factoextra::get_dist(t(input_dat), method = 'pearson'))
      count_matrix_dist[is.na(count_matrix_dist)] <- 0
      
      for (i in 4:clusterNum) {
        alpha_var <- as.character(alpha)
        var <- paste0("alpha_", alpha_var, "_k", i)
        clus_stats <- fpc::cluster.stats(d = count_matrix_dist, clustering = snf_output[[alpha_var]][[var]]$group, silhouette = FALSE)
        clus_stats <- lapply(clus_stats, FUN = function(x) {
          if (length(x) > 1 & !is.null(x)) {
            x <- toString(round(x, 2))
          } else {
            x <- round(as.numeric(x), digits = 2)
          }
        })
        clus_stats <- t(as.data.frame(unlist(clus_stats)))
        clus_stats <- clus_stats %>%
          as.data.frame() %>%
          dplyr::mutate(alpha = alpha,
                        k = i,
                        dataset = j) %>%
          dplyr::select(alpha, k, dataset, cluster.size, sindex, noisen, diameter, average.distance, median.distance, average.between, average.within, within.cluster.ss, dunn, dunn2, entropy, wb.ratio)
        output_df <- rbind(output_df, clus_stats)
      }
    }
  }
  
  write.table(output_df, file = file.path(output_dir, "snfcc_clusterstats.tsv"), sep = "\t", row.names = FALSE)
  
  max_k <- output_df %>% 
    group_by(alpha, k) %>% 
    dplyr::summarise(n = sum(as.numeric(average.between)) - sum(as.numeric(average.within))) %>%
    filter(n == max(n)) %>%
    pull(alpha,k)
  
  sil_score_df <- data.frame()
  if (length(max_k) > 1) {
    for (i in 1:length(max_k)) {
      fit <- snf_output[[as.character(max_k[i])]]
      summary <- CancerSubtypes::silhouette_SimilarityMatrix(fit$group, fit$distanceMatrix) %>% summary
      average_width <- summary$sil.width
      tmp <- data.frame(k = max_k[i], sil_score = average_width)
      sil_score_df <- rbind(sil_score_df, tmp)
    }
    max_k <- sil_score_df[which.max(sil_score_df$sil_score), "k"]
  }
  
  # Save best fit to file
  saveRDS(object = snf_output[[as.character(max_k)]], file = file.path(output_dir, "snfcc_best_fit.rds"))
  
  # Return the SNF output corresponding to the most optimal k for downstream analyses
  return(snf_output[[as.character(max_k)]])
}
