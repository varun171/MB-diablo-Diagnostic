# This script is to calculate cluster stats for each k
suppressPackageStartupMessages({
  library(CancerSubtypes)
  library(tidyverse)
})

run_clusterstats <- function(dat, wt = NULL, output_dir, clusterNum){
  
  # if weight is not assigned, use default
  if(is.null(wt)){
    wt = if(is.list(dat)) rep(1,length(dat)) else 1
  }
  
  # do this for each cluster
  # run SNFCC of Multiple data 
  snf_fname <- file.path(output_dir, "snfcc_fit_all.rds")
  
  if(!file.exists(snf_fname)){
    snf_output <- list() # snfcc output
    for(i in 4:clusterNum){
      print(i)
      var <- paste0("k", i)
      snf_output[[var]] <- ExecuteSNF.CC(dat, clusterNum = i, K = 20, 
                                         alpha = 0.5, t = 50, 
                                         maxK = 20, pItem= 0.8, reps=500, 
                                         title = "MB", 
                                         plot = "png", finalLinkage ="average")
    }
    # save final nmf output for all k-values
    saveRDS(snf_output, file = snf_fname)
  } else {
    snf_output <- readRDS(snf_fname)
  }
  
  # do for each cluster
  output_fname <- file.path(output_dir, "snfcc_clusterstats.tsv")
  output_df <- data.frame() # output df for clusterstats
  for(j in 1:length(dat)){
    print(j)
    input_dat <- dat[[j]]
    count_matrix_dist <- as.matrix(factoextra::get_dist(t(input_dat), method = 'pearson'))
    count_matrix_dist[is.na(count_matrix_dist)] <- 0
    for(i in 4:clusterNum){
      print(i)
      var <- paste0("k", i)
      clus_stats <- fpc::cluster.stats(d = count_matrix_dist, clustering = snf_output[[var]]$group, silhouette = FALSE)
      clus_stats <- lapply(clus_stats, FUN = function(x){
        if(length(x) > 1 & !is.null(x)) {
          x <- toString(round(x, 2))
        } else {
          x <- round(as.numeric(x), digits = 2)
        }
      })
      clus_stats <- t(as.data.frame(unlist(clus_stats)))
      clus_stats <- clus_stats %>%
        as.data.frame() %>%
        dplyr::mutate(k = var,
                      dataset = j) %>%
        dplyr::select(k, dataset, cluster.size, sindex, noisen, diameter, average.distance, median.distance, average.between, average.within, within.cluster.ss, dunn, dunn2, entropy, wb.ratio)
      # combine with other k values
      output_df <- rbind(output_df, clus_stats)
    }
  }
  
  # save dataframe with cluster stats for all k-values
  write_tsv(output_df, file = output_fname)
  
  max_k <- output_df %>% 
    group_by(k) %>% 
    dplyr::summarise(n = sum(as.numeric(average.between)) - sum(as.numeric(average.within))) %>%
    filter(n == max(n)) %>%
    pull(k)
  
  # If there are ties
  sil_score_df <- data.frame()
  if(length(max_k) > 1){
    for(i in 1:length(max_k)){
      fit <- snf_output[[max_k[i]]]
      summary <- CancerSubtypes::silhouette_SimilarityMatrix(fit$group, fit$distanceMatrix) %>% summary
      average_width <- summary$sil.width
      tmp <- data.frame(k = max_k[i], sil_score = average_width)
      sil_score_df <- rbind(sil_score_df, tmp)
    }
    # get the k with max average silhouette width
    max_k <- sil_score_df[which.max(sil_score_df$sil_score),"k"]
  }
  
  # save best fit to file
  saveRDS(object = snf_output[[max_k]], file = file.path(output_dir, "snfcc_best_fit.rds"))
  
  # return the nmf output corresponding to the most optimal k for downstream analyses
  return(snf_output[[max_k]])
}
