# This script is to calculate cluster stats for each k
suppressPackageStartupMessages({
  library(IntNMF)
  library(tidyverse)
})

run_clusterstats <- function(dat, wt = NULL, output_dir, k_value){
  
  # if weight is not assigned, use default
  if(is.null(wt)){
    wt = if(is.list(dat)) rep(1,length(dat)) else 1
  }
  
  nmf_fname <- file.path(output_dir, "intnmf_fit_all.rds")
  if(!file.exists(nmf_fname)){
    nmf_output <- list() # nmf mnnals output
    for(i in 4:k_value){
      print(i)
      var <- paste0("k", i)
      nmf_output[[var]] <- nmf.mnnals(dat = dat, maxiter = 200, st.count = 20, wt = wt, k = i, seed = TRUE)
    }
    saveRDS(nmf_output, file = nmf_fname)
  } else {
    nmf_output <- readRDS(nmf_fname)
  }
  
  output_fname <- file.path(output_dir, "intnmf_clusterstats.tsv")
  output_df <- data.frame() # output df for clusterstats
  
  for (m in 4:k_value) {
    var <- paste0("k", m)
    fit <- nmf_output[[var]]
    distance <- 1 - fit$consensus
    summary <- cluster::silhouette(fit$clusters, dmatrix = distance) %>% summary
    average_width <- summary$avg.width
    
    output_df <- rbind(output_df, data.frame(k = var, m = m, avg_sil_width = average_width))
  }
  
  write_tsv(output_df, file = output_fname)
  
  max_k <- output_df[which.max(output_df$avg_sil_width), ]
  
  if (!is.null(max_k)) {
    best_fit_fname <- file.path(output_dir, "nmf_best_fit.rds")
    saveRDS(object = nmf_output[[max_k$k]], file = best_fit_fname)
    return(nmf_output[[max_k$k]])
  } else {
    return("No optimal k found.")
  }
}
