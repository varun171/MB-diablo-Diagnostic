# function to compute v.test
v_test <- function(x, clustering_col){
  
  # get unique clusters i.e. 1-n clusters
  clusters <- unique(x[[clustering_col]])
  
  out <- data.frame()
  # iterate over each cluster
  for(i in 1:length(clusters)){
    
    # subset x to cluster i
    y = x %>% 
      filter(get(clustering_col) == clusters[i])
    
    # mean expression per cluster per gene  
    mean.gene.cluster = unique(y$cluster_gene_mean_score)
    
    # mean expression per gene (i.e. global mean)
    mean.gene = unique(y$gene_mean_score)
    
    # calculate numerator
    num = mean.gene.cluster - mean.gene 
    
    # total sample size
    n = nrow(x)
    
    # cluster sample size
    ng = nrow(y)
    
    # variance of expression per gene (i.e. global variance)
    var.gene = unique(x$gene_variance) 
    
    # calculate denominator
    denom = (n-ng/n-1)*(var.gene/ng)
    denom = sqrt(denom)
    
    # calculate vscore
    v = num/denom
    out[i,'cluster'] <- clusters[i]
    out[i,'v_score'] <- v
  }
  return(out)
}

# function to compute means and variance global and per group
compute_v_stats <- function(df_with_clusters, cluster_col, n_features){
  df_with_clusters  <- df_with_clusters %>%
    group_by(get(cluster_col), gene) %>%
    dplyr::mutate(cluster_gene_mean_score = mean(expression)) %>% # mean of gene expression per cluster per gene
    ungroup() %>%
    group_by(gene) %>%
    dplyr::mutate(gene_mean_score = mean(expression),
                  gene_variance = var(expression)) # global mean & variance of expression per gene
  
  # apply v-test function per gene
  out <- plyr::ddply(.data = df_with_clusters, 
                     .variables = "gene", 
                     .fun = function(x) v_test(x, clustering_col = cluster_col))
  
  # only select n_features per cluster
  best_features_per_cluster <- out %>%
    group_by(gene) %>%
    filter(v_score == max(v_score)) %>%
    group_by(cluster) %>%
    dplyr::arrange(desc(v_score)) %>%
    dplyr::slice_head(n = n_features)
  
  # create a matrix of v-test scores for top 10 features per cluster
  out <- out %>%
    spread(cluster, v_score) %>%
    filter(gene %in% best_features_per_cluster$gene) %>%
    column_to_rownames("gene")
  
  # get the max cluster and add as a column
  out$best_fit <- colnames(out[,-1])[max.col(out[,-1], ties.method = "first")]
  
  return(out)
}
