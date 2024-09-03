suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
})

fgsea_barplot <- function(input_dat, gene_set, output_file){
  pdf(file = output_file, width = 14, height = 10, onefile = T)
  # run for each cluster and save to one pdf file
  cluster_max <- input_dat %>% pull(cluster) %>% max()
  for(i in 1:cluster_max){
    print(i)
    input_dat_cluster <- input_dat %>%
      filter(cluster == i)
    ranks <- input_dat_cluster$value
    names(ranks) <- input_dat_cluster$feature
    
    # fgsea using reactome dataset
    set.seed(42)
    res_cluster <- fgsea(gene_set, ranks, minSize = 5, maxSize = 1500, eps = 0.0)
    res_cluster$cluster <- i
    
    # filter to significant results
    res_cluster <- res_cluster %>%
      mutate(direction = ifelse(ES > 0, "Up", "Down")) %>%
      arrange(direction, padj) %>%
      filter(padj < 0.05) %>%
      dplyr::select(pathway, padj, ES, direction, cluster) %>%
      unique()
    
    # pull top 50 pathways per direction by adjusted p-value to reduce barplots size
    res_cluster <- res_cluster %>% 
      group_by(direction) %>%
      dplyr::arrange(padj) %>% 
      slice_head(n = 50)
    
    # set levels
    print(dim(res_cluster))
    res_cluster$pathway <- factor(res_cluster$pathway, levels = unique(res_cluster$pathway))
    res_cluster$direction <- factor(res_cluster$direction, levels = c("Up", "Down"))
    
    # barplot
    if(nrow(res_cluster) > 0){
      p <- ggplot(res_cluster, aes(pathway, y = (-1)*log10(padj), fill = direction)) + 
        geom_bar(stat="identity") + coord_flip() + theme_bw() +
        xlab("") + 
        ylab("-log10 Adj. P-Value") + 
        scale_fill_manual(name = "Direction", values = c("Down" = "forest green", "Up" = "red")) +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
        scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
        ggtitle(paste("Enrichment for Cluster:", i))
      print(p)
    }
  }
  dev.off()
}
