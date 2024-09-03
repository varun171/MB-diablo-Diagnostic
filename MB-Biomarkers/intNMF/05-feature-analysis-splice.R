# script to perform downstream feature level analysis for splice data

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggupset)
  library(wordcloud)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "multimodal_clustering")

# input directory
input_dir <- file.path(analysis_dir, "results")

# output directory
output_dir <- file.path(analysis_dir, "results", "feature_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots", "feature_analysis")
dir.create(plots_dir, showWarnings = F, recursive = T)

# read full output
nmf_output <- readRDS(file.path(input_dir, "intnmf_best_fit.rds"))

# read splice data and extract features with non-zero NMF weights per cluster
splice_selected <- nmf_output$H$H5
splice_selected <- splice_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  filter(value > 0) %>%
  mutate(feature = gsub("_.*", "", feature)) %>% # convert splice variant to gene name
  mutate(cluster = as.numeric(cluster)) # convert cluster to numeric

# get max cluster
cluster_max <- splice_selected %>% 
  pull(cluster) %>% 
  max()

# enrichment output in a list
go_enrich <- list()
for(i in 1:cluster_max){
  print(i)
  splice_genes <- splice_selected %>%
    filter(cluster == i) %>%
    pull(feature) %>% 
    unique()
  
  # perform over-representation analysis for each cluster based on gene terms with associated visuals (code below)
  # GSOA
  set.seed(100)
  go_enrich[[i]] <- clusterProfiler::enrichGO(gene = splice_genes,
                                         OrgDb = org.Hs.eg.db, 
                                         keyType = 'SYMBOL',
                                         readable = T,
                                         ont = "BP",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)
}

# Upset plot
pdf(file = file.path(plots_dir, "splice_gsoa.pdf"), width = 12, onefile = TRUE)
for(i in 1:cluster_max){
  go_enrich_cluster <- go_enrich[[i]] %>% 
    as.data.frame()
  if(nrow(go_enrich_cluster) > 0) {
    print(enrichplot::upsetplot(x = go_enrich[[i]]))
  } else {
    plot.new()
  }
}
dev.off()

# Word cloud 
pdf(file = file.path(plots_dir, "splice_wordcloud.pdf"), onefile = TRUE)
for(i in 1:cluster_max){
  go_enrich_cluster <- go_enrich[[i]] %>% 
    as.data.frame()
  if(nrow(go_enrich_cluster) > 0){
    wcdf <- read.table(text = go_enrich_cluster$GeneRatio, sep = "/")[1]
    wcdf$term <- go_enrich_cluster[,2]
    set.seed(100)
    wordcloud(words = wcdf$term, 
              freq = wcdf$V1, scale = (c(1, .1)), 
              colors = brewer.pal(8, "Dark2"), max.words = 25)
  } else {
    plot.new()
  }
}
dev.off()

#  ggplot by p-value
pdf(file = file.path(plots_dir, "splice_barplot.pdf"), onefile = TRUE, width = 10, height = 8)
for(i in 1:cluster_max){
  print(i)
  p.input <- go_enrich[[i]]@result %>% 
    dplyr::arrange(p.adjust) %>% 
    dplyr::filter(p.adjust < 0.05)
  if(nrow(p.input) > 0){
    p <- ggplot(p.input, aes(x = reorder(Description, p.adjust), y = Count, fill = p.adjust)) + 
      geom_bar(stat="identity") + coord_flip() + theme_bw(base_size = 10) +
      xlab("") + 
      ylab("Count") + 
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))  +
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
      ggtitle(paste0("Enrichment in Cluster", i))
    print(p)
  } else {
    plot.new()
  }
  
}
dev.off()

# Barplot: Graphics package
pdf(file = file.path(plots_dir, "splice_barplot_graphics.pdf"), onefile = TRUE, width = 8, height = 6)
for(i in 1:cluster_max){
  print(i)
  go_enrich_cluster <- go_enrich[[i]] %>% 
    as.data.frame()
  if(nrow(go_enrich_cluster) > 0){
    p <- barplot(go_enrich[[i]], 
            drop = TRUE, 
            showCategory = 10, 
            title = paste0("GO Biological Pathways in Cluster",i),
            font.size = 8)
    print(p)
  } else {
    plot.new()
  }
}
dev.off()

# Enrichplot 
pdf(file = file.path(plots_dir, "splice_goplot.pdf"), onefile = TRUE, width = 14, height = 12)
for(i in 1:cluster_max){
  print(i)
  go_enrich_cluster <- go_enrich[[i]] %>% 
    as.data.frame()
  if(nrow(go_enrich_cluster) > 0){
    p <- goplot(x = go_enrich[[i]], showCategory = 10)
    print(p)
  } else {
    plot.new()
  }
}
dev.off()
