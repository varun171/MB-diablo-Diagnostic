# Author: Komal S. Rathi
# Function: Connectivity Analysis using LINCS

suppressPackageStartupMessages({
  library(signatureSearch)
  library(tidyverse)
  library(dplyr)
  library(optparse)
  library(ExperimentHub)
  library(rhdf5)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(visNetwork)
  library(igraph)
  library(cowplot)
  library(gridExtra)
  library(ggpubr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "multimodal_clustering", "lincs_analysis")
data_dir <- file.path(analysis_dir, "input")

# source functions
source(file.path(analysis_dir, "utils", "network_to_file.R"))
source(file.path(analysis_dir, "utils", "drug_barplots.R"))

# get LINCS data
eh <- ExperimentHub()
lincs <- eh[["EH3226"]]

# read oncology drugs of interest
# drugs_of_interest <- read_tsv(file.path(data_dir, "oncology_drugs_for_LINCS.tsv"), col_names = F)
# ref_trts <- drugs_of_interest$X1

# read touchstone compounds
touchstone_drugs <- read_tsv(file.path(data_dir, "touchstone_data.txt"))

# run connectivity analysis on each cluster
lincs_connectivity <- function(input, num_features = 2000, num_sets = 25, method = c("LINCS", "Cor"), wtcs_fdr_cutoff = 0.05, trend_val = c("up", "down"), drug_type = c("reversers", "mimickers"), cor_score_cutoff = 0, prefix, output_dir, plots_dir){
  
  # directories
  dir.create(output_dir, showWarnings = F, recursive = T)
  dir.create(plots_dir, showWarnings = F, recursive = T)
  
  # up and down genes
  cluster_upset <- input %>%
    filter(direction == "up") %>%
    arrange(log2FoldChange) %>%
    slice_head(n = num_features)
  cluster_downset <- input %>%
    filter(direction == "down") %>%
    arrange(log2FoldChange) %>%
    slice_head(n = num_features)
  
  if(method == "LINCS"){
    # map to ENTREZ identifiers
    if(nrow(cluster_upset) > 1){
      upset = mapIds(org.Hs.eg.db, keys = cluster_upset$genes, column = "ENTREZID", keytype = "SYMBOL")
    } else {
      upset <- ""
    }
    if(nrow(cluster_downset) > 1){
      downset = mapIds(org.Hs.eg.db, keys = cluster_downset$genes, column = "ENTREZID", keytype = "SYMBOL")
    } else {
      downset <- ""
    }
    
    # LINCS-based similarity metric
    qSig_output <- qSig(query = list(upset = upset, downset = downset), gess_method = "LINCS", refdb = lincs)
    qSig_output <- gess_lincs(qSig = qSig_output, sortby = "NCS", tau = T, workers = 4)
    qSig_output <- result(qSig_output)
    print("Before filtering:")
    print(dim(qSig_output))
    
    # filter drugs on trend and FDR
    qSig_output <- qSig_output %>%
      filter(trend == trend_val,
             WTCS_FDR < wtcs_fdr_cutoff)

    # filter based on NCS value
    if(drug_type == "mimickers"){
      qSig_output <- qSig_output %>%
        filter(NCS > 0)
    } else if(drug_type == "reversers"){
      qSig_output <- qSig_output %>%
        filter(NCS < 0)
    }

  } else if(method == "Cor"){
    # correlation-based similarity
    query_mat <- cluster_upset %>%
      plyr::rbind.fill(cluster_downset) %>% 
      mutate(id = mapIds(org.Hs.eg.db, keys = geneSymbol, column = "ENTREZID", keytype = "SYMBOL")) %>% 
      arrange(id) %>% 
      filter(!is.na(id)) %>%
      column_to_rownames('id') %>%
      dplyr::select(score) %>% as.matrix()
    qSig_output <- qSig(query = query_mat, gess_method = "Cor", refdb = lincs)
    qSig_output <- gess_cor(qSig = qSig_output, method = "spearman", workers = 4)
    qSig_output <- result(qSig_output)
    print("Before filtering:")
    print(dim(qSig_output))
    
    # filter drugs
    qSig_output <- qSig_output %>%
      filter(cor_score < cor_score_cutoff)
  }
  
  print("After filtering:")
  print(dim(qSig_output))
  drugs <- unique(qSig_output$pert)
  
  # add touchstone annotation
  qSig_output <- qSig_output %>%
    mutate(Touchstone = ifelse(pert %in% touchstone_drugs$Name, TRUE, FALSE))
  
  # write output
  fname <- file.path(output_dir, paste0(prefix, "_qSig_output.txt"))
  write.table(qSig_output, file = fname, quote = F, sep = "\t", row.names = F) 
  
  # default plots are set to NULL and only updated where data is available
  p1 <- ggplot()
  p2 <- ggplot()
  p3 <- ggplot()
  p4 <- ggplot()
  if(length(drugs) != 0){
    
    # Query signature drug barplot
    p1 <- drug_barplots(dat = qSig_output, 
                        xlab = "pert", ylab = "WTCS_FDR",
                        top = 20, fill_var = NULL,
                        title = "Query Signature")
    
    # Top drug barplot 
    p2 <- drug_barplots(dat = qSig_output %>% filter(Touchstone), 
                        xlab = "pert", ylab = "WTCS_FDR",
                        top = 20, fill_var = NULL,
                        title = "Touchstone Signature")
    
    # hypergeometric TSEA using Reactome 
    tsea_reactome <- tsea_dup_hyperG(drugs = drugs, 
                                     type = "Reactome", 
                                     pvalueCutoff = 0.5, 
                                     dt_anno = 'DrugBank',
                                     qvalueCutoff = 0.5, readable = TRUE)
    if(!is.null(tsea_reactome)){
      tsea_reactome_df <- result(tsea_reactome)
      fname <- file.path(output_dir, paste0(prefix, "_tsea_reactome_output.txt"))
      write.table(tsea_reactome_df %>% mutate(Description = gsub("\r:", ":", Description)), file = fname, quote = F, sep = "\t", row.names = F) 
      
      if(nrow(tsea_reactome_df) > 0){
        # TSEA drug barplot
        tsea_reactome_df <- tsea_reactome_df %>% 
          mutate(Description = gsub("Homo sapiens\r: ", "", Description))
        p3 <- drug_barplots(dat = tsea_reactome_df, 
                            xlab = "Description", ylab = "p.adjust",
                            top = 20, fill_var = NULL,
                            title = "TSEA Reactome")
      }
    }
    
    # hypergeometric DSEA using GO MF
    dsea_go_mf <- dsea_hyperG(drugs = drugs, type = "GO", ont = "MF")
    if(!is.null(dsea_go_mf)){
      dsea_go_mf_df <- result(dsea_go_mf)
      fname <- file.path(output_dir, paste0(prefix, "_dsea_go_mf_output.txt"))
      write.table(dsea_go_mf_df, file = fname, quote = F, sep = "\t", row.names = F)
      
      if(nrow(dsea_go_mf_df) > 0){
        # network plot
        top_set <- dsea_go_mf_df %>% 
          arrange(p.adjust) %>% 
          slice_head(n = num_sets) %>% 
          pull(ID)
        if(length(top_set) > 0){
          dnet_object <- dtnetplot(drugs = drugs(dsea_go_mf), set = top_set, ont = "MF")
          # don't save html output
          # fname <- file.path(plots_dir, paste0(prefix, "_dsea_go_mf_output.html"))
          # visSave(dnet_object, fname) 
          fname <- file.path(plots_dir, paste0(prefix, "_dsea_go_mf_output.pdf"))
          network_to_file(dnet_object = dnet_object, filename = fname)
        }
        
        # DSEA drug barplot
        p4 <- drug_barplots(dat = dsea_go_mf_df, 
                            xlab = "Description", ylab = "p.adjust",
                            top = 20, fill_var = NULL,
                            title = "DSEA GO MF")
      }
    } 
  } 
  
  # return list of plots per comparison
  p <- list("Query Signature" = p1, 
            "Touchstone Signature" = p2,
            "TSEA Reactome" = p3, 
            "DSEA GO MF" = p4)
  
  # output barplots
  pdf(file = file.path(plots_dir, paste0(prefix, "_drug_pathways_barplot.pdf")))
  print(p)
  dev.off()
}

