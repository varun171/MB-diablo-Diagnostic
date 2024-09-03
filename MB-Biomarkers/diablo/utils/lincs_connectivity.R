# run connectivity analysis on each cluster
# adapted from Komal Rathi's script in d3b-patient-report-analysis

suppressPackageStartupMessages({
  library(signatureSearch)
  library(AnnotationDbi)
  library(visNetwork)
  library(dplyr)
  library(tidyr)
  library(plyr)
  library(data.table)
})

lincs_connectivity <-
  function(input,
           num_features = 2000,
           method = c("LINCS", "Cor"),
           wtcs_fdr_cutoff = 1,
           trend_val = c("up", "down"),
           cor_score_cutoff = 0,
           output_dir,
           plots_dir,
           db_path = lincs) 
    # up and down genes
    cluster_upset <- input %>%
      dplyr::filter(diff_means > 0) %>%
      dplyr::arrange(desc(diff_means)) %>%
      dplyr::slice_head(n = num_features)
    cluster_upset = cluster_upset %>%
      plyr::rename(c('diff_means' = 'score')) %>%
      dplyr::select(protein_id, score)
    
    cluster_downset <- input %>%
      dplyr::filter(diff_means < 0) %>%
      dplyr::arrange(diff_means) %>%
      dplyr::slice_head(n = num_features)
    cluster_downset = cluster_downset %>%
      plyr::rename(c('diff_means' = 'score')) %>%
      dplyr::select(protein_id, score)
    
    
    if (method == "LINCS") {
      # map to ENTREZ identifiers
      upset = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = cluster_upset$protein_id,
        column = "ENTREZID",
        keytype = "SYMBOL"
      )
      downset = AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = cluster_downset$protein_id,
        column = "ENTREZID",
        keytype = "SYMBOL"
      )
      upset = as.character(upset)
      downset = as.character(downset)
      
      # LINCS-based similarity metric
      qSig_output <-
        signatureSearch::qSig(
          query = list(upset = upset, downset = downset),
          gess_method = "LINCS",
          refdb = lincs
        )
      qSig_output <-
        signatureSearch::gess_lincs(
          qSig = qSig_output,
          sortby = "NCS",
          tau = T,
          workers = 4
        )
      qSig_output <- result(qSig_output)
      
      # filter drugs
      qSig_output <- qSig_output %>%
        dplyr::filter(trend == trend_val &
                        WTCS_FDR < wtcs_fdr_cutoff) %>%
        arrange(trend, WTCS_FDR)
      drugs <- unique(qSig_output$pert)
    } else if (method == "Cor") {
      # correlation-based similarity
      query_mat <- cluster_upset %>%
        plyr::rbind.fill(cluster_downset) %>%
        mutate(id = mapIds(
          org.Hs.eg.db,
          keys = geneSymbol,
          column = "ENTREZID",
          keytype = "SYMBOL"
        )) %>%
        arrange(id) %>%
        filter(!is.na(id)) %>%
        column_to_rownames('id') %>%
        dplyr::select(score) %>% as.matrix()
      qSig_output <-
        qSig(query = query_mat,
             gess_method = "Cor",
             refdb = lincs)
      qSig_output <-
        gess_cor(qSig = qSig_output,
                 method = "spearman",
                 workers = 4)
      qSig_output <- result(qSig_output)
      
      # filter drugs
      qSig_output <- qSig_output %>%
        filter(cor_score < cor_score_cutoff)
      drugs <- unique(qSig_output$pert)
    }
    fname <-
      file.path(output_dir,  "qSig_output.txt")
    write.table(
      qSig_output,
      file = fname,
      quote = F,
      sep = "\t",
      row.names = F
    )
    
    # default plots are set to NULL and only updated where data is available
    p1 <- NULL
    if (length(drugs) != 0) {
      # Query signature drug barplot
      p1 <- drug_barplots(
        dat = qSig_output,
        xlab = "pert",
        ylab = "WTCS_FDR",
        top = 20,
        fill_var = NULL,
        title = "Query Signature"
      )
      # Save the drug barplot to PDF file
      plot_filename <-
        file.path(plots_dir, "drug_barplot.pdf")
      pdf(file = plot_filename,
          width = 8.5,
          height = 11)
      print(p1)
      dev.off()
    }
    return(qSig_output)
