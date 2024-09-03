suppressPackageStartupMessages({
  library(signatureSearch)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(plyr)
  library(data.table)
  library(ExperimentHub)
  library(rhdf5)
})

lincs_analysis <- function(input,
                           num_features = 2000,
                           wtcs_fdr_cutoff = 1,
                           trend_val = c("up", "down"),
                           cor_score_cutoff = 0,
                           output_dir,
                           plots_dir,
                           db_path = NULL) {
  
  upset <- input %>%
  arrange(log2FoldChange) %>%
  slice_head(n = num_features) %>%
  column_to_rownames(var = "V1")   
  
  downset <- input %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = num_features) %>%
  column_to_rownames(var = "V1")   
  
  upset <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rownames(upset),
    column = "ENTREZID",
    keytype = "SYMBOL"
  )
  downset <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rownames(downset),
    column = "ENTREZID",
    keytype = "SYMBOL"
  )
  
  genes <-
    list('upset' = na.omit(upset), 'downset' = na.omit(downset))
  
  qsig_lincs <- signatureSearch::qSig(query = genes,
                                      gess_method = "LINCS",
                                      refdb = lincs)
  
  qSig_output <-
    signatureSearch::gess_lincs(
      qSig = qsig_lincs,
      sortby = "NCS",
      tau = F,
      workers = 8,
      ref_trts = NULL,
      addAnnotations = F
    )
  
  qSig_output <- result(qSig_output)
  
  fname <-
    file.path(output_dir,  "qSig_output.txt")
  write.table(
    qSig_output,
    file = fname,
    quote = F,
    sep = "\t",
    row.names = F
  )
  
  # filter drugs
  qSig_output <- qSig_output %>%
    dplyr::filter(trend %in% trend_val &
                    WTCS_FDR < wtcs_fdr_cutoff) %>%
    arrange(trend, WTCS_FDR)
  drugs <- unique(qSig_output$pert)
  
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
}
