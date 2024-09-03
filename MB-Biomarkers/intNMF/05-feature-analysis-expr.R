# script to perform downstream feature level analysis for expression data

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "multimodal_clustering")
utils_dir <- file.path(analysis_dir, "utils")

# source functions
source(file.path(utils_dir, "fgsea_barplot.R"))

# input directory
input_dir <- file.path(analysis_dir, "results")

# output directory
output_dir <- file.path(analysis_dir, "results", "feature_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots", "feature_analysis")
dir.create(plots_dir, showWarnings = F, recursive = T)

# genesets 
gene_set <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") 
gene_set <- gene_set %>% dplyr::select(gs_name, human_gene_symbol)
gene_set <- base::split(gene_set$human_gene_symbol, list(gene_set$gs_name))

# read full output
nmf_output <- readRDS(file.path(input_dir, "intnmf_best_fit.rds"))

# read expression and extract features with non-zero NMF weights per cluster
expr_selected <- nmf_output$H$H1
expr_selected <- expr_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  filter(value > 0) %>%
  mutate(cluster = as.numeric(cluster))

# expression data: fgsea barplot
fgsea_barplot(input_dat = expr_selected, 
              gene_set = gene_set, 
              output_file = file.path(plots_dir, "expression_fgsea.pdf"))
