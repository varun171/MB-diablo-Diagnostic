# run multi-modal clustering with non-SNFCC
# Script adapted from @KomalRathi MB MM-clustering using IntNMF
suppressPackageStartupMessages({
  library(SNFtool)
  library(tidyverse)
  library(CancerSubtypes)
  library(scater)
  library(scran)
})

# define directories
ot_dir <- file.path("~/Documents", "GitHub", "OpenPedCan-analysis", "data", "v12")
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC", "data")
analysis_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir,"results","UMAP")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = F, recursive = T)

# read data
count_data <- read_tsv(file.path(input_dir, "norm_counts1.tsv")) %>% column_to_rownames() 
methyl_data <- read_tsv(file.path(input_dir, "methyl_data1.tsv")) %>% column_to_rownames() 
#snv_data <- read_tsv(file.path(input_dir, "snv_data.tsv")) %>% column_to_rownames() %>% t()
#cnv_data <- read_tsv(file.path(input_dir, "cnv_data.tsv")) %>% column_to_rownames()%>% t()
#splice_data <- read_tsv(file.path(input_dir, "splice_data1.tsv")) %>% column_to_rownames()%>% t()
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv")) 


# combine into a list (took almost 6 hours)
#dat <- list(as.matrix(count_data), as.matrix(methyl_data), as.matrix(snv_data), as.matrix(cnv_data), as.matrix(splice_data))
dat <- list(as.matrix(count_data), as.matrix(methyl_data))

#Practice
count_data <- read_tsv(file.path(input_dir, "norm_counts1.tsv")) %>% column_to_rownames() 

dat <- runMultiUMAP(dat, metric = "euclidean", name = "MultiUMAP")
wt = if(is.list(dat)) rep(1,length(dat)) else 1
# source function
source(file.path(analysis_dir, "utils", "dbscan_clustering.R"))
source(file.path(analysis_dir, "utils", "intnmf_clustering.R"))
source(file.path(analysis_dir, "utils", "lspline_clustering.R"))
source(file.path(analysis_dir, "utils", "run_ccp.R"))
source(file.path(analysis_dir, "utils", "get_cdf_datapoints.R"))
source(file.path(analysis_dir, "utils", "nbmclust_clustering.R"))
source(file.path(analysis_dir, "utils", "final_composite_score.R"))

dbsacn_output <- run_dbscan(dat = dat, 
                               wt = wt, 
                               output_dir = output_dir, 
                               k_value = 15)

expr_mat <- read_tsv(file.path(input_dir, "norm_counts1.tsv")) %>% column_to_rownames()%>% t()


# dbscan: Issue with kneedle package not avaliable for this version of R (t worked)
dbsacn_output <- dbscan_clustering (expr_mat = expr_mat,
                              filter_expr = TRUE, dispersion_percentile_val = 0.2,
                              protein_coding_only = TRUE, gencode_version = 27,
                              feature_selection = "variance",
                              var_prop = 10, transformation_type = "log2",
                              eps = NULL, minpts_val = NULL, output_dir)

expr_mat <- read_tsv(file.path(input_dir, "norm_counts1.tsv"), locale = locale(encoding = "UTF-8")) %>%
  column_to_rownames() %>% t()

# IntNMF
IntNMF_output <- intnmf_clustering (expr_mat = expr_mat,
                                    filter_expr = TRUE, dispersion_percentile_val = 0.2,
                                    protein_coding_only = TRUE, gencode_version = 27,
                                    feature_selection = "variance",
                                    var_prop = 10, transformation_type = "log2",
                                    max_k = 6, output_dir)


expr_mat <- read_tsv(file.path(input_dir, "norm_counts1.tsv"), locale = locale(encoding = "UTF-8")) %>%
  column_to_rownames() %>% t()

#lspline # No hist and No Cox survival - using sample ID and not specimen ID
lspline_output <- lspline_clustering (expr_mat, hist_file, 
                               algorithms = c("hc", "pam", "km"), 
                               distances = c("pearson", "spearman", "euclidean", "manhattan", "binary", "maximum", "canberra", "minkowski"),
                               filter_expr = TRUE, dispersion_percentile_val = 0.2,
                               protein_coding_only = TRUE, gencode_version = 27,
                               feature_selection = "variance", 
                               var_prop = 10, min_n = NULL, transformation_type = "log2",
                               max_k = 6, coef_cutoff = 0.5,
                               min_cluster_size_prop = NULL, max_cluster_size_prop = NULL,
                               compute_all_equal = TRUE, output_dir)

expr_mat <- read_tsv(file.path(input_dir, "norm_counts1.tsv"), locale = locale(encoding = "UTF-8")) %>%
  column_to_rownames() %>% t()

# nbmclust
nbmclust_output <- nbmclust_clustering (expr_mat, filter_expr = TRUE, dispersion_percentile_val = 0.2,
                                        protein_coding_only = TRUE, gencode_version = 27,
                                        feature_selection = "variance",
                                        var_prop = 10,
                                        max_k = 10,
                                        output_dir)

