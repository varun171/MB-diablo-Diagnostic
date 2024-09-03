# Function: DGE analysis by DESeq2

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(DESeq2)
})

# parse command line options
option_list <- list(
  make_option(c("--expr_mat"), type = "character",
              help = "expression data matrix, preferably counts (.rds) "),
  make_option(c("--gtf_file"), type = "character",
              help = "gencode gtf file"),
  make_option(c("--cluster_file"), type = "character",
              help = "path to cluster annotation file"),
  make_option(c("--sample_map_file"), type = "character",
              help = "path to sample mapping file"),
  make_option(c("--kegg_medicus_file"), type = "character",
              help = "path to KEGG MEDICUS pathway file (.gmt)"),
  make_option(c("--results_dir"), type = "character",
              help = "path to results directory"),
  make_option(c("--plots_dir"), type = "character",
              help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "multimodal_clustering", "dge_pathway_analysis")
utils_dir <- file.path(analysis_dir, "utils")
source(file.path(utils_dir, "perform_enrichment_gsea.R"))

# results directory
results_dir <- opt$results_dir
dir.create(results_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = opt$gtf_file) %>%
  as.data.frame() %>% 
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# count data
expr_mat <- readRDS(opt$expr_mat)

# get sample map information to assign sample ids
samples_map <- read_tsv(file.path(opt$sample_map_file))
expr_mat <- expr_mat %>%
  dplyr::select(samples_map$Kids_First_Biospecimen_ID_RNA)
colnames(expr_mat) <- samples_map$sample_id

# filter expression count file to contain only protein coding gene
expr_mat <- expr_mat %>%
  filter(rownames(expr_mat) %in% gencode_gtf$gene_name)

# now read cluster information for these samples
intnmf_clusters <- read_tsv(file.path(opt$cluster_file))
intnmf_clusters <- intnmf_clusters %>%
  dplyr::arrange(intnmf_cluster) %>%
  dplyr::mutate(intnmf_cluster = paste0("cluster_", intnmf_cluster)) 
expr_mat <- expr_mat %>%
  dplyr::select(intnmf_clusters$sample_id)
stopifnot(identical(colnames(expr_mat), intnmf_clusters$sample_id))

# use DCGA to filter out low count, low variance features
expr_mat <- DGCA::filterGenes(inputMat = expr_mat, 
                               filterTypes = c("central", "dispersion"),
                               filterDispersionType = "cv", 
                               filterDispersionPercentile = 0.2,
                               sequential = TRUE)

# use a for-loop 
clusters <- unique(intnmf_clusters$intnmf_cluster)
output_df <- data.frame()
for(i in 1:length(clusters)){
  
  # prefix for output files and plots
  prefix <- paste0(clusters[i], "_vs_rest")
  print(prefix)
  
  # cluster of interest
  intnmf_clusters$group <- ifelse(intnmf_clusters$intnmf_cluster == clusters[i], "COI", "Others")
  group <- as.factor(intnmf_clusters$group)
  
  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = round(expr_mat),
                                colData = intnmf_clusters,
                                design = ~group)
  dds <- DESeq(dds)
  deseq_results <- results(dds, contrast = c("group", "COI", "Others"))
  deseq_output <- deseq_results %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    dplyr::rename("log2FC" = "log2FoldChange") %>%
    mutate(direction = ifelse(log2FC > 0, "up", "down"),
           comparison = prefix) %>%
    rownames_to_column("genes") 
  
  # combine with dataframe
  output_df <- rbind(output_df, deseq_output)
  
  # pathway enrichment using REACTOME
  reactome_pathways <- msigdbr::msigdbr(category = "C2", subcategory = "CP:REACTOME")
  reactome_pathways <- reactome_pathways %>% 
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename("term" = "gs_name",
                  "gene" = "gene_symbol")
  perform_enrichment_gsea(diffexpr_res = deseq_output,
                          pathways = reactome_pathways,
                          minGSSize = 10,
                          maxGSSize = 150,
                          prefix = prefix,
                          plots_dir = file.path(plots_dir, "reactome"),
                          results_dir = file.path(results_dir, "reactome"))
  
  # pathway enrichment using KEGG MEDICUS
  kegg_medicus_pathways <- clusterProfiler::read.gmt(gmtfile = opt$kegg_medicus_file)
  perform_enrichment_gsea(diffexpr_res = deseq_output,
                          pathways = kegg_medicus_pathways,
                          minGSSize = 10,
                          maxGSSize = 500,
                          prefix = prefix,
                          plots_dir = file.path(plots_dir, "kegg_medicus"),
                          results_dir = file.path(results_dir, "kegg_medicus"))
  
}

# write output to tsv
write_tsv(x = output_df, file = file.path(results_dir, "diffexpr_output_per_cluster.tsv"))
