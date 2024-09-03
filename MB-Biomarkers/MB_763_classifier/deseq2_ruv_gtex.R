# Load Required Libraries
library(dplyr)
library(tidyr)
library(readr)
library(DESeq2)
library(Biobase)
library(rtracklayer)
library(tibble)
library(RUVSeq)
library(ggplot2)
library(ggfortify)

# Directories
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
result_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/plots")

# Ensure output directories exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Step 1: Load the data
count_cbtn_data <- readRDS(file.path(data_dir, "gene-counts-rsem-expected_count_medullo-collapsed.rds"))

# All Medulloblastoma samples
sample_cbtn_map <- read_tsv(file.path(data_dir, "histologies_medullo.tsv"))

# Filter the data where experimental_strategy is RNA-Seq and select specific columns
cbtn_rna_map <- sample_cbtn_map %>%
  filter(
    experimental_strategy == "RNA-Seq",
    composition == "Solid Tissue",
    !is.na(molecular_subtype)
  ) %>%
  select(sample_id, Kids_First_Biospecimen_ID, molecular_subtype)

# Remove duplicates based on sample_id, keeping the first occurrence
cbtn_rna_map_unique <- cbtn_rna_map %>%
  distinct(sample_id, .keep_all = TRUE)

# Only methyl MB samples with predicted subtype information - use this to get the subtype for RNA-Seq using sample_id
sample_cbtn_methyl_map <- read_tsv(file.path(result_dir, "sample_cbtn_map.tsv"))

# Select specific columns
sample_cbtn_methyl_map_selected <- sample_cbtn_methyl_map %>%
  select(sample_id, Kids_First_Biospecimen_ID, molecular_subtype, rf_predicted_subtype)

# Remove duplicates based on sample_id, keeping the first occurrence
sample_cbtn_methyl_map_unique <- sample_cbtn_methyl_map_selected %>%
  distinct(sample_id, .keep_all = TRUE)

# Join the two data frames on the sample_id column
tumor_hist_df <- cbtn_rna_map_unique %>%
  inner_join(sample_cbtn_methyl_map_unique %>% select(sample_id, rf_predicted_subtype), by = "sample_id")

# Filter count_cbtn_data to keep only the columns matching Kids_First_Biospecimen_Id in merged_data
tumor_count_data <- count_cbtn_data[, colnames(count_cbtn_data) %in% tumor_hist_df$Kids_First_Biospecimen_ID]

# Select 3 samples from each brain region for the normal brain data
set.seed(123)  # For reproducibility
gtex_hist_df <- file.path(data_dir, "histologies.tsv") %>%
  read_tsv() %>%
  filter(
    grepl("GTEx", cohort),
    grepl("Brain", gtex_group),
    experimental_strategy == "RNA-Seq"
  ) %>%
  mutate(molecular_subtype = "Normal_Brain") %>%
  mutate(short_histology = "Normal_Brain") %>%
  group_by(gtex_subgroup) %>%
  sample_n(size = min(3, n()), replace = FALSE) %>%
  ungroup() %>%
  select(sample_id, Kids_First_Biospecimen_ID, short_histology, molecular_subtype) %>%
  dplyr::rename(rf_predicted_subtype = short_histology)

# Combine histology data
final_hist_df <- bind_rows(tumor_hist_df, gtex_hist_df)

# Read gene expression data for GTEx
gtex_count_data <- readRDS(file.path(data_dir, "gtex_gene-counts-rsem-expected_count-collapsed.rds"))

# Filter gtex_count_data to include only the columns that match Kids_First_Biospecimen_ID in gtex_hist_df
gtex_count_data <- gtex_count_data[, colnames(gtex_count_data) %in% gtex_hist_df$Kids_First_Biospecimen_ID]

# Identify common genes
common_genes <- intersect(rownames(tumor_count_data), rownames(gtex_count_data))

# Subset both datasets to include only common genes
tumor_count_data <- tumor_count_data[common_genes, ]
gtex_count_data <- gtex_count_data[common_genes, ]

# Combine expression data
combined_expr_data <- cbind(
  tumor_count_data[, colnames(tumor_count_data) %in% tumor_hist_df$Kids_First_Biospecimen_ID],
  gtex_count_data[, colnames(gtex_count_data) %in% gtex_hist_df$Kids_First_Biospecimen_ID]
)

# Round combined expression data to nearest integer
combined_expr_data <- round(combined_expr_data)

# Filter for protein-coding genes
gencode_gtf <- rtracklayer::import(con = file.path(data_dir, "gencode.v39.primary_assembly.annotation.gtf.gz")) %>%
  as.data.frame() %>%
  select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

combined_expr_data <- combined_expr_data %>%
  rownames_to_column('gene') %>%
  filter(gene %in% gencode_gtf$gene_name) %>%
  column_to_rownames('gene')

# Apply RUVSeq for removing unwanted variation only to GTEx (normal brain) data
gtex_set <- newSeqExpressionSet(as.matrix(gtex_count_data),
                                phenoData = data.frame(gtex_hist_df, row.names = colnames(gtex_count_data)))

# Function to plot PCA
plot_pca <- function(data, metadata, title, color_by) {
  pca_res <- prcomp(t(data))
  pca_df <- data.frame(pca_res$x)
  pca_df$color_by <- as.factor(metadata[[color_by]])
  ggplot(pca_df, aes(x = PC1, y = PC2, color = color_by)) +
    geom_point() +
    ggtitle(title) +
    theme_minimal() +
    labs(color = color_by)
}

# Plot PCA before RUV normalization
pca_plot_before <- plot_pca(combined_expr_data, final_hist_df, "PCA Before RUV", "rf_predicted_subtype")
print(pca_plot_before)
ggsave(file.path(plots_dir, "PCA_Before_RUV.pdf"), plot = pca_plot_before)

# Explore different values of k for RUVg normalization and plot PCA
for (k_val in 1:5) {
  set_k <- RUVg(gtex_set, rownames(gtex_count_data), k = k_val)
  normalized_counts_gtex_k <- exprs(set_k)
  
  # Combine normalized GTEx data with tumor data
  combined_normalized_counts_k <- cbind(
    tumor_count_data,
    normalized_counts_gtex_k
  )
  
  # Round combined normalized counts to nearest integer
  combined_normalized_counts_k <- round(combined_normalized_counts_k)
  
  # Plot PCA after RUV normalization
  pca_plot_after_ruv <- plot_pca(combined_normalized_counts_k, final_hist_df, paste("PCA After RUV (k =", k_val, ")"), "rf_predicted_subtype")
  print(pca_plot_after_ruv)
  ggsave(file.path(plots_dir, paste0("PCA_After_RUV_gtex_k_", k_val, ".pdf")), plot = pca_plot_after_ruv)
}

# Apply chosen RUVg normalization (example with k=1)
set_gtex <- RUVg(gtex_set, rownames(gtex_count_data), k = 1)
normalized_counts_gtex <- exprs(set_gtex)

# Combine normalized GTEx data with tumor data
combined_normalized_counts <- cbind(
  tumor_count_data,
  normalized_counts_gtex
)

# Round combined normalized counts to nearest integer
combined_normalized_counts <- round(combined_normalized_counts)

# Plot PCA after chosen RUV normalization (k=1)
pca_plot_after <- plot_pca(combined_normalized_counts, final_hist_df, "PCA After RUV (k=1)", "rf_predicted_subtype")
print(pca_plot_after)
ggsave(file.path(plots_dir, "PCA_After_RUV_k_1.pdf"), plot = pca_plot_after)

# Differential expression analysis using normalized counts
unique_histologies <- unique(final_hist_df$rf_predicted_subtype)
unique_histologies <- unique_histologies[unique_histologies != "Normal_Brain"]

for (histology in unique_histologies) {
  final_hist_df_subset <- final_hist_df %>%
    filter(rf_predicted_subtype %in% c(histology, "Normal_Brain"))
  
  expr_data_subset <- combined_normalized_counts[, final_hist_df_subset$Kids_First_Biospecimen_ID]
  
  # Ensure sample names match between expr_data_subset and final_hist_df_subset
  expr_data_subset <- expr_data_subset[, order(colnames(expr_data_subset))]
  final_hist_df_subset <- final_hist_df_subset[order(final_hist_df_subset$Kids_First_Biospecimen_ID), ]
  
  # Ensure rf_predicted_subtype is a factor with unique levels
  unique_subtypes <- unique(final_hist_df_subset$rf_predicted_subtype)
  unique_subtypes <- unique_subtypes[unique_subtypes != "Normal_Brain"]
  final_hist_df_subset <- final_hist_df_subset %>%
    mutate(rf_predicted_subtype = factor(rf_predicted_subtype, levels = c("Normal_Brain", unique_subtypes))) %>%
    column_to_rownames(var = "Kids_First_Biospecimen_ID")
  
  # Create a DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = expr_data_subset,
                                colData = final_hist_df_subset,
                                design = ~ rf_predicted_subtype)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results for the specific histology vs Normal_Brain
  res <- results(dds, contrast = c("rf_predicted_subtype", histology, "Normal_Brain"))
  
  # Order results by adjusted p-value
  res <- res[order(res$padj), ]
  
  # Save results to file
  output_file <- file.path(result_dir, paste0("DESeq2_ruv_results_", histology, "_vs_Normal_Brain.csv"))
  write.csv(as.data.frame(res), file = output_file)
  
  # Print top differentially expressed genes
  print(paste("Top differentially expressed genes for", histology, "vs Normal_Brain:"))
  print(head(res))
}

grp3_alpha_DE <- read_csv(file.path(result_dir, "DESeq2_ruv_results_subtype: Group3_alpha_vs_Normal_Brain.csv"))
grp3_beta_DE <- read_csv(file.path(result_dir, "DESeq2_ruv_results_subtype: Group3_beta_vs_Normal_Brain.csv"))
grp3_gamma_DE <- read_csv(file.path(result_dir, "DESeq2_ruv_results_subtype: Group3_gamma_vs_Normal_Brain.csv"))
