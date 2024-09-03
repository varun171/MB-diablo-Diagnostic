suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(datawizard)
  library(reshape2)
  library(survival)
  library(survminer)
  library(ggsankey)
  library(MOFA2)
  library(MOFAdata)
  library(tibble)
})

ot_dir <- file.path("~/Documents/v15-Pedcan")
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier")
result_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/plots")

# read data
# Read methyl cbtn data
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds")) 
methyl_cbtn_data <- as.data.frame(t(methyl_cbtn_data))

count_cbtn_data <- readRDS(file.path(data_dir, "gene-counts-rsem-expected_count_medullo-collapsed.rds")) 

# All Medulloblastoma samples
sample_cbtn_map <- read_tsv(file.path(data_dir, "histologies_medullo.tsv"))

# Filter the data where experimental_strategy is RNAseq and select specific columns
cbtn_rna_map <- sample_cbtn_map %>%
  filter(
    experimental_strategy == "RNA-Seq",
    composition == "Solid Tissue",
    !is.na(molecular_subtype)
  )%>%
  select(sample_id, Kids_First_Biospecimen_ID, molecular_subtype)

# Remove duplicates based on sample_id, keeping the first occurrence
cbtn_rna_map_unique <- cbtn_rna_map %>%
  distinct(sample_id, .keep_all = TRUE)

# Only methyl MB samples with predicted subtype information - use this to get the subtype for RNAseq using sample_id
sample_cbtn_methyl_map <- read_tsv(file.path(result_dir, "sample_cbtn_map.tsv"))

# Select specific columns
sample_cbtn_methyl_map_selected <- sample_cbtn_methyl_map %>%
  select(sample_id, Kids_First_Biospecimen_ID, molecular_subtype, rf_predicted_subtype)

#############################
# Find sample_id with more than one unique rf_predicted_subtype
sample_ids_with_multiple_rf_predicted_subtype <- sample_cbtn_methyl_map_selected %>%
  group_by(sample_id) %>%
  filter(n_distinct(rf_predicted_subtype) > 1) %>%
  distinct(sample_id)

# Display the sample_ids with more than one unique rf_predicted_subtype
print(sample_ids_with_multiple_rf_predicted_subtype)
#7316-8817 more than predicted subtype
########################

# Remove duplicates based on sample_id, keeping the first occurrence
sample_cbtn_methyl_map_unique <- sample_cbtn_methyl_map_selected %>%
  distinct(sample_id, .keep_all = TRUE)

# Join the two data frames on the sample_id column - keep the duplicate entries
merged_data <- cbtn_rna_map_unique %>%
  inner_join(sample_cbtn_methyl_map_unique %>% select(sample_id, rf_predicted_subtype), by = "sample_id")

# Filter count_cbtn_data to keep only the columns matching Kids_First_Biospecimen_Id in merged_data
filtered_count_data <- count_cbtn_data[, colnames(count_cbtn_data) %in% merged_data$Kids_First_Biospecimen_ID]

#######

# group 3 only for count data

# Filter the data where molecular_subtype is equal to "group3"
filtered_merged_data <- merged_data %>%
  filter(molecular_subtype == "MB, Group3" & !str_detect(rf_predicted_subtype, "4"))


# Filter count_cbtn_data to keep only the columns matching Kids_First_Biospecimen_Id in merged_data
filtered_count_data <- filtered_count_data[, colnames(filtered_count_data) %in% filtered_merged_data$Kids_First_Biospecimen_ID]


########

# CNV data for feature identifications -- did not work

# CNV data
cnv_data <- data.table::fread(file.path(ot_dir, "20231205_release.All.gainloss.txt.gz"))

cnv_dat <- cnv_data %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "BS_ID") %>%
  filter(Kids_First_Biospecimen_ID %in% sample_cbtn_map$Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Biospecimen_ID, gene, log2) %>%
  inner_join(
    sample_cbtn_map %>% dplyr::select(
      Kids_First_Biospecimen_ID,
      sample_id,
      tumor_ploidy,
      tumor_fraction
    ),
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  unique()

# get unique biospecimens per sample id
unique_ids <- cnv_dat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
cnv_dat <- cnv_dat %>%
  filter(
    sample_id %in% unique_ids$sample_id,
    Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID
  )
cnv_samples <- cnv_dat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
# function to adjust copy number and status
adjust_cn <- function(x) {
  # get tumor fraction and ploidy
  tumor_fraction <- unique(na.omit(x$tumor_fraction))
  tumor_ploidy <- unique(na.omit(x$tumor_ploidy))
  
  if (length(tumor_fraction) == 1 & length(tumor_ploidy) == 1) {
    # calculate adjusted copy number if tumor fraction and ploidy info is available
    x$adjusted_cn <- (((2 ^ (x$log2) - (
      1 - tumor_fraction
    )) * tumor_ploidy) / tumor_fraction) - 0.5
    x$adjusted_cn <- round(x$adjusted_cn)
    x$adjusted_status <- ifelse(x$adjusted_cn == 0,
                                "Complete Loss",
                                ifelse(
                                  x$adjusted_cn == 1,
                                  "Loss",
                                  ifelse(
                                    x$adjusted_cn %in% c(tumor_ploidy + 1:9),
                                    "Gain",
                                    ifelse(x$adjusted_cn >= 10, "Amplification", "Neutral")
                                  )
                                ))
    
    # replace old columns with new ones
    x <- x %>%
      dplyr::rename("status" = "adjusted_status", # rename new columns
                    "copy_number" = "adjusted_cn")
    
  }
}

# apply function to all samples in the consensus file
cnv_dat <- plyr::ddply(
  .data = cnv_dat,
  .variables = "sample_id",
  .fun = function(x)
    adjust_cn(x = x)
)
cnv_dat <- cnv_dat %>%
  filter(status != "Neutral") %>%
  dplyr::select(-c(log2)) %>%
  unique()


# source functions
source(file.path(data_dir, "filter_cnv.R"))

# cancer genes
cancer_genes <- readRDS(file.path(ot_dir, "cancer_gene_list.rds"))

# subset to cancer genes
cnv_dat <- cnv_dat %>%
  dplyr::rename("hgnc_symbol" = "gene") %>%
  filter_cnv(myCancerGenes = cancer_genes)
cnv_dat <- acast(
  cnv_dat,
  sample_id ~ hgnc_symbol,
  value.var = "copy_number",
  fill = 0,
  fun.aggregate = max
)

# print dimensions
print(dim(cnv_dat))

# Join the two data frames on the sample_id column - keep the duplicate entries
sample_cbtn_cnv_map <- sample_cbtn_map %>% 
                          filter(sample_id %in% rownames(cnv_dat)) %>%
  select(sample_id, Kids_First_Biospecimen_ID, molecular_subtype) %>%
  distinct(sample_id, .keep_all = TRUE)
  

# Join the two data frames on the sample_id column - keep the duplicate entries
merged_cnv_map <- sample_cbtn_cnv_map %>%
  inner_join(sample_cbtn_methyl_map_unique %>% select(sample_id, rf_predicted_subtype), by = "sample_id")

# Convert to data frame for filtering if it's a matrix
cnv_dat <- as.data.frame(cnv_dat)

# Filter the data frame based on row names
cnv_dat <- cnv_dat %>% filter(rownames(cnv_dat) %in% merged_cnv_data$sample_id)

# Convert back to matrix if needed
cnv_dat <- as.matrix(cnv_dat)

# MOFA - did not work

# Function to train MOFA and plot feature weights
train_and_plot_mofa <- function(subtype, subtype_data, plots_dir) {
  # Transpose the data for MOFA
  single_omics_data <- list(omics_data = t(subtype_data))
  
  # Create MOFA object
  MOFAobject <- create_mofa(single_omics_data)
  
  # Define data options
  data_opts <- get_default_data_options(MOFAobject)
  
  # Define model options
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- 5
  
  # Define training options
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- "slow"
  train_opts$seed <- 40
  
  # Prepare and train the MOFA model
  MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts)
  MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)
  
  # Get the number of active factors
  num_active_factors <- MOFAobject@dimensions$K
  
  # Plotting the feature weights for each active factor to identify important features
  for (factor in 1:num_active_factors) {
    plot_title <- paste("Top Features for Factor", factor, "in Subtype", subtype)
    print(paste("Plotting feature weights for Factor", factor, "in Subtype", subtype))
    
    p <- plot_weights(MOFAobject,
                      view = "omics_data",
                      factor = factor,
                      nfeatures = 10,   # Number of top features to plot
                      scale = TRUE) +
      ggtitle(plot_title)
    
    # Save the plot
    ggsave(filename = file.path(plots_dir, paste0("Subtype_", subtype, "_Factor_", factor, "_Top_Features.pdf")), plot = p)
  }
}

# Filter data by each subtype and train MOFA model
for (subtype in unique(merged_cnv_map$rf_predicted_subtype)) {
  subtype_data <- cnv_dat[rownames(cnv_dat) %in% merged_cnv_map$sample_id[merged_cnv_map$rf_predicted_subtype == subtype], ]
  
  print(paste("Training MOFA model for Subtype:", subtype))
  
  # Train MOFA and plot feature weights for the current subtype
  train_and_plot_mofa(subtype, subtype_data, plots_dir)
}

##################################################################################################

# direct copy number using adjusted copy number for subtype validation - not 100% sure

# Required Libraries

# Convert to data frame for filtering if it's a matrix
cnv_dat <- as.data.frame(t(cnv_dat))

# Step 1: Filter cnv_dat for the Genes of Interest
genes_of_interest <- c("MYC", "MYCN", "OTX2")
cnv_filtered <- cnv_dat[rownames(cnv_dat) %in% genes_of_interest, ]

# Convert to data frame for easier manipulation
cnv_filtered_df <- as.data.frame(cnv_filtered)

# Add gene names as a column for merging later
cnv_filtered_df$gene <- rownames(cnv_filtered_df)

# Step 2: Add rf_predicted_subtype to the Filtered Data
# Ensure sample_id is a column in merged_cnv_map
merged_cnv_map <- merged_cnv_map %>% select(sample_id, rf_predicted_subtype)

# Gather the cnv_filtered_df to long format
cnv_long <- cnv_filtered_df %>%
  gather(key = "sample_id", value = "copy_number", -gene)

# Merge with subtype information
cnv_long <- cnv_long %>%
  inner_join(merged_cnv_map, by = "sample_id")

# Step 3: Determine Status Based on Copy Number
calculate_status <- function(copy_number, tumor_ploidy) {
  adjusted_status <- ifelse(copy_number == 0, "Complete Loss",
                            ifelse(copy_number == 1, "Loss",
                                   ifelse(copy_number %in% (tumor_ploidy + 1):(tumor_ploidy + 9), "Gain",
                                          ifelse(copy_number >= 10, "Amplification", "Neutral"))))
  return(adjusted_status)
}

tumor_ploidy <- 2  # Example value, replace with actual value if available

cnv_long <- cnv_long %>%
  mutate(adjusted_status = calculate_status(copy_number, tumor_ploidy))

# Step 4: Summarize by Subtype
summary_by_subtype <- cnv_long %>%
  group_by(rf_predicted_subtype, gene, adjusted_status) %>%
  summarise(count = n()) %>%
  spread(key = adjusted_status, value = count, fill = 0) %>%
  arrange(rf_predicted_subtype, gene)

# Display the summary
print(summary_by_subtype)

 #################################

# direct copy number using log2 values for subtype validation

# Step 1: Filter cnv_dat for the Genes of Interest
genes_of_interest <- c("MYC", "MYCN", "OTX2")
cnv_filtered <- cnv_dat %>% filter(gene %in% genes_of_interest)

# Step 2: Use merged_cnv_map to filter the sample_ids
# Ensure sample_id is a column in merged_cnv_map
merged_cnv_map <- merged_cnv_map %>% select(sample_id, rf_predicted_subtype)

# Merge cnv_filtered with merged_cnv_map to get the rf_predicted_subtype for each sample
cnv_with_subtypes <- cnv_filtered %>%
  inner_join(merged_cnv_map, by = "sample_id")

# Step 3: Calculate the average log2 for each predicted subtype for each gene
average_log2_by_subtype <- cnv_with_subtypes %>%
  group_by(rf_predicted_subtype, gene) %>%
  summarise(average_log2 = mean(log2, na.rm = TRUE)) %>%
  arrange(rf_predicted_subtype, gene)

# Display the summary
print(average_log2_by_subtype)
 ######################################################

# differential expression using deseq2

# Load Required Libraries
library(dplyr)
library(tidyr)
library(readr)
library(DESeq2)
library(Biobase)
library(rtracklayer)
library(tibble)

ot_dir <- file.path("~/Documents/v15-Pedcan")
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier")
result_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/plots")


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
  tumor_count_data %>%
    select(any_of(tumor_hist_df$Kids_First_Biospecimen_ID)),
  gtex_count_data %>%
    select(any_of(gtex_hist_df$Kids_First_Biospecimen_ID))
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

# Differential expression analysis
unique_histologies <- unique(final_hist_df$rf_predicted_subtype)
unique_histologies <- unique_histologies[unique_histologies != "Normal_Brain"]

for (histology in unique_histologies) {
  final_hist_df_subset <- final_hist_df %>%
    filter(rf_predicted_subtype %in% c(histology, "Normal_Brain"))
  
  expr_data_subset <- combined_expr_data %>%
    select(any_of(final_hist_df_subset$Kids_First_Biospecimen_ID))
  
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
  output_file <- file.path(result_dir, paste0("DESeq2_results_", histology, "_vs_Normal_Brain.csv"))
  write.csv(as.data.frame(res), file = output_file)
  
  # Print top differentially expressed genes
  print(paste("Top differentially expressed genes for", histology, "vs Normal_Brain:"))
  print(head(res))
}

grp3_alpha_DE <- read_csv(file.path(result_dir, "DESeq2_results_subtype: Group3_alpha_vs_Normal_Brain.csv"))
grp3_beta_DE <- read_csv(file.path(result_dir, "DESeq2_results_subtype: Group3_beta_vs_Normal_Brain.csv"))
grp3_gamma_DE <- read_csv(file.path(result_dir, "DESeq2_results_subtype: Group3_gamma_vs_Normal_Brain.csv"))
