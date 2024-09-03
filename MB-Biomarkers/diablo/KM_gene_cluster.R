suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(datawizard)
  library(reshape2)
  library(survival)
  library(survminer)
  library(ggsankey)
})

data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo")
input_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/input")
nmf_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/results")
output_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/plots")

# read data
count_data <- read_tsv(file.path(input_dir, "norm_counts.tsv")) %>% column_to_rownames() 
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(nmf_dir, "intnmf_clusters.tsv")) 
mm_clusters_diablo <- mm_clusters %>% inner_join(samples_map)

# subset cluster specific samples and genes
tpm_data <- readRDS(file.path(input_dir, "gene-expression-rsem-tpm-collapsed-v12.rds"))
tpm_data <- tpm_data %>% dplyr::select(any_of(mm_clusters_diablo$Kids_First_Biospecimen_ID_RNA)) 
tpm_data <- tpm_data %>% filter(rownames(tpm_data) %in% colnames(count_data))
rownames(tpm_data) <- make.names(rownames(tpm_data))

# combine new clusters with survival data
anno_file <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, EFS_days, age_last_update_days, age_at_diagnosis_days,molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters_diablo, by = "Kids_First_Biospecimen_ID_RNA")

#colnames(anno_file_rna)

anno_file <- anno_file %>%
  mutate(OS_days = as.numeric(OS_days),
         PFS_days = as.numeric(EFS_days))
anno_file <- anno_file %>%
  mutate(EFS_status = case_when(OS_status == "LIVING" & EFS_days == OS_days ~ 1,
                                OS_status %in% c("LIVING", "DECEASED") & is.na(EFS_days) ~ 0,
                                OS_status %in% c("LIVING", "DECEASED") & EFS_days < OS_days ~ 1,
                                OS_status == "DECEASED" & EFS_days == OS_days ~ 0
  ))

anno_file <- anno_file %>%
  mutate(EFS_days = ifelse(is.na(EFS_days), age_last_update_days, EFS_days)) 

# Encode OS status
anno_file$OS_status <-
  ifelse(anno_file$OS_status == "DECEASED", 1, 0)


gene_of_interest <- colnames(t(tpm_data))
tpm_surv <- tpm_data %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Kids_First_Biospecimen_ID_RNA")
anno_df <- anno_file %>% dplyr::select(Kids_First_Biospecimen_ID_RNA, OS_status, OS_days, intnmf_cluster)
# Filter out incomplete cases
anno_df <- anno_df[complete.cases(anno_df), ]

# Now, proceed with your survival analysis using `anno_df`

# Create an empty list to store cluster-specific gene signatures
cluster_specific_gene_signatures <- list()

# Loop through each cluster
for (cluster_id in unique(anno_df$intnmf_cluster)) {
  
  # Subset annotation dataframe for the current cluster
  anno_df_cluster <- anno_df %>% filter(intnmf_cluster == cluster_id)
  
  # Merge expression data with annotation data for the current cluster
  merged_df_cluster <- merge(anno_df_cluster, tpm_surv, by = "Kids_First_Biospecimen_ID_RNA")
  
  # Create an empty list to store KM results for the current cluster
  km_results <- list()
  
  # Loop through each gene of interest
  for (gene in gene_of_interest) {
    
    # Check if the gene exists in the merged dataframe
    if (!(gene %in% colnames(merged_df_cluster))) {
      cat("Warning: Gene", gene, "not found in cluster", cluster_id, "\n")
      next  # Skip to the next gene
    }
    
    # Calculate median expression for the current gene
    median_expression <- median(merged_df_cluster[[gene]], na.rm = TRUE)
    
    # Create a new column indicating high or low expression based on median
    merged_df_cluster$expression_status <- ifelse(merged_df_cluster[[gene]] >= median_expression, "High", "Low")
    
    # Check if there are both High and Low expression statuses
    if (length(unique(merged_df_cluster$expression_status)) != 2) {
      cat("Warning: Unable to perform survival analysis for gene", gene, "due to insufficient data.\n")
      next  # Skip to the next gene
    }
    
    # Perform Kaplan-Meier analysis for the current gene in the cluster
    km_fit <- survfit(Surv(OS_days, OS_status) ~ expression_status, data = merged_df_cluster)
    
    # Perform log-rank test and get p-value
    p_value <- survdiff(Surv(OS_days, OS_status) ~ expression_status, data = merged_df_cluster)$chisq[1]
    
    # Store KM results for the current gene
    km_results[[gene]] <- list(
      km_fit = km_fit,
      p_value = p_value
    )
  }
  
  # Filter genes with significant p-value for the current cluster
  significant_genes <- Filter(function(x) x$p_value < 0.05, km_results)
  
  # Skip cluster if no significant genes found
  if (length(significant_genes) == 0) {
    cat("No significant genes found for cluster", cluster_id, "\n")
    next
  }
  
  # Save significant genes for the current cluster
  # Extract gene names
  significant_gene_names <- names(significant_genes)
  significant_gene_p_values <- unlist(lapply(significant_genes, function(x) x$p_value))
  # Create a data frame to store significant genes and p-values
  significant_genes_df <- data.frame(Cluster_ID = cluster_id, Gene = significant_gene_names, P_Value = significant_gene_p_values)
  
   # Save significant genes to a TSV file
  write.table(significant_genes_df, file = file.path(output_dir, paste("significant_genes_cluster_", cluster_id, ".tsv", sep = "")), sep = "\t", row.names = FALSE)
}

# Print or further process cluster-specific gene signatures with significant p-values
print(cluster_specific_gene_signatures)
