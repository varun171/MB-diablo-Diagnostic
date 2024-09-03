# Load required libraries
library(survival)
library(dplyr)

# Define the function to flatten nested lists
flatten_list <- function(lst) {
  data.frame(
    Cluster_ID = rep(names(lst), lengths(lst)),
    Gene = unlist(lst)
  )
}

# Set file paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering/data")
analysis_dir <- file.path(root_dir, "multimodal_clustering/diablo")
input_dir <- file.path(root_dir, "multimodal_clustering/intNMF/input")
nmf_dir <- file.path(root_dir, "multimodal_clustering/intNMF/results")
output_dir <- file.path(root_dir, "multimodal_clustering/diablo/results/intnmf/survival")
plots_dir <- file.path(root_dir, "multimodal_clustering/diablo/plots/intnmf/survival")

# Read data
count_data <- read_tsv(file.path(input_dir, "norm_counts.tsv")) %>% column_to_rownames() 
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(analysis_dir, "consensus_clusters.tsv")) %>% dplyr::rename("sample_id" = "samID")
mm_clusters_diablo <- mm_clusters %>% inner_join(samples_map)

# Subset cluster-specific samples and genes
tpm_data <- readRDS(file.path(input_dir, "gene-expression-rsem-tpm-collapsed-v12.rds"))
tpm_data <- tpm_data %>% dplyr::select(any_of(mm_clusters_diablo$Kids_First_Biospecimen_ID_RNA)) 
tpm_data <- tpm_data %>% filter(rownames(tpm_data) %in% colnames(count_data))
rownames(tpm_data) <- make.names(rownames(tpm_data))
# Combine new clusters with survival data
anno_file <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, EFS_days, age_last_update_days, age_at_diagnosis_days,molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters_diablo, by = "Kids_First_Biospecimen_ID_RNA")

# Process survival data
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
anno_file$OS_status <- ifelse(anno_file$OS_status == "DECEASED", 1, 0)

gene_of_interest <- colnames(t(tpm_data))
tpm_surv <- tpm_data %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Kids_First_Biospecimen_ID_RNA")
anno_df <- anno_file %>% dplyr::select(Kids_First_Biospecimen_ID_RNA, OS_status, OS_days, clust)
# Filter out incomplete cases
anno_df <- anno_df[complete.cases(anno_df), ]
###################################################################################################
# Extract gene names
# Initialize list to store results
cluster_specific_gene_signatures <- list()

# Loop through each cluster
for (cluster_id in unique(anno_df$clust)) {
  # Subset data for the current cluster
  merged_df_cluster <- merge(anno_df %>% filter(clust == cluster_id), tpm_surv, by = "Kids_First_Biospecimen_ID_RNA")
  
  # Get gene names for the current cluster (excluding the first four columns)
  gene_names <- colnames(merged_df_cluster)[-c(1:4)]
  
  # Remove unexpected symbols from gene names
  gene_names_cleaned <- gene_names
  
  # Update column names in the merged dataframe
  colnames(merged_df_cluster)[-c(1:4)] <- gene_names_cleaned
  
  # Perform Cox proportional hazards regression
  formula <- as.formula(paste("Surv(OS_days, OS_status) ~ ", paste(gene_names_cleaned, collapse = " + ")))
  cox_model <- try(coxph(formula, data = merged_df_cluster), silent = TRUE)
  
  # Check if Cox model was successfully created
  if (class(cox_model) == "try-error") {
    cat("Error: Unable to create Cox model for cluster", cluster_id, "\n")
  } else {
    # Store Cox results for the current cluster
    cox_results <- list(
      cox_model = cox_model,
      p_values = summary(cox_model)$coefficients[, "Pr(>|z|)"]
    )
    cluster_specific_gene_signatures[[paste("Cluster", cluster_id)]] <- cox_results
  }
}

# Print or further process cluster-specific gene signatures with significant p-values
print(cluster_specific_gene_signatures)
