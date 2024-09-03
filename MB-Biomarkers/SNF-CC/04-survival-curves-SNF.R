# create survival curves for IntNMF clusters, RNA and methyl-derived subtypes 

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC","data")
analysis_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results")

# plots directory
plots_dir <- file.path(analysis_dir, "plots", "survival")
dir.create(plots_dir, showWarnings = F, recursive = T)

# combine samples map with IntNMF derived clusters
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
samples_map <- samples_map %>%
  arrange(sample_id)

mm_clusters <- read_tsv(file.path(output_dir, "SNF_clusters.tsv"), col_names = FALSE)
samples_map <- samples_map %>%
  mutate(clusterID = mm_clusters$X1)

# combine SNFCC clusters with RNA-derived molecular subtypes
anno_file_rna <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(samples_map, by = "Kids_First_Biospecimen_ID_RNA")

# combine SNFCC clusters with methylation-derived subclass
anno_file_methyl <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, dkfz_v11_methylation_subclass, dkfz_v12_methylation_subclass) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(samples_map, by = "Kids_First_Biospecimen_ID_Methyl")

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  inner_join(anno_file_methyl)
anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)

# format survival data
surv_data <- anno_file
surv_data$OS_status <- ifelse(surv_data$OS_status == "DECEASED", 1, 0)

# survival curves stratified by RNA-derived molecular subtype 
surv_data$molecular_subtype <- factor(surv_data$molecular_subtype, levels = sort(unique(surv_data$molecular_subtype)))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ molecular_subtype, data = surv_data)
pdf(file = file.path(plots_dir, "MB", "survival_molsubtype.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "Survival stratified by RNA-derived molecular subtypes",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.coord = c(6000, 0.75),
                pval.method = TRUE,
                pval.method.coord = c(6000, 0.80),
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

# survival curves stratified by Methylation-derived molecular subtype (v11)
surv_data$dkfz_v11_methylation_subclass <- factor(surv_data$dkfz_v11_methylation_subclass, levels = sort(unique(surv_data$dkfz_v11_methylation_subclass)))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ dkfz_v11_methylation_subclass, data = surv_data)
pdf(file = file.path(plots_dir,"MB", "survival_methyl_subtype_v11.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "Survival stratified by Methylation-derived molecular subtypes",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.coord = c(6000, 0.75),
                pval.method = TRUE,
                pval.method.coord = c(6000, 0.80),
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

# survival curves stratified by Methylation-derived molecular subtype (v12)
surv_data$dkfz_v12_methylation_subclass <- factor(surv_data$dkfz_v12_methylation_subclass, levels = sort(unique(surv_data$dkfz_v12_methylation_subclass)))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ dkfz_v12_methylation_subclass, data = surv_data)
pdf(file = file.path(plots_dir, "MB", "survival_methyl_subtype_v12.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "Survival stratified by Methylation-derived molecular subtypes",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.coord = c(6000, 0.75),
                pval.method = TRUE,
                pval.method.coord = c(6000, 0.80),
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

# survival curves stratified by SNFCC subtypes
surv_data$clusterID <- factor(surv_data$clusterID, levels = sort(as.numeric(unique(surv_data$clusterID))))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ clusterID, data = surv_data)
pdf(file = file.path(plots_dir, "MB", "survival_Snfcc_clusters.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "Survival stratified by Newly derived clusters",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.coord = c(6000, 0.75),
                pval.method = TRUE,
                pval.method.coord = c(6000, 0.80),
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  +
  guides(colour = guide_legend(ncol = 4)) 
print(p)
dev.off()
