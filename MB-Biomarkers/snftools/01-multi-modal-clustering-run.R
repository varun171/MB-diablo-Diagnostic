# run multi-modal clustering
suppressPackageStartupMessages({
  library(SNFtool)
  library(tidyverse)
  library(survival)
  library(survminer)
})

# define directories
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/snf")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = F, recursive = T)


# read data
count_data <- read_tsv(file.path(input_dir, "norm_counts.tsv")) %>% column_to_rownames()
methyl_data <- read_tsv(file.path(input_dir, "methyl_data.tsv")) %>% column_to_rownames()
snv_data <- read_tsv(file.path(input_dir, "snv_data.tsv")) %>% column_to_rownames()
cnv_data <- read_tsv(file.path(input_dir, "cnv_data.tsv")) %>% column_to_rownames()
splice_data <- read_tsv(file.path(input_dir, "splice_data.tsv")) %>% column_to_rownames()

# data Normalization
count_data = standardNormalization(count_data)
methyl_data = standardNormalization(methyl_data)
snv_data = standardNormalization(snv_data)
cnv_data = standardNormalization(cnv_data)
splice_data = standardNormalization(splice_data)

## Calculate the pair-wise distance;
## If the data is continuous, we recommend to use the function "dist2" as follows
count_data = (dist2(as.matrix(count_data),as.matrix(count_data)))^(1/2)
methyl_data = (dist2(as.matrix(methyl_data),as.matrix(methyl_data)))^(1/2)
splice_data = (dist2(as.matrix(splice_data),as.matrix(splice_data)))^(1/2)
cnv_data = (dist2(as.matrix(cnv_data),as.matrix(cnv_data)))^(1/2)
snv_data = (dist2(as.matrix(snv_data),as.matrix(snv_data)))^(1/2)

# Define the range for K, alpha, and T
K_range <- seq(10, 30, by = 10)
alpha_range <- seq(0.3, 0.8, by = 0.1)
T_range <- seq(10, 20, by = 10)

# Initialize a list to store estimation results for each combination of parameters
results <- list()

# Nested loops to iterate over K, alpha, and T
for (K in K_range) {
  for (alpha in alpha_range) {
    for (T in T_range) {
      # Compute affinity matrices for different data sources
      W1 <- affinityMatrix(count_data, K, alpha)
      W2 <- affinityMatrix(methyl_data, K, alpha)
      W3 <- affinityMatrix(snv_data, K, alpha)
      W4 <- affinityMatrix(cnv_data, K, alpha)
      W5 <- affinityMatrix(splice_data, K, alpha)
      
      # Fuse affinity matrices using SNF
      W <- SNF(list(W1, W2, W3, W4, W5), K, T)
      
      # Estimate the number of clusters
      estimationResult <- estimateNumberOfClustersGivenGraph(W, 2:8)
      
      # Store the results along with the parameter combination
      result_entry <- list(K = K, alpha = alpha, T = T, estimationResult = estimationResult)
      results <- c(results, list(result_entry))
    }
  }
}

results_df <- do.call(rbind, lapply(results, as.data.frame))

# Write the results to a TSV file
write.table(results_df, file.path(output_dir, "results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

################################################################
# for 4 clusters
K = 10; alpha = 0.3; T = 10; C = 4;
W1 = affinityMatrix(count_data, K, alpha)
W2 = affinityMatrix(methyl_data, K, alpha)
W3 = affinityMatrix(snv_data, K, alpha)
W4 = affinityMatrix(cnv_data, K, alpha)
W5 = affinityMatrix(splice_data, K, alpha)
W = SNF(list(W1,W2,W3,W4,W5), K, T)
group = spectralClustering(W,C)

snf_clusters <- data.frame(sample_id = rownames(count_data), snf_cluster = group)

NMI_scores <- rankFeaturesByNMI(list(count_data, methyl_data, snv_data, cnv_data, splice_data), W)

#######################################################################
# for 5 clusters
K = 10; alpha = 0.5; T = 10; C = 5;
W1 = affinityMatrix(count_data, K, alpha)
W2 = affinityMatrix(methyl_data, K, alpha)
W3 = affinityMatrix(snv_data, K, alpha)
W4 = affinityMatrix(cnv_data, K, alpha)
W5 = affinityMatrix(splice_data, K, alpha)
W = SNF(list(W1,W2,W3,W4,W5), K, T)
group = spectralClustering(W,C)

snf_clusters <- data.frame(sample_id = rownames(count_data), snf_cluster = group)

NMI_scores <- rankFeaturesByNMI(list(count_data, methyl_data, snv_data, cnv_data, splice_data), W)



# 3) output clusters with per sample  
write_tsv(snf_clusters, file = file.path(output_dir, "snf_clusters.tsv"))

#4) Concordance Plot
# combine samples map with IntNMF derived clusters
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(output_dir, "snf_clusters.tsv"))
mm_clusters <- mm_clusters %>%
  inner_join(samples_map)

# combine IntNMF clusters with RNA-derived molecular subtypes
anno_file_rna <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# 3) generate balloon plot with at least 5 samples in a group
# SNF clusters vs RNA-derived molecular subtypes
dat <- anno_file_rna %>%
  filter(!is.na(molecular_subtype)) %>%
  group_by(molecular_subtype, snf_cluster)  %>%
  dplyr::summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = snf_cluster, value = n, fill = 0) %>%
  column_to_rownames("molecular_subtype")
pdf(file = file.path(plots_dir, "snf_5_clusters_vs_molsubtype_balloonplot.pdf"), width = 10)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "SNF clusters vs Molecular subtypes",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()

# format survival data
surv_data <- anno_file_rna
surv_data$OS_status <- ifelse(surv_data$OS_status == "DECEASED", 1, 0)

# 4) survival curves stratified by RNA-derived molecular subtype 
surv_data$intnmf_cluster <- factor(surv_data$snf_cluster, levels = sort(as.numeric(unique(surv_data$snf_cluster))))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ snf_cluster, data = surv_data)
pdf(file = file.path(plots_dir, "survival_SNF_clusters.pdf"), height = 8, width = 10, onefile = FALSE)
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

# 5) survival curves stratified by RNA-derived molecular subtype 
surv_data$molecular_subtype <- factor(surv_data$molecular_subtype, levels = sort(unique(surv_data$molecular_subtype)))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ molecular_subtype, data = surv_data)
pdf(file = file.path(plots_dir, "survival_molsubtype.pdf"), height = 8, width = 10, onefile = FALSE)
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
