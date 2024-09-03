suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(datawizard)
  library(reshape2)
  library(survival)
  library(survminer)
  library(ggsankey)
  library(DESeq2)
  library(AnnotationHub)   
  library(signatureSearch)
  library(ExperimentHub)
  library(rhdf5)
  library(PubChemR)
  library(data.table)
  library(NOISeq)
})


data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo")
input_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/input")
nmf_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/results")
output_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/plots")

# read data

cluster1_genes <- read_tsv(file.path(output_dir, "significant_genes_cluster_1.tsv"))  
cluster2_genes <- read_tsv(file.path(output_dir, "significant_genes_cluster_2.tsv")) 
cluster3_genes <- read_tsv(file.path(output_dir, "significant_genes_cluster_3.tsv"))  
cluster4_genes <- read_tsv(file.path(output_dir, "significant_genes_cluster_4.tsv")) 
cluster5_genes <- read_tsv(file.path(output_dir, "significant_genes_cluster_5.tsv"))  

final_cluster <- rbind(cluster1_genes,cluster2_genes, cluster3_genes, cluster4_genes, cluster5_genes)
final_cluster <- final_cluster[!duplicated(final_cluster$Gene), ]
#final_cluster <- cluster1_genes
###########################################################################################################

# histology

samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(nmf_dir, "intnmf_clusters.tsv")) 
mm_clusters_diablo <- mm_clusters %>%
  inner_join(samples_map)

mm_cluster_diablo <- mm_clusters_diablo %>% arrange(Kids_First_Biospecimen_ID_RNA) 

tumor_hist_df <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, gtex_group, gtex_subgroup) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  arrange(Kids_First_Biospecimen_ID_RNA) %>%  # Sort the rows by Kids_First_Biospecimen_ID_RNA
  semi_join(mm_cluster_diablo, by = "Kids_First_Biospecimen_ID_RNA") %>%
  mutate(gtex_group = "MB", gtex_subgroup = mm_cluster_diablo$intnmf_cluster)


hist_df <- file.path(input_dir, "20231115_release.annotated_histologies.tsv") %>%
  read_tsv()

# Filter GTex
gtex_hist_df <- hist_df %>%
  filter(
    grepl("GTEx", cohort),
    grepl("Brain", gtex_group),
    experimental_strategy == "RNA-Seq") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, gtex_group, gtex_subgroup) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") 

final_hist_df <- rbind(tumor_hist_df, gtex_hist_df)

tpm_data <- readRDS(file.path(input_dir, "gene-expression-rsem-tpm-collapsed-v12.rds"))

# genes of interest
#genes_of_interest <- c("FZD1", "OTX2", "GABRB3", "EZH2", "BAG6", "ZBTB18", "RPS7", "RPS17","SFRP1", "UNC5D", "SFRP1","PDLIM3","HNRNPH3", "WIF1")
#genes_of_interest <- c("FZD1", "OTX2", "GABRB3", "GLCE", "WIF1", "CDK5R1", "UNC5D", "SFRP1","PDLIM3","HNRNPH3")
#genes_of_interest <- c("HDAC2", "EPHB2", "PRKDC", "CACNA2D1", "CDK6", "FZD1", "TNNI1", "SIPA1L2","AKT3", "KIF5C", "IGFBP5","PEG3")

genes_of_interest <- final_cluster$Gene
# Assuming 'genes_of_interest' vector is defined and 'tpm' is your data frame
tpm_data <- tpm_data[rownames(tpm_data) %in% genes_of_interest, ]

# filtered exper_tpm data for filter to genes of interest
tpm_data <- log2(tpm_data + 1)
tpm_data <- tpm_data %>%
  rownames_to_column("gene_symbol") %>%
  filter(gene_symbol %in% genes_of_interest) %>% column_to_rownames("gene_symbol")%>% t() %>% as.data.frame() %>% rownames_to_column("Kids_First_Biospecimen_ID_RNA")

# join with histology
expr_tpm <- tpm_data %>%
  inner_join(final_hist_df, by= "Kids_First_Biospecimen_ID_RNA")

#gene_data <- expr_tpm %>%
#  gather(key = "gene", value = "expression", genes_of_interest)

# Check if any genes_of_interest exist in expr_tpm
#missing_genes <- genes_of_interest[!(genes_of_interest %in% colnames(expr_tpm))]

# Proceed with gather operation for existing genes
gene_data <- expr_tpm %>%
  gather(key = "gene", value = "expression", any_of(genes_of_interest))

# Create separate t-test results for each gene
t_test_results <- gene_data %>%
  group_by(gene) %>%
  summarize(p_value = t.test(expression ~ gtex_group)$p.value)

# Save the t-test results as a text file
write.table(t_test_results, file = file.path(output_dir, "t_test_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Create plots for each gene and save in a single PDF file
pdf(file.path(plots_dir, "gene_expression_violin_plots_KM_IntNMF_1.pdf"), height = 10, width = 15)
for (gene_name in unique(gene_data$gene)) {
  p <- ggplot(subset(gene_data, gene == gene_name), aes(x = gtex_subgroup, y = expression, fill = gtex_subgroup)) +
    geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.1, fill = "white", position = position_dodge(width = 0.9)) +
    labs(x = "Subgroup", y = "Gene Expression") +
    facet_wrap(~gene, scales = "free") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
  print(p)
}
dev.off()
############################################################################################
# sankey plot with IntNMF 14 clusters

samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))

intnmf_clusters <- read_tsv(file.path(nmf_dir, "intnmf_clusters.tsv"))
movics_clusters <- read_tsv(file.path(analysis_dir, "consensus_clusters.tsv")) %>% dplyr::rename("sample_id" = "samID", "movics_cluster" = "clust")
intnmf_movics_clusters <- intnmf_clusters %>% inner_join(movics_clusters)
intnmf_movics_clusters <- intnmf_movics_clusters %>% inner_join(samples_map)
intnmf_movics_clusters <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(intnmf_movics_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# Assuming intnmf_movics_clusters is your data frame
# Summarize data
links <- intnmf_movics_clusters %>%
  make_long(intnmf_cluster, molecular_subtype, movics_cluster) 


# Convert next_node to character type
links$next_node <- as.character(links$next_node)

# Create Sankey plot
pl <- ggplot(links, aes(x = x,                        
                     next_x = next_x,                                     
                     node = node,
                     next_node = next_node,        
                     fill = factor(node),
                     label = node))             # This Creates a label for each node

pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(Size = 3, 
                             color = "black", 
                             fill = "white") # This specifies the Label format for ea

# Save the plot as a PDF
ggsave("sankey_plot.pdf", plot = pl, width = 8, height = 6)



#############################################################################################

################################################################################################
# methyl data expression # memory exhausted

# Read the methyl data file
methyl_data <- readRDS("~/Documents/GitHub/OpenPedCan-analysis/data/v12/methyl-m-values.rds")

tumor_hist_df <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, gtex_group, gtex_subgroup) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  arrange(Kids_First_Biospecimen_ID_Methyl) %>%  # Sort the rows by Kids_First_Biospecimen_ID_RNA
  semi_join(mm_cluster_diablo, by = "Kids_First_Biospecimen_ID_Methyl") %>%
  mutate(gtex_group = "MB", gtex_subgroup = mm_cluster_diablo$intnmf_cluster)

methyl_data <- methyl_data %>% as.data.frame() %>% column_to_rownames(var = "Probe_ID") %>%  
  dplyr::select(any_of(tumor_hist_df$Kids_First_Biospecimen_ID_methyl))
#methyl_data <- methyl_data[complete.cases(methyl_data),]
write_rds(methyl_data, file = file.path(input_dir, "mb_methyl_m-value.rds"))

# rename column names

################################################################################
# differential expression between CHOP-MB count data with Intnmf clusters and Gtex normal brain 
# histology

data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo")
input_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/input")
nmf_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/results")
output_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/results/CHOP-DEG")
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/plots/CHOP-DEG")
dir.create(plots_dir, showWarnings = F, recursive = T)


samples_map <- read_tsv(file.path(data_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(nmf_dir, "intnmf_clusters.tsv")) 
mm_clusters$intnmf_cluster <-
  paste("cluster", mm_clusters$intnmf_cluster, sep = "")
mm_clusters_diablo <- mm_clusters %>%
  inner_join(samples_map)

mm_cluster_diablo <- mm_clusters_diablo %>% arrange(Kids_First_Biospecimen_ID_RNA) 

tumor_hist_df <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, gtex_group, gtex_subgroup) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  arrange(Kids_First_Biospecimen_ID_RNA) %>%  # Sort the rows by Kids_First_Biospecimen_ID_RNA
  semi_join(mm_cluster_diablo, by = "Kids_First_Biospecimen_ID_RNA") %>%
  mutate(gtex_group = "MB", gtex_subgroup = mm_cluster_diablo$intnmf_cluster)


hist_df <- file.path(data_dir, "histologies.tsv") %>%
  read_tsv()

# Filter GTex
gtex_hist_df <- hist_df %>%
  filter(
    grepl("GTEx", cohort),
    grepl("Brain", gtex_group),
    experimental_strategy == "RNA-Seq") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, gtex_group, gtex_subgroup) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  mutate(gtex_subgroup = "Normal_Brain")

final_hist_df <- rbind(tumor_hist_df, gtex_hist_df)

count_data <- readRDS(file.path(data_dir, "gene-counts-rsem-expected_count_medullo-collapsed.rds"))
gtex_data <-  readRDS(file.path(data_dir, "gtex_gene-counts-rsem-expected_count-collapsed.rds"))
final_count_data <- count_data[, intersect(colnames(count_data), final_hist_df$Kids_First_Biospecimen_ID_RNA)]
final_gtex_data <- gtex_data[, intersect(colnames(gtex_data), final_hist_df$Kids_First_Biospecimen_ID_RNA)]

# Find common row names
common_row_names <- intersect(rownames(final_count_data), rownames(final_gtex_data))

# Subset the data frames to include only the common row names
final_count_df_common <- final_count_data[common_row_names, ]
final_gtex_df_common <- final_gtex_data[common_row_names, ]

# Combine the two data frames using cbind
combined_df <- cbind(final_count_df_common, final_gtex_df_common)

# Round the count data to the nearest integer
combined_df <- round(combined_df)

# filter expression count file to contain only protein coding gene
# read gtf and filter to protein coding 
gencode_gtf <- read_tsv(file = file.path(data_dir, "gencode.v39.primary_assembly.annotation.pc.tsv")) %>%
  as.data.frame() %>% 
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

combined_df <- combined_df %>%
  filter(rownames(combined_df) %in% gencode_gtf$gene_name)

#################################################################################
#Deseq2
# PCA to check variability
myPCA <- pca(combined_df, method = "svd", center = TRUE, scale = "uv")
loadings_matrix <- loadings(myPCA)
loadings_df <- as.data.frame(loadings_matrix)

# Add intnmf_cluster as a column
loadings_df$subgroup <- final_hist_df$gtex_subgroup

# Plot using ggplot2
p <- ggplot(loadings_df, aes(x = PC1, y = PC2, color = subgroup)) +
  geom_point() +
  labs(title = "PCA Plot with intnmf_cluster", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_discrete(name = "subgroup")
ggsave(filename = file.path(plots_dir, "PCA_combat_plot.pdf"), plot = p, width = 10, height = 10)

#combined_df <- read.table(file.path(output_dir, "Normalized_count_data.tsv"), header = TRUE, sep = "\t") %>%
# column_to_rownames(var = "X")

# use DCGA to filter out low count, low variance features
combined_df <- DGCA::filterGenes(inputMat = combined_df, 
                              filterTypes = c("central", "dispersion"),
                              filterDispersionType = "cv", 
                              filterDispersionPercentile = 0.2,
                              sequential = TRUE)


# Convert 'gtex_subgroup' to factor
final_hist_df$gtex_subgroup <- as.factor(final_hist_df$gtex_subgroup)

# Round the count data to the nearest integer
#combined_df <- round(combined_df)

# Create DESeq2 object for all subtypes
dds <- DESeqDataSetFromMatrix(countData = combined_df,
                              colData = final_hist_df,
                              design = ~ gtex_subgroup)

# Run DESeq2 normalization and analysis
dds <- DESeq(dds)

# Loop through each subtype and compare it against "Normal_Brain" subtype
for (subtype in unique(final_hist_df$gtex_subgroup)) {
  if (subtype != "Normal_Brain") {
    # Specify the contrast
    contrast <- c("gtex_subgroup", subtype, "Normal_Brain")
    
    # Get differential expression results for the current subtype vs. "Normal_Brain"
    res <- results(dds, contrast = contrast)
    
    # Perform further analysis or save results as needed
    # For example, you can save results to a file
    write.csv(as.data.frame(res), file = file.path(output_dir, paste0(subtype, "_vs_Normal_Brain_deseq2_results_normalized.csv")))
  }
}
##################################################################
# Noiseq analyssi to double check DE

# Format data for NOIseq
mydata <- readData(data = combined_df, factors = final_hist_df %>% column_to_rownames("Kids_First_Biospecimen_ID_RNA"))

# PCA before batch correction
myPCA <- dat(mydata, type = "PCA")
pdf(file = file.path(plots_dir, "noiseq_pca.pdf"), width = 10, height = 10, onefile = TRUE)
explo.plot(myPCA, factor = "gtex_subgroup")
dev.off()

# Batch correction
mydata_corr <- ARSyNseq(mydata, factor = "gtex_subgroup", batch = FALSE, norm = "uqua", logtransf = FALSE)

# PCA after batch correction
myPCA_corr <- dat(mydata_corr, type = "PCA")
pdf(file = file.path(plots_dir, "noiseq_pca_after_correction.pdf"), width = 10, height = 10, onefile = TRUE)
explo.plot(myPCA_corr, factor = "gtex_subgroup")
dev.off()

# Format data for NOIseq
mydata_corr <- readData(data = exprs(mydata_corr), factors = final_hist_df %>% column_to_rownames("Kids_First_Biospecimen_ID_RNA"))
#write.table(data, file.path(output_dir, "Normalized_count_data.tsv"), sep = "\t", quote = FALSE, col.names = NA)
# Read the TSV file back into R
#combined_df <- read.table(file.path(output_dir, "Normalized_count_data.tsv"), header = TRUE, sep = "\t") %>%
#  column_to_rownames(var = "X")

# run though each cluster vs Normal Brain
for (subtype in unique(final_hist_df$gtex_subgroup)) {
  if (subtype != "Normal_Brain") {
  
    # Subset the data for the current subtype
    subset_data <- mydata_corr[, final_hist_df$gtex_subgroup %in% c(subtype, "Normal_Brain")]
    
    # Run noiseqbio
    mynoiseqbio <- noiseqbio(input = subset_data,
                             k = 0.5,
                             norm = "uqua",
                             factor = "gtex_subgroup",
                             conditions = c(subtype, "Normal_Brain"),
                             r = 20,
                             adj = 1.5,
                             plot = FALSE,
                             a0per = 0.9,
                             random.seed = 12345,
                             filter = 2)
    
    # Get DE genes
    noiseq_output <- degenes(mynoiseqbio, q = 0.8, M = NULL)
    
     # Add direction column
    noiseq_output <- noiseq_output %>%
      mutate(direction = ifelse(log2FC > 0, "up", "down"),
             comparison = paste(subtype, "vs_Normal_Brain")) %>%
      rownames_to_column("gene") %>%
      dplyr::select(comparison, gene, prob, log2FC, direction)
    
    # Save results to a file
    write.csv(as.data.frame(noiseq_output), file = file.path(output_dir, paste0(subtype, "_vs_Normal_Brain_noiseq_results.csv")))
  }
}

###################################################################
# Lincs analysis

# source function
source(file.path(analysis_dir, "utils", "drug_barplots.R"))
source(file.path(analysis_dir, "utils", "lincs_analysis.R"))

input_df1 <- fread(file.path(output_dir, "cluster3_vs_Normal_Brain_deseq2_results.csv"))
#input_df2 <- fread(file.path(output_dir, "cluster3_vs_Normal_Brain_deseq2_results_normalized.csv"))
#input_df3 <- fread(file.path(output_dir, "cluster3_vs_Normal_Brain_noiseq_results.csv"))
#sig_gene_df <- input_df[input_df$log2FoldChange > 5 & input_df$padj < 0.005, ]

# get LINCS data
eh <- ExperimentHub()
lincs <- eh[["EH3226"]]

# LINCS analysis (reverse)
result <- lincs_analysis(input = input_df,
                   num_features = 200,
                   trend_val = "down",
                   wtcs_fdr_cutoff = 0.05,
                   output_dir = output_dir,
                   plots_dir = plots_dir,
                   db_path = lincs)

#######################################################################################
#reterive smiles from pubchem
qsig_df <- fread(file.path(output_dir, "qSig_output.txt"))
touchstone_df <- fread(file.path(data_dir, "touchstone_data.txt"))
qsig <- merge(qsig_df, touchstone_df, by.x = "pert", by.y = "Name", all = FALSE)
qsig <- unique(qsig, by = "pert")
qsig <- qsig[qsig$WTCS_FDR < 0.05 & qsig$NCS < 0, ]
unique_pert <- unique(qsig$pert)

cids <- get_cids(unique_pert)
group3_cids <- cids[cids$CID != "No CID", ]
final_cids <- group3_cids$CID
writeLines(final_cids, file.path(output_dir, "final_cids.txt"))
