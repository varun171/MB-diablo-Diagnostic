# prepare files for multi-modal clustering
suppressPackageStartupMessages({
  library(tidyverse)
  library(datawizard)
  library(reshape2)
})

# define directories
# define directories
ot_dir <- file.path("~/Documents/v15-Pedcan")
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF")


# output directory
output_dir <- file.path(analysis_dir, "input")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
#source(file.path(root_dir, "utils", "filter_cnv.R"))

# cancer genes
#cancer_genes <- readRDS(file.path(data_dir, "cancer_gene_list.rds"))

# read histology file for tumor fraction and tumor ploidy
histology_file <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() 

# 1) read count data
count_file <- file.path(data_dir, "gene-counts-rsem-expected_count_medullo-collapsed.rds")
count_mat <- readRDS(count_file)
count_mat <- count_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
count_mat <- melt(as.matrix(count_mat), varnames = c("Gene", "Kids_First_Biospecimen_ID")) 

# read gtf and filter to protein coding 
gtf_file <- file.path(data_dir, "gencode.v39.primary_assembly.annotation.pc.tsv")
gencode_gtf <- read_tsv(gtf_file)

# filter expression count file to contain only protein coding genes
count_mat <- count_mat %>%
  filter(Gene %in% c(gencode_gtf$gene_name))

# combine with histology to map sample id
count_mat <- count_mat %>%
  inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id), by = 'Kids_First_Biospecimen_ID')

# get unique biospecimens per sample id
unique_ids <- count_mat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
count_mat <- count_mat %>% 
  filter(sample_id %in% unique_ids$sample_id,
         Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID)
count_samples <- count_mat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
count_mat <- count_mat %>%
  acast(Gene ~ sample_id, value.var = "value")



# 4) Methylation
# read beta-values
methyl_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds"))
methyl_data <- methyl_data %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
methyl_data <- methyl_data[complete.cases(methyl_data),]

# convert to long format
methyl_data <- melt(as.matrix(methyl_data), varnames = c("Probe_ID", "Kids_First_Biospecimen_ID"), value.name = "value")

# # reduce to protein coding genes and gene symbol by taking mean/median value of probe id per gene symbol
# methyl_annot <- data.table::fread(file.path(data_dir, "infinium.gencode.v39.probe.annotations.tsv.gz"))
# methyl_annot <- methyl_annot %>%
#   filter(Gene_symbol %in% gencode_gtf$gene_name) %>% # filter to protein coding genes to reduce features
#   dplyr::select(Gene_symbol, Probe_ID) %>%
#   unique()
# methyl_data <- methyl_data %>%
#   filter(Probe_ID %in% methyl_annot$Probe_ID) %>%
#   inner_join(methyl_annot, by = "Probe_ID")
methyl_data <- methyl_data %>%
  inner_join(histology_file %>% 
               filter(experimental_strategy == "Methylation") %>% 
               dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>% 
               unique(), by = 'Kids_First_Biospecimen_ID') 

# get unique biospecimens per sample id
unique_ids <- histology_file %>%
  filter(Kids_First_Biospecimen_ID %in% methyl_data$Kids_First_Biospecimen_ID,
         sample_id %in% methyl_data$sample_id) %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
methyl_data <- methyl_data %>% 
  filter(sample_id %in% unique_ids$sample_id,
         Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID)
methyl_samples <- methyl_data %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# # summarise mean values
# methyl_data <- methyl_data %>%
#   group_by(sample_id, Gene_symbol) %>%
#   dplyr::summarise(mean_val = mean(value),
#                    median_val = median(value))

# convert into sample id-probe matrix
methyl_data <- methyl_data %>%
  acast(sample_id ~ Probe_ID, value.var = "value")


# subset to samples of interest
samples_of_interest <- intersect(colnames(count_mat), rownames(methyl_data))

# now final filter/transformation on samples of interest

# 1) RNA
count_mat <- t(count_mat) %>% as.data.frame()
count_mat <- count_mat[samples_of_interest,]
# remove genes with 0 counts across all samples
count_mat <- count_mat[,colSums(count_mat) > 0]  

# top 1000 most variable genes
num_genes = 1000
keep <- apply(count_mat, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
count_mat <- count_mat[,colnames(count_mat) %in% keep]

# rank transformation
count_mat <- ranktransform(t(count_mat) %>% as.data.frame())
count_mat <- t(count_mat) %>% as.data.frame()
print(dim(count_mat)) # 1000
write_tsv(as.data.frame(count_mat) %>% rownames_to_column(), file = file.path(data_dir, "norm_counts.tsv"))

# 4) Methylation
methyl_data <- methyl_data[samples_of_interest,]
methyl_data <- methyl_data[,colSums(is.na(methyl_data)) < nrow(methyl_data)]

# top 1000 most variable features
num_genes = 1000
keep <- apply(methyl_data, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
methyl_data <- methyl_data[,colnames(methyl_data) %in% keep]
print(dim(methyl_data)) # 1000
write_tsv(as.data.frame(methyl_data) %>% rownames_to_column(), file = file.path(data_dir, "methyl_data.tsv"))

# final sample map
rna_samples <- count_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID")
methyl_samples <- methyl_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") 
 
sample_map <- rna_samples %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_RNA) %>% 
  inner_join(methyl_samples) 
write_tsv(sample_map, file = file.path(data_dir, "samples_map.tsv"))
