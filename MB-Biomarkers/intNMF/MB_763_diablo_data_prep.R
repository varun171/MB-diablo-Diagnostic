suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(datawizard)
  library(reshape2)
  library(survival)
  library(survminer)
  library(ggsankey)
  library(mixOmics)
  library(data.table)
})

data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo")
input_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/input")
nmf_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/results")
output_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/results/MB763")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/plots/MB763")



# Read the train data (Cavali et al. 2017)
count_data <- fread("~/Documents/cancer-cell-MB/GSE85217_M_exp_763_MB.txt") 
methyl_data <- fread("~/Documents/cancer-cell-MB/GSE85212_Methylation_763_beta_values.txt.gz") 
methyl_data <- methyl_data %>% dplyr::rename(Probe_ID = V1)

MB_Expr_Anno <- fread("~/Documents/cancer-cell-MB/GSE85217_series_matrix.txt", skip =30) 
MB_Anno <- MB_Expr_Anno %>% filter(row_number() %in% c(10, 11)) %>% 
                     t() %>% 
                     as.data.frame() %>% 
                     rownames_to_column(var = "sample_id") %>%
                     `colnames<-`(c("sample_id", "subgroup", "subtype"))
samples_map <- MB_Anno[-1, ] 
rownames(samples_map) <- NULL
# Save MB_anno as a TSV file
write.table(samples_map, file = file.path(data_dir, "samples_763map.tsv"), sep = "\t", row.names = FALSE)


#MB_Methyl_Anno <- fread("~/Documents/cancer-cell-MB/GSE85212_series_matrix.txt", skip =30) 

#gencode_v39_gtf <- file.path("gencode.v39.primary_assembly.annotation.gtf") %>%
#  rtracklayer::import() %>%
#  as.data.frame() %>%
#  dplyr::select(gene_id, gene_name, gene_type) %>%
#  filter(!grepl("ENSG", gene_name), gene_type == "protein_coding") %>%
#  unique()
# write_tsv(gencode_v39_gtf, file = "gencode.v39.primary_assembly.annotation.pc.tsv")


# read gtf and filter to protein coding 
gencode_gtf <- fread(file = file.path(data_dir, "gencode.v39.primary_assembly.annotation.pc.tsv"))
 
# filter expression count file to contain only protein coding genes
count_data <- count_data %>%
  dplyr::filter(HGNC_symbol_from_ensemblv77 %in% c(gencode_gtf$gene_name)) %>%
  dplyr::select(-1, -2, -3, -5) %>%
  `colnames<-`(c("gene_name", colnames(.)[-1]))%>%
  distinct(gene_name, .keep_all = TRUE)

# filter methyl data to contain only protein coding genes
methyl_annot <- data.table::fread(file.path(data_dir, "infinium.gencode.v39.probe.annotations.tsv"))
methyl_annot <- methyl_annot %>%
   filter(Gene_symbol %in% gencode_gtf$gene_name) %>% # filter to protein coding genes to reduce features
   dplyr::select(Gene_symbol, Probe_ID) %>%
   unique()
methyl_data <- methyl_data %>%
filter(Probe_ID %in% methyl_annot$Probe_ID) %>%
  distinct(Probe_ID, .keep_all = TRUE)

# feature selections and the save the data sets 
count_data <- count_data %>% column_to_rownames(var = "gene_name") %>% t() %>% as.data.frame()
write.table(count_data, file = file.path(data_dir, "count_763data.tsv"), sep = "\t", row.names = FALSE)

# top 1000 most variable genes
num_genes = 1000
keep <- apply(count_data, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
count_data <- count_data[,colnames(count_data) %in% keep]

# rank transformation
count_data <- ranktransform(t(count_data) %>% as.data.frame())
count_data <- t(count_data) %>% as.data.frame()
print(dim(count_data)) # 1000
write_tsv(as.data.frame(count_data) %>% rownames_to_column(), file = file.path(data_dir, "norm_763counts.tsv"))


methyl_data <- methyl_data %>% column_to_rownames(var = "Probe_ID") %>% t() %>% as.data.frame()
saveRDS(methyl_data, file = file.path(data_dir, "methyl_763data.rds"))

num_genes = 1000
keep <- apply(methyl_data, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
methyl_data <- methyl_data[,colnames(methyl_data) %in% keep]
print(dim(methyl_data)) # 1000
write_tsv(as.data.frame(methyl_data) %>% rownames_to_column(), file = file.path(data_dir, "methyl_763data.tsv"))


