# prepare files for multi-modal clustering
#Script adapted from @KomalRathi MB MM-clustering using IntNMF
suppressPackageStartupMessages({
  library(tidyverse)
  library(datawizard)
  library(reshape2)
})

# define directories
ot_dir <- file.path("~/Documents", "GitHub", "OpenPedCan-analysis", "data", "v12")
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC", "data")
analysis_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC")

# output directory
output_dir <- file.path("input")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(analysis_dir, "utils", "filter_cnv.R"))

# cancer genes
cancer_genes <- readRDS(file.path(ot_dir, "cancer_gene_list.rds"))

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

# 2) CNV
cnv_dat <- file.path(data_dir, "All.gainloss_medullo.tsv") %>%
  data.table::fread()
cnv_dat <- cnv_dat %>% 
  dplyr::rename("Kids_First_Biospecimen_ID" = "BS_ID") %>%
  filter(Kids_First_Biospecimen_ID %in% histology_file$Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Biospecimen_ID, gene, log2) %>%
  inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id, tumor_ploidy, tumor_fraction), by = "Kids_First_Biospecimen_ID") %>%
  unique()

# get unique biospecimens per sample id
unique_ids <- cnv_dat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
cnv_dat <- cnv_dat %>% 
  filter(sample_id %in% unique_ids$sample_id,
         Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID)
cnv_samples <- cnv_dat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
# function to adjust copy number and status
adjust_cn <- function(x){
  
  # get tumor fraction and ploidy
  tumor_fraction <- unique(na.omit(x$tumor_fraction))
  tumor_ploidy <- unique(na.omit(x$tumor_ploidy))
  
  if(length(tumor_fraction) == 1 & length(tumor_ploidy) == 1){
    # calculate adjusted copy number if tumor fraction and ploidy info is available
    x$adjusted_cn <- (((2 ^ (x$log2) - (1 - tumor_fraction)) * tumor_ploidy) / tumor_fraction) - 0.5
    x$adjusted_cn <- round(x$adjusted_cn)
    x$adjusted_status <- ifelse(x$adjusted_cn == 0, "Complete Loss",
                                ifelse(x$adjusted_cn == 1, "Loss",
                                       ifelse(x$adjusted_cn %in% c(tumor_ploidy + 1:9), "Gain",
                                              ifelse(x$adjusted_cn >= 10, "Amplification", "Neutral"))))
    
    # replace old columns with new ones
    x <- x %>% 
      dplyr::rename("status" = "adjusted_status", # rename new columns
                    "copy_number" = "adjusted_cn")
    
  } 
}

# apply function to all samples in the consensus file
cnv_dat <- plyr::ddply(.data = cnv_dat, 
                       .variables = "sample_id", 
                       .fun = function(x) adjust_cn(x = x))
cnv_dat <- cnv_dat %>% 
  filter(status != "Neutral") %>% 
  dplyr::select(-c(log2)) %>%
  unique()

# subset to cancer genes
cnv_dat <- cnv_dat %>%
  dplyr::rename("hgnc_symbol" = "gene") %>%
  filter_cnv(myCancerGenes = cancer_genes)
cnv_dat <- acast(cnv_dat, sample_id ~ hgnc_symbol, value.var = "copy_number", fill = 0, fun.aggregate = max)

# 3) SNV
snv_dat <- file.path(data_dir, "snv-consensus-plus-hotspots_medullo.maf.tsv") %>%
  data.table::fread()
maf_nonsynonymous <- c("Missense_Mutation", "Frame_Shift_Del", "In_Frame_Ins",
                       "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", 
                       "In_Frame_Del", "Nonstop_Mutation", "Translation_Start_Site")
snv_dat <- snv_dat %>%
  filter(Variant_Classification %in% maf_nonsynonymous,
         Tumor_Sample_Barcode %in% histology_file$Kids_First_Biospecimen_ID)

# combine with histology to map sample id
snv_dat <- snv_dat %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") %>%
  inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id), by = 'Kids_First_Biospecimen_ID')

# get unique biospecimens per sample id
unique_ids <- snv_dat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
snv_dat <- snv_dat %>% 
  filter(sample_id %in% unique_ids$sample_id,
         Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID)
snv_samples <- snv_dat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix # Add fun.aggregate fun here -VK
snv_dat <- acast(snv_dat, sample_id ~ Hugo_Symbol, value.var = "Variant_Classification")
snv_dat[snv_dat > 0] <- 1 # anything > 0 should be coded as 1 


# 4) Methylation
# read beta-values
methyl_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds"))
methyl_data <- methyl_data %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
methyl_data <- methyl_data[complete.cases(methyl_data),]
# convert to long format
methyl_data <- melt(as.matrix(methyl_data), varnames = c("Probe_ID", "Kids_First_Biospecimen_ID"), value.name = "value")

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

# convert into sample id-probe matrix
methyl_data <- methyl_data %>%
  acast(sample_id ~ Probe_ID, value.var = "value")

# 5) Splice dataset
# read splice data
splice_file <- file.path(data_dir, "splice_data_medullo.rds")
splice_mat <- readRDS(splice_file)
splice_mat <- splice_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
splice_mat <- melt(as.matrix(splice_mat), varnames = c("Splice_Variant", "Kids_First_Biospecimen_ID")) 

# combine with histology to map sample id
splice_mat <- splice_mat %>%
  inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id), by = 'Kids_First_Biospecimen_ID')

# get unique biospecimens per sample id
unique_ids <- splice_mat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
splice_mat <- splice_mat %>% 
  filter(sample_id %in% unique_ids$sample_id,
         Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID)
splice_samples <- splice_mat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
splice_mat <- splice_mat %>%
  acast(Splice_Variant ~ sample_id, value.var = "value")

# subset to samples of interest
samples_of_interest <- intersect(colnames(count_mat), rownames(cnv_dat))
samples_of_interest <- intersect(samples_of_interest, rownames(snv_dat))
samples_of_interest <- intersect(samples_of_interest, rownames(methyl_data))
samples_of_interest <- intersect(samples_of_interest, colnames(splice_mat)) # 152 samples vs 91 from before

# now final filter/transformation on samples of interest

# 1) RNA
count_mat <- t(count_mat) %>% as.data.frame()
count_mat <- count_mat[samples_of_interest,]
# remove genes with 0 counts across all samples
#count_mat1 <- count_mat[,colSums(count_mat) > 0]  
write_tsv(as.data.frame(count_mat) %>% rownames_to_column(), file = file.path(output_dir, "norm_counts.tsv"))

# 2) CNV
cnv_dat <- cnv_dat[samples_of_interest,]
print(dim(cnv_dat)) # 1193
write_tsv(as.data.frame(cnv_dat) %>% rownames_to_column(), file = file.path(output_dir, "cnv_data.tsv"))

# 3) SNV
# filter out genes with low mutation rate 
snv_dat <- snv_dat[samples_of_interest,]
print(dim(snv_dat)) # 170
write_tsv(as.data.frame(snv_dat) %>% rownames_to_column(), file = file.path(output_dir, "snv_data.tsv"))

# 4) Methylation
methyl_data <- methyl_data[samples_of_interest,]
#methyl_data1 <- methyl_data[,colSums(is.na(methyl_data)) < nrow(methyl_data)]
print(dim(methyl_data)) 
write_tsv(as.data.frame(methyl_data) %>% rownames_to_column(), file = file.path(output_dir, "methyl_data.tsv"))

# 5) Splicing
splice_mat <- t(splice_mat) %>% as.data.frame()
splice_mat <- splice_mat[samples_of_interest,]
#splice_mat1 <- splice_mat[,colSums(splice_mat) > 0] 
print(dim(splice_mat)) # 1000
write_tsv(as.data.frame(splice_mat) %>% rownames_to_column(), file = file.path(output_dir, "splice_data.tsv"))

# final sample map
rna_samples <- count_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID")
methyl_samples <- methyl_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") 
cnv_samples <- cnv_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_CNV" = "Kids_First_Biospecimen_ID") 
snv_samples <- snv_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_SNV" = "Kids_First_Biospecimen_ID") 
splice_samples <- splice_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>% 
  dplyr::rename("Kids_First_Biospecimen_ID_Splice" = "Kids_First_Biospecimen_ID") 
sample_map <- rna_samples %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_RNA) %>% 
  inner_join(cnv_samples) %>% 
  inner_join(snv_samples) %>% 
  inner_join(methyl_samples) %>% 
  inner_join(splice_samples) 
write_tsv(sample_map, file = file.path(output_dir, "samples_map.tsv"))
