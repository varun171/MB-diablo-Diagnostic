# script to perform downstream feature level analysis for cnv data

suppressPackageStartupMessages({
  library(tidyverse)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "multimodal_clustering")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results", "feature_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# read full output
nmf_output <- readRDS(file.path(analysis_dir, "results", "intnmf_best_fit.rds"))

# cnv data: pull top 10 features by NMF weight per cluster 
cnv_selected <- nmf_output$H$H4
cnv_selected <- cnv_selected %>%
  as.data.frame() %>%
  rownames_to_column("intnmf_cluster") %>%
  gather(key = "feature", value = "intnmf_weight", -c(intnmf_cluster)) %>% 
  mutate(intnmf_cluster = as.numeric(intnmf_cluster)) %>%
  group_by(intnmf_cluster) %>%
  dplyr::arrange(intnmf_cluster, desc(intnmf_weight)) %>%
  slice_head(n = 10)

# get sample level data for cnv (needed to map molecular subtype information to cluster)
cnv_data <- read_tsv(file.path(input_dir, "cnv_data.tsv")) %>% 
  column_to_rownames() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  gather(key = "feature", value = "copy_number_value", -c(sample_id))

# add cluster information to sample level data
mm_clusters <- read_tsv(file.path(analysis_dir, "results", "intnmf_clusters.tsv"))
cnv_data <- cnv_data %>%
  inner_join(mm_clusters) 

# add molecular subtype to sample level data
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
anno_file <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, dkfz_v11_methylation_subclass, dkfz_v12_methylation_subclass) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(samples_map, by = "Kids_First_Biospecimen_ID_Methyl") %>%
  dplyr::select(sample_id, dkfz_v12_methylation_subclass)
cnv_data <- cnv_data %>%
  inner_join(anno_file)

# combine with top 10 IntNMF features per cluster information
cnv_data <- cnv_data %>%
  inner_join(cnv_selected, by = c("feature", "intnmf_cluster"))

# now filter to group3/4 molecular subtypes only
cnv_data <- cnv_data %>%
  filter(grepl("G34", dkfz_v12_methylation_subclass))

# filter to features that have an amplification (>= 10 is for amplification)
cnv_data <- cnv_data %>%
  filter(copy_number_value >= 10) %>%
  dplyr::select(-c(copy_number_value))

# collapse sample identifiers to comma-separated values to reduce table
cnv_data <- cnv_data %>% 
  group_by(dkfz_v12_methylation_subclass, intnmf_cluster, intnmf_weight, feature) %>% 
  dplyr::summarise(sample_id = toString(sample_id))

# Using the summary from Figure 8 from PMID: 31076851, compare gene lists from new group 3/4 subgroups from intNMF to the gene lists from G3/4 subgroups reported in Figure 8.
# figure 8 amplifications
fig_8_amp <- data.frame(group = c("MB_G34_I", "MB_G34_II",  "MB_G34_III", "MB_G34_IV", 
                                  "MB_G34_V", "MB_G34_VI", "MB_G34_VII", "MB_G34_VIII"), 
                        feature = c("OTX2", "MYC", "MYC, MYCN", NA, 
                                    "MYCN", "MYCN", NA, NA))
fig_8_amp <- fig_8_amp %>%
  mutate(feature = strsplit(feature, ", ")) %>% 
  unnest(feature) # split comma-sep values into new rows 

# combine with cnv data
df <- cnv_data %>%
  left_join(fig_8_amp %>% mutate(fig8_overlap = TRUE), by = c("feature", "dkfz_v12_methylation_subclass" = "group")) %>%
  dplyr::select(intnmf_cluster, dkfz_v12_methylation_subclass, feature, intnmf_weight, sample_id, fig8_overlap) %>%
  dplyr::arrange(fig8_overlap, intnmf_cluster)

# save as tsv file
df %>%
  write_tsv(file = file.path(output_dir, "cnv_top10_feature_group34_comparison.tsv"))
