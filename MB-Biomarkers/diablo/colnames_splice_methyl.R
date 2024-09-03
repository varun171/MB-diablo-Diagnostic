# change the column names of splice data and methyl data

library(tidyverse)
library(dplyr)

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering/data")
input_dir <-
  file.path(root_dir, "multimodal_clustering/intNMF/input")

# read data

methyl_data <-
  read_tsv(file.path(input_dir, "methyl_data.tsv")) %>% column_to_rownames()
splice_data <-
  read_tsv(file.path(input_dir, "splice_data.tsv")) %>% column_to_rownames()

# Splice data colnames

original_colnames <- colnames(splice_data)

# Generate new column names
new_colnames <-
  gsub("_(.*)$", "", original_colnames)  # Remove everything after the first underscore
new_colnames <-
  paste0(new_colnames, "_AT_", seq_along(new_colnames))  # Add _AT_ and incrementing numbers

# Update the column names of the dataframe
colnames(splice_data) <- new_colnames

# Create a data frame to map original column names to new column names
colname_mapping <-
  data.frame(Original = original_colnames, New = new_colnames)

# Save the mapping to a CSV file (optional)
write.csv(
  colname_mapping,
  file = file.path(input_dir, "column_name_mapping.csv"),
  row.names = FALSE
)

write_tsv(
  as.data.frame(splice_data) %>% rownames_to_column(),
  file = file.path(input_dir, "splice_data1.tsv")
)

methyl_annot <-
  data.table::fread(file.path(data_dir, "infinium.gencode.v39.probe.annotations.tsv"))
gtf_file <-
  file.path(data_dir, "gencode.v39.primary_assembly.annotation.pc.tsv")
gencode_gtf <- read_tsv(gtf_file)

methyl_annot <- methyl_annot %>%
  filter(Gene_symbol %in% gencode_gtf$gene_name) %>% # filter to protein coding genes to reduce features
  dplyr::select(Gene_symbol, Probe_ID) %>%
  unique()

# Get the matching Probe_IDs
matching_probe_ids <-
  methyl_annot$Probe_ID[methyl_annot$Probe_ID %in% colnames(methyl_data)]

# Filter methyl_data to keep only matching columns
methyl_data <- methyl_data[, matching_probe_ids, drop = FALSE]

# Get the matching Gene_symbols for the selected columns
matching_symbols <-
  methyl_annot$Gene_symbol[methyl_annot$Probe_ID %in% matching_probe_ids]

# Update the column names with matching symbols
colnames(methyl_data) <-
  paste(colnames(methyl_data), matching_symbols, sep = "_")

# save the file
write_tsv(
  as.data.frame(methyl_data) %>% rownames_to_column(),
  file = file.path(input_dir, "methyl_data1.tsv")
)
