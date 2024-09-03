# script to create subset matrices containing only MB data from OT v15

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# directory
ot_dir <- file.path("~/Documents/v15-Pedcan")
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF")
cancer_group_of_interest <- "Medulloblastoma"

# histology (adapt)
hist_df <- file.path(ot_dir, "histologies.tsv") %>% # test ot_dir instead of data_dir
  read_tsv() %>%
  filter(short_histology %in% cancer_group_of_interest)
write_tsv(hist_df, file = file.path(data_dir, "histologies_medullo.tsv"))

# RNA-seq counts (v15)
rna_seq_counts <- file.path(ot_dir, "gene-counts-rsem-expected_count-collapsed.rds") %>%
  readRDS() %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))
saveRDS(rna_seq_counts, file = file.path(data_dir, "gene-counts-rsem-expected_count_medullo-collapsed.rds"))


# Gencode v39 GTF subset to protein coding genes only
gencode_v39_gtf <- file.path(ot_dir, "gencode.v39.primary_assembly.annotation.gtf.gz") %>%
  rtracklayer::import() %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(!grepl("ENSG", gene_name),
         gene_type == "protein_coding") %>%
  unique()
write_tsv(gencode_v39_gtf, file = file.path(data_dir, "gencode.v39.primary_assembly.annotation.pc.tsv"))

# Methylation beta-values (v15)
methyl_beta_values <- file.path(ot_dir, "methyl-beta-values.rds") %>%
  readRDS() %>%
  distinct(Probe_ID, .keep_all = TRUE) %>%
  column_to_rownames("Probe_ID") %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))
saveRDS(methyl_beta_values, file = file.path(data_dir, "methyl-beta-values_medullo.rds"))

