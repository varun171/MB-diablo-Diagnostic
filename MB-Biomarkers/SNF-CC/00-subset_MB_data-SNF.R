# script to create subset matrices containing only MB data from OT v12
# Script adapted from @KomalRathi MB MM-clustering using IntNMF

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})

# directory
ot_dir <- file.path("~/Documents", "GitHub", "OpenPedCan-analysis", "data", "v12")
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC", "data")
cancer_group_of_interest <- "Medulloblastoma"

# histology (adapt)
hist_df <- file.path(ot_dir, "histologies.tsv") %>% # test ot_dir instead of data_dir
  read_tsv() %>%
  filter(short_histology %in% cancer_group_of_interest)
write_tsv(hist_df, file = file.path(data_dir, "histologies_medullo.tsv"))

# RNA-seq counts (v12)
rna_seq_counts <- file.path(ot_dir, "gene-counts-rsem-expected_count-collapsed.rds") %>%
  readRDS() %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))
saveRDS(rna_seq_counts, file = file.path(data_dir, "gene-counts-rsem-expected_count_medullo-collapsed.rds"))

# SNV (v12)
snv_data <- file.path(ot_dir, "snv-consensus-plus-hotspots.maf.tsv.gz") %>%
  data.table::fread() %>%
  filter(Tumor_Sample_Barcode %in% hist_df$Kids_First_Biospecimen_ID)
write_tsv(snv_data, file = file.path(data_dir, "snv-consensus-plus-hotspots_medullo.maf.tsv"))

# Gencode v39 GTF subset to protein coding genes only
gencode_v39_gtf <- file.path(ot_dir, "gencode.v39.primary_assembly.annotation.gtf.gz") %>%
  rtracklayer::import() %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(!grepl("ENSG", gene_name),
         gene_type == "protein_coding") %>%
  unique()
write_tsv(gencode_v39_gtf, file = file.path(data_dir, "gencode.v39.primary_assembly.annotation.pc.tsv"))

# Methylation m-values (v12)
methyl_m_values <- file.path(ot_dir, "methyl-m-values.rds") %>%
  readRDS() %>%
  column_to_rownames("Probe_ID") %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))
saveRDS(methyl_m_values, file = file.path(data_dir, "methyl-m-values_medullo.rds"))

# Methylation beta-values (v12)
methyl_beta_values <- file.path(ot_dir, "methyl-beta-values.rds") %>%
  readRDS() %>%
  column_to_rownames("Probe_ID") %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))
saveRDS(methyl_beta_values, file = file.path(data_dir, "methyl-beta-values_medullo.rds"))

# Splice data (filtered to functional sites from Ammar)
splice_data <- file.path(ot_dir, "splice_events_pan_cancer_functional_filtered.rds") %>%
  readRDS() %>%
  column_to_rownames("Splice_ID") %>%
  dplyr::select(any_of(hist_df$Kids_First_Biospecimen_ID))
saveRDS(splice_data, file = file.path(data_dir, "splice_data_medullo.rds"))

# CNV (20230309)
cnv_data <- file.path(ot_dir, "20230721_release.All.gainloss.txt") %>%
  data.table::fread() %>%
  filter(BS_ID %in% hist_df$Kids_First_Biospecimen_ID)
write_tsv(cnv_data, file = file.path(data_dir, "All.gainloss_medullo.tsv"))
