# script to create subset matrices containing only MB data from OT v15

# load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# directory
ot_dir <- file.path("~/Documents/v15-Pedcan")
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_All_IntNMF")
cancer_group_of_interest <- "Medulloblastoma"

# histology (adapt)
hist_df <- file.path(ot_dir, "histologies.tsv") %>% # test ot_dir instead of data_dir
  read_tsv() %>%
  filter(short_histology %in% cancer_group_of_interest)
#write_tsv(hist_df, file = file.path(data_dir, "histologies_medullo.tsv"))


# SNV consensus + hotspots
snv_data <- data.table::fread(file.path(ot_dir, "snv-consensus-plus-hotspots.maf.tsv.gz"))

snv_data <- snv_data %>%
  filter(Tumor_Sample_Barcode %in% hist_df$Kids_First_Biospecimen_ID)
data.table::fwrite(
  snv_data,
  file = file.path(ot_dir, "mb_snv-consensus-plus-hotspots.maf.tsv.gz"),
  sep = "\t",
  na = NA,
  quote = FALSE
)

---------------------------------------------------------------------------------------
  
# Filter hist_df to get the IDs of MB, Group3
group3_ids <- hist_df %>%
  filter(molecular_subtype == "MB, Group3") %>%
  pull(Kids_First_Biospecimen_ID)

# Filter snv_data to include only the entries with Tumor_Sample_Barcode in group3_ids
snv_data_group3 <- snv_data %>%
  filter(Tumor_Sample_Barcode %in% group3_ids)

# Write the filtered data to a new file
data.table::fwrite(
  snv_data_group3,
  file = file.path(ot_dir, "mb_snv-consensus-plus-hotspots_group3.maf.tsv.gz"),
  sep = "\t",
  na = NA,
  quote = FALSE
)

kbtbd4_records <- snv_data_group3 %>%
  filter(Hugo_Symbol == "KBTBD4")
-----------------------------------------------------------------------------------------
  # CNV data
  cnv_data <- data.table::fread(file.path(ot_dir, "20231205_release.All.gainloss.txt.gz"))

cnv_data <- cnv_data %>%
  filter(BS_ID %in% hist_df$Kids_First_Biospecimen_ID)
data.table::fwrite(
  cnv_data,
  file = file.path(ot_dir, "MB.All.gainloss.txt.gz"),
  sep = "\t",
  na = NA,
  quote = FALSE
)
  # Filter hist_df to get the IDs of MB, Group3
  group3_ids <- hist_df %>%
  filter(molecular_subtype == "MB, Group3") %>%
  pull(Kids_First_Biospecimen_ID)

# Filter snv_data to include only the entries with Tumor_Sample_Barcode in group3_ids
cnv_data_group3 <- cnv_data %>%
  filter(BS_ID %in% group3_ids)

# Write the filtered data to a new file
data.table::fwrite(
  cnv_data_group3,
  file = file.path(ot_dir, "mb_cnv_grp3.tsv"),
  sep = "\t",
  na = NA,
  quote = FALSE
)

myc_records <- cnv_data_group3 %>%
  filter(gene == "MYC")
  
CDK5_records <- cnv_data_group3 %>%
  filter(gene == "CDK5")  
  
OTX2_records <- cnv_data_group3 %>%
  filter(gene == "OTX2")  
  
 
----------------------------------------------------------------------------------------

# Splice data (filtered to functional sites from Ammar Naqvi)
splice_data <- readRDS(opt$splice_file)
splice_data <- splice_data %>%
  column_to_rownames("Splice_ID") %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
saveRDS(splice_data,
        file = file.path(
          output_dir,
          "splice_events_pan_cancer_functional_filter.rds"
        ))
rm(splice_data)

# CNV gainloss file
cnv_data <- data.table::fread(opt$cnv_file)
cnv_data <- cnv_data %>%
  filter(BS_ID %in% histology_file$Kids_First_Biospecimen_ID)
data.table::fwrite(
  cnv_data,
  file = file.path(output_dir, "All.gainloss.txt.gz"),
  sep = "\t",
  na = NA,
  quote = FALSE
)
rm(cnv_data)


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
#------------------------------------------------------
# Load necessary packages
library(ggplot2)
library(ggsankey)
library(dplyr)

# Prepare the data
data <- data.frame(
  IntNMF_Cluster = c(9, 10, 10, 11, 12, 12, 6, 13),
  Generation2_Cluster = c("IV", "II", "III", "I", "V", "VI", "VII", "VIII")
)

# Create the links data frame suitable for ggsankey
links <- data %>%
  mutate(IntNMF_Cluster = as.character(IntNMF_Cluster),
         Generation2_Cluster = as.character(Generation2_Cluster)) %>%
  make_long(IntNMF_Cluster, Generation2_Cluster)

# Create the Sankey plot using ggplot2 and ggsankey
pl <- ggplot(links, aes(x = x,                        
                        next_x = next_x,                                     
                        node = node,
                        next_node = next_node,        
                        fill = factor(node),
                        label = node)) +             
  geom_sankey(flow.alpha = 0.5,                 
              node.color = "black",              
              show.legend = TRUE) +              
  geom_sankey_label(size = 3, 
                    color = "black", 
                    fill = "white") +            
  theme_void() +                                 
  theme(legend.position = "none") +
  geom_text(aes(x = "IntNMF_Cluster", y = Inf, label = "IntNMF Cluster"), 
            color = "black", size = 4, vjust = 1.5, hjust = 0.5) +
  geom_text(aes(x = "Generation2_Cluster", y = Inf, label = "Generation 2 Cluster"), 
            color = "black", size = 4, vjust = 1.5, hjust = 0.5)
# Display the plot
print(pl)

# Save the Sankey plot to a PDF file
ggsave("sankey_plot.pdf", plot = pl, device = "pdf", width = 8, height = 6)
