# script to create plots for each data modality arranged and annotated by clusters
# to determine how the clusters relate to subtypes

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(gridExtra)
  library(randomcoloR)
  library(gplots)
  library(ggplotify)
  library(corrplot)
})

# define directories
ot_dir <- file.path("~/Documents", "GitHub", "OpenPedCan-analysis", "data", "v12")
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC", "data")
analysis_dir <- file.path(root_dir, "multimodal_clustering", "SNF-CC")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = F, recursive = T)

# combine samples map with SNFCC derived clusters
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
samples_map <- samples_map %>%
  arrange(sample_id)

mm_clusters <- read_tsv(file.path(output_dir, "SNF_clusters.tsv"), col_names = FALSE)
samples_map <- samples_map %>%
  mutate(clusterID = mm_clusters$X1)

# combine SNF clusters with RNA-derived molecular subtypes
anno_file_rna <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(samples_map, by = "Kids_First_Biospecimen_ID_RNA")

  
# combine SNF clusters with methylation-derived subclass
anno_file_methyl <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, dkfz_v11_methylation_subclass, dkfz_v12_methylation_subclass) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(samples_map, by = "Kids_First_Biospecimen_ID_Methyl")

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  inner_join(anno_file_methyl)
anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)

# 1) adjusted rand index between various subtypes and SNF-CC derived clusters
sink(file.path(output_dir, "snf-cc_vs_subtypes.txt"), type = c("output"))
print("SNF-CC clusters vs RNA-derived molecular_subtypes")
mclust::adjustedRandIndex(anno_file$molecular_subtype, anno_file$clusterID)

print("SNF-CC clusters vs dkfz_v11_methylation_subclass")
mclust::adjustedRandIndex(anno_file$dkfz_v11_methylation_subclass, anno_file$clusterID)

print("SNF-CC clusters vs dkfz_v12_methylation_subclass")
mclust::adjustedRandIndex(anno_file$dkfz_v12_methylation_subclass, anno_file$clusterID)
sink()

# 2) chi-square test of independence between various subtypes and SNF-CC derived clusters
sink(file.path(output_dir, "chisq_snf_vs_subtypes.txt"), type = c("output"))
print("SNF-CC clusters vs RNA-derived molecular_subtypes")
table(anno_file$molecular_subtype, anno_file$clusterID)
chisq.test(x = anno_file$molecular_subtype, y = anno_file$clusterID)

print("Snfcc clusters vs dkfz_v11_methylation_subclass")
table(anno_file$dkfz_v11_methylation_subclass, anno_file$clusterID)
chisq.test(x = anno_file$dkfz_v11_methylation_subclass, y = anno_file$clusterID)

print("Snfcc clusters vs dkfz_v12_methylation_subclass")
table(anno_file$dkfz_v12_methylation_subclass, anno_file$clusterID)
chisq.test(x = anno_file$dkfz_v12_methylation_subclass, y = anno_file$clusterID)
sink()

# 3) generate balloon plot with at least 5 samples in a group
# SNF-CC clusters vs RNA-derived molecular subtypes
dat <- anno_file %>%
  filter(!is.na(molecular_subtype)) %>%
  group_by(molecular_subtype, clusterID)  %>%
  summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = clusterID, value = n, fill = 0) %>%
  column_to_rownames("molecular_subtype")
pdf(file = file.path(plots_dir,"MB", "Snfcc_clusters_vs_molsubtype_balloonplot.pdf"), width = 10)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "Snfcc clusters vs Molecular subtypes",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# SNF-CC clusters vs RNA-derived molecular subtypes
chisq <- chisq.test(dat)
pdf(file = file.path(plots_dir, "MB", "snfcc_clusters_vs_molsubtype_corrplot.pdf"))
corrplot(chisq$residuals, is.cor = FALSE, tl.srt = 360, tl.offset = 1, mar = c(1, 2, 1, 1),
         title = "Snfcc clusters vs Molecular subtypes")
dev.off()

# 5) generate balloon plot with at least 5 samples in a group
# SNF clusters vs methylation-derived dkfz_v11_methylation_subclass
dat <- anno_file %>%
  filter(!is.na(dkfz_v11_methylation_subclass)) %>%
  group_by(dkfz_v11_methylation_subclass, clusterID)  %>%
  summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = clusterID, value = n, fill = 0) %>%
  column_to_rownames("dkfz_v11_methylation_subclass")
pdf(file = file.path(plots_dir,"MB", "snfcc_clusters_vs_dkfz_v11_methylation_subclass_balloonplot.pdf"), width = 14)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "Snfcc clusters vs dkfz_v11_methylation_subclass",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()

# 6) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# SNF clusters vs methylation-derived dkfz_v11_methylation_subclass
chisq <- chisq.test(dat)
pdf(file = file.path(plots_dir,"MB", "snfcc_clusters_vs_dkfz_v11_methylation_subclass_corrplot.pdf"))
corrplot(chisq$residuals, is.cor = FALSE, tl.srt = 360, tl.offset = 1, mar = c(1, 2, 1, 1),
         title = "Snfcc clusters vs dkfz_v11_methylation_subclass")
dev.off()

# 7) generate balloon plot with at least 5 samples in a group
# SNF clusters vs methylation-derived dkfz_v12_methylation_subclass
dat <- anno_file %>%
  filter(!is.na(dkfz_v12_methylation_subclass)) %>%
  group_by(dkfz_v12_methylation_subclass, clusterID)  %>%
  summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = clusterID, value = n, fill = 0) %>%
  column_to_rownames("dkfz_v12_methylation_subclass")
pdf(file = file.path(plots_dir,"MB", "snfcc_clusters_vs_dkfz_v12_methylation_subclass_balloonplot.pdf"), width = 14)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "Snfcc clusters vs dkfz_v12_methylation_subclass",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()

# 8) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# SNF clusters vs methylation-derived dkfz_v12_methylation_subclass
chisq <- chisq.test(dat)
pdf(file = file.path(plots_dir,"MB", "snfcc_clusters_vs_dkfz_v12_methylation_subclass_corrplot.pdf"))
corrplot(chisq$residuals, is.cor = FALSE, tl.srt = 360, tl.offset = 1, mar = c(1, 2, 1, 1),
         title = "Snfcc clusters vs dkfz_v12_methylation_subclass")
dev.off()

.................................................................................
# 9) feature-level heatmaps using top 10 most representative features per modality
# Need to discuss with Komal and Adam on feature selection () - Not done yet
# Output heatmap
pheatmap::pheatmap(count_data,
                   annotation_col = samples_map$clusterID,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   color = bluered(256),
                   width = 12,
                   height = 8,
                   show_colnames = FALSE,  
                   fontsize_row = 5,
                   show_rownames = TRUE,
                   legend = FALSE,
                   filename = file.path(plots_dir, "gsva_heatmap.pdf"))


# 10) sample-level heatmaps : Discuss with Komal and Adam about how to get cluster specific features
# heatmaps ordered by sample cluster, by molecular subtype, methylation derived subclass - not done yet
# define distinct color palette for each annotation
set.seed(100)
cluster_palette <- distinctColorPalette(length(unique(anno_file$clusterID)))
names(cluster_palette) <- sort(unique(anno_file$clusterID))

set.seed(100)
mol_subtype_palette <- distinctColorPalette(length(unique(anno_file$molecular_subtype)))
names(mol_subtype_palette) <- sort(unique(anno_file$molecular_subtype))

set.seed(100)
dkfz_v11_methylation_subclass_palette <- distinctColorPalette(length(unique(anno_file$dkfz_v11_methylation_subclass)))
names(dkfz_v11_methylation_subclass_palette) <- sort(unique(anno_file$dkfz_v11_methylation_subclass))

set.seed(100)
dkfz_v12_methylation_subclass_palette <- distinctColorPalette(length(unique(anno_file$dkfz_v12_methylation_subclass)))
names(dkfz_v12_methylation_subclass_palette) <- sort(unique(anno_file$dkfz_v12_methylation_subclass))

# list of annotation and their colors
annots_colors <- list(cluster = cluster_palette, 
                      molecular_subtype = mol_subtype_palette,
                      dkfz_v11_methylation_subclass = dkfz_v11_methylation_subclass_palette,
                      dkfz_v12_methylation_subclass = dkfz_v12_methylation_subclass_palette)

# arrange
annots <- anno_file %>% 
  column_to_rownames('sample_id') %>%
  dplyr::select(dkfz_v12_methylation_subclass, dkfz_v11_methylation_subclass, molecular_subtype, clusterID) %>%
  dplyr::arrange(clusterID, molecular_subtype, dkfz_v11_methylation_subclass, dkfz_v12_methylation_subclass)
annots$clusterID <- as.character(annots$clusterID)

# expression
pdf(file.path(plots_dir, "sample_level_heatmaps.pdf"), width = 11, height = 12, onefile = T)
count_data <- read_tsv(file.path(input_dir, "norm_counts1.tsv")) %>% column_to_rownames()
p1 <- count_data[rownames(annots),] %>% 
  dplyr::select(colnames(expr_selected)) %>%
  t() %>%
  pheatmap::pheatmap(annotation_col = annots, 
                     annotation_colors = annots_colors,
                     fontsize = 5, cellwidth = 2, cellheight = 5,
                     show_colnames = F, 
                     cluster_cols = F, cluster_rows = T,
                     scale = "row", 
                     color = bluered(256),
                     silent = T,
                     main = paste0("Expression Data\nTop 10 features per cluster"))
p1 <- p1 %>% as.ggplot
print(p1)

# methylation
methyl_data <- read_tsv(file.path(input_dir, "methyl_data1.tsv")) %>% column_to_rownames()
p2 <- methyl_data[rownames(annots),] %>% 
  dplyr::select(colnames(methyl_selected)) %>%
  t() %>%
  pheatmap::pheatmap(annotation_col = annots, 
                     annotation_colors = annots_colors,
                     fontsize = 5, cellwidth = 2, cellheight = 5,
                     show_colnames = F, 
                     cluster_cols = F, cluster_rows = T,
                     scale = "row", 
                     color = bluered(256),
                     silent = T,
                     main = paste0("Methylation Data\nTop 10 features per cluster"))
p2 <- p2 %>% as.ggplot
print(p2)

# snv
snv_data <- read_tsv(file.path(input_dir, "snv_data.tsv")) %>% column_to_rownames()
p3 <- snv_data[rownames(annots),] %>% 
  dplyr::select(colnames(snv_selected)) %>%
  t() %>%
  pheatmap::pheatmap(annotation_col = annots,
                     annotation_colors = annots_colors,
                     fontsize = 5, cellwidth = 2, cellheight = 5,
                     show_colnames = F, 
                     cluster_cols = F, cluster_rows = T,
                     color = colorpanel(2, low = "grey95", high = "black"),
                     legend_breaks = seq(0, 1, by = 1),
                     silent = T,
                     main = paste0("Mutation Data\nTop 10 features per cluster")) 
p3 <- p3 %>% as.ggplot
print(p3)

# cnv
cnv_data <- read_tsv(file.path(input_dir, "cnv_data.tsv")) %>% column_to_rownames()
p4 <- cnv_data[rownames(annots),] %>% 
  dplyr::select(colnames(cnv_selected)) %>%
  t() %>%
  pheatmap::pheatmap(annotation_col = annots, 
                     annotation_colors = annots_colors,
                     fontsize = 5, cellwidth = 2, cellheight = 5,
                     show_colnames = F, 
                     cluster_cols = F, cluster_rows = T,
                     scale = "row", 
                     color = bluered(256),
                     silent = T, 
                     main = paste0("CNV Data\nTop 10 features per cluster")) 
p4 <- p4 %>% as.ggplot
print(p4)

# splicing
splice_data <- read_tsv(file.path(input_dir, "splice_data.tsv")) %>% column_to_rownames()
p5 <- splice_data[rownames(annots),] %>% 
  dplyr::select(colnames(splice_selected)) %>%
  t() %>%
  pheatmap::pheatmap(annotation_col = annots, 
                     annotation_colors = annots_colors,
                     fontsize = 5, cellwidth = 2, cellheight = 5,
                     show_colnames = F, 
                     cluster_cols = F, cluster_rows = T,
                     scale = "row", 
                     color = bluered(256),
                     silent = T,
                     main = paste0("Splice Data\nTop 10 features per cluster"))
p5 <- p5 %>% as.ggplot
print(p5)
dev.off()
