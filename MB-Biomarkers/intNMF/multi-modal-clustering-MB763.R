# run multi-modal clustering
suppressPackageStartupMessages({
  library(IntNMF)
  library(tidyverse)
  library(pheatmap)
  library(gridExtra)
  library(randomcoloR)
  library(gplots)
  library(ggplotify)
  library(corrplot)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# define directories
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results/MB763")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots/MB763")
dir.create(plots_dir, showWarnings = F, recursive = T)

# read data
count_data <- read_tsv(file.path(data_dir, "norm_763counts.tsv")) %>% column_to_rownames()
methyl_data <- read_tsv(file.path(data_dir, "methyl_763data.tsv")) %>% column_to_rownames()


# combine into a list (took almost 6 hours)
dat <- list(as.matrix(count_data), as.matrix(methyl_data))
wt = if(is.list(dat)) rep(1,length(dat)) else 1

nmf_opt <- nmf.opt.k(dat = dat, n.runs = 30, n.fold = 5, k.range = 2:8, result = TRUE,
                     make.plot = TRUE, progress = TRUE, st.count = 10, maxiter = 100,
                     wt=if(is.list(dat)) rep(1,length(dat)) else 1)

nmf_output <- nmf.mnnals(dat = dat, k = 4, maxiter = 200, st.count = 10, n.ini = 10, ini.nndsvd = TRUE, seed = TRUE,wt=if(is.list(dat)) rep(1,length(dat)) else 1)
saveRDS(nmf_output, file.path(output_dir, "intnmf_output.rds"))

# 1) ConsensusMatPlot
# Given the integrative NMF fit object, the function creates image plot of the consensus matrix ordered
# according to clusters groups. Cleaner block structure indicates stronger clusters
pdf(file = file.path(plots_dir, "intnmf_consensus_plot.pdf"), width = 10, height = 10, onefile = F)
ConsensusMatPlot(fit = nmf_output, rowLab = TRUE, colLab = TRUE)
dev.off()

# 2) SilhouettePlot 
# Silhouette width plot is returned together with mean silhouette width for each group, overall silhouette width and summary statistics.
pdf(file = file.path(plots_dir, "intnmf_silhouette_plot.pdf"), width = 10, height = 10, onefile = F)
SilhouettePlot(fit = nmf_output, cluster.col = NULL)
dev.off()

# 3) output clusters with per sample  
df <- data.frame(sample_id = names(nmf_output$clusters), intnmf_cluster = nmf_output$clusters)
write_tsv(df, file = file.path(output_dir, "intnmf_clusters.tsv"))

#4) Concordance Plot
# combine samples map with IntNMF derived clusters
samples_map <- read_tsv(file.path(data_dir, "samples_763map.tsv"))
mm_clusters <- read_tsv(file.path(output_dir, "intnmf_clusters.tsv"))
mm_clusters <- mm_clusters %>%
  inner_join(samples_map)


# 3) generate balloon plot with at least 5 samples in a group
# IntNMF clusters vs RNA-derived molecular subtypes
dat <- mm_clusters %>%
  filter(!is.na(subgroup)) %>%
  group_by(subgroup, intnmf_cluster)  %>%
  dplyr::summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = intnmf_cluster, value = n, fill = 0) %>%
  column_to_rownames("subgroup")
pdf(file = file.path(plots_dir, "intnmf_clusters_vs_763subgroup_balloonplot.pdf"), width = 10)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "IntNMF clusters vs Molecular subtypes",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()

# 4) generate balloon plot with at least 5 samples in a group
# IntNMF clusters vs RNA-derived molecular subtypes
dat <- mm_clusters %>%
  filter(!is.na(subtype)) %>%
  group_by(subtype, intnmf_cluster)  %>%
  dplyr::summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = intnmf_cluster, value = n, fill = 0) %>%
  column_to_rownames("subtype")
pdf(file = file.path(plots_dir, "intnmf_clusters_vs_763subtype_balloonplot.pdf"), width = 10)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "IntNMF clusters vs Molecular subtypes",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()


# 5) chi-square test of independence between various subtypes and IntNMF derived clusters
sink(file.path(output_dir, "chisq_intnmf_vs_763subgroup.txt"), type = c("output"))
print("IntNMF clusters vs 763 molecular_subgroup")
table(mm_clusters$subgroup, mm_clusters$intnmf_cluster)
chisq.test(x = mm_clusters$subgroup, y = mm_clusters$intnmf_cluster)
sink()

# 6) chi-square test of independence between various subtypes and IntNMF derived clusters
sink(file.path(output_dir, "chisq_intnmf_vs_763subtype.txt"), type = c("output"))
print("IntNMF clusters vs 763 molecular_subtype")
table(mm_clusters$subtype, mm_clusters$intnmf_cluster)
chisq.test(x = mm_clusters$subtype, y = mm_clusters$intnmf_cluster)
sink()

# 7) feature-level heatmaps using top 10 most representative features per modality
# read nmf output for best fit
nmf_output <- readRDS(file.path(output_dir, "intnmf_output.rds"))

# expression
expr_selected <- nmf_output$H$H1
expr_selected <- expr_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  slice_head(n = 20) %>% 
  spread(key = "feature", value = "value", fill = 0) %>%
  column_to_rownames("cluster")
p1 <- pheatmap::pheatmap(mat = t(expr_selected), 
                         fontsize = 5, cellwidth = 8, cellheight = 5,
                         scale = "row",
                         angle_col = 0,
                         silent = T,
                         main = paste0("Expression Data\nTop 20 features per cluster"))

# methylation
methyl_selected <- nmf_output$H$H2
methyl_selected <- methyl_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  slice_head(n = 20) %>% 
  spread(key = "feature", value = "value", fill = 0) %>%
  column_to_rownames("cluster")
p2 <- pheatmap::pheatmap(mat = t(methyl_selected), 
                         fontsize = 5, cellwidth = 8, cellheight = 5,
                         scale = "row",
                         angle_col = 0,
                         silent = T,
                         main = paste0("Methylation Data\nTop 20 features per cluster"))


# combine all in one
pdf(file.path(plots_dir, "feature_level_heatmaps.pdf"), width = 12, height = 23)
grid.arrange(arrangeGrob(grobs = list(p1$gtable, p2$gtable), ncol = 3))
dev.off()

# 8) sample-level heatmaps
# heatmaps ordered by sample cluster, by molecular subtype, methylation derived subclass
# define distinct color palette for each annotation
set.seed(100)
cluster_palette <- distinctColorPalette(length(unique(mm_clusters$intnmf_cluster)))
names(cluster_palette) <- sort(unique(mm_clusters$intnmf_cluster))

set.seed(100)
mol_subtype_palette <- distinctColorPalette(length(unique(mm_clusters$subgroup)))
names(mol_subtype_palette) <- sort(unique(mm_clusters$subgroup))

# list of annotation and their colors
annots_colors <- list(cluster = cluster_palette, 
                      molecular_subtype = mol_subtype_palette)

# arrange
annots <- mm_clusters %>% 
  column_to_rownames('sample_id') %>%
  dplyr::select(subgroup, intnmf_cluster) %>%
  dplyr::arrange(intnmf_cluster, subgroup)
annots$intnmf_cluster <- as.character(annots$intnmf_cluster)

# expression
pdf(file.path(plots_dir, "sample_level_heatmaps.pdf"), width = 15, height = 15, onefile = T)
count_data <- read_tsv(file.path(data_dir, "norm_763counts.tsv")) %>% column_to_rownames()
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
methyl_data <- read_tsv(file.path(data_dir, "methyl_763data.tsv")) %>% column_to_rownames()
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
dev.off()