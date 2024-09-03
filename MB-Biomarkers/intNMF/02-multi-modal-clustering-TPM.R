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
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = F, recursive = T)

# source function
#source(file.path(analysis_dir, "utils", "nmf_opt_k.R"))
#source(file.path(analysis_dir, "utils", "run_clusterstats.R"))
source(file.path(analysis_dir, "utils", "run_NMF_vk.R"))
# read data
tpm_data <- read_tsv(file.path(data_dir, "norm_tpm.tsv")) %>% column_to_rownames() 
methyl_data <- read_tsv(file.path(data_dir, "methyl_mdata.tsv")) %>% column_to_rownames() 

## Also rescale the datasets so that they are comparable.
if (!all(tpm_data>=0)) tpm_data <- pmax(tpm_data + abs(min(tpm_data)), .Machine$double.eps)
tpm_data <- tpm_data/max(tpm_data)
if (!all(methyl_data>=0)) methyl_data <- pmax(methyl_data + abs(min(methyl_data)), .Machine$double.eps)
methyl_data <- methyl_data/min(methyl_data)

# combine into a list (took almost 6 hours)
#dat <- list(as.matrix(dat1), as.matrix(dat2), as.matrix(dat3), as.matrix(dat4), as.matrix(dat5))
dat <- list((tpm_data), (methyl_data))
#dat <- list(as.matrix(dat1), as.matrix(dat2))

source(file.path(analysis_dir, "utils", "run_NMF_vk.R"))

wt = if(is.list(dat)) rep(1,length(dat)) else 1

nmf_output <- run_clusterstats(dat = dat, 
                               wt = wt, 
                               output_dir = output_dir, 
                               k_value = 8)
# output_fname <- file.path(output_dir, "cpi_output.rds")
# if(!file.exists(output_fname)){
#   cpi_output <- nmf_opt_k(dat = dat, 
#                           wt = wt,
#                           n.runs = 30, # default
#                           maxiter = 100, # default
#                           st.count = 10, # default
#                           n.fold = 5, # default
#                           k.range = 2:15, 
#                           result = TRUE, 
#                           make.plot = TRUE,
#                           progress = TRUE,
#                           output_dir = plots_dir)
#   saveRDS(cpi_output, file = output_fname)
# } 

#nmf_opt <- nmf.opt.k(dat = dat, n.runs = 30, n.fold = 5, k.range = 2:8, result = TRUE,
#          make.plot = TRUE, progress = TRUE, st.count = 10, maxiter = 100,
#          wt=if(is.list(dat)) rep(1,length(dat)) else 1)


#nmf_output <- nmf.mnnals(dat = dat, k = 6, maxiter = 200, st.count = 10, n.ini = 10, ini.nndsvd = TRUE, seed = TRUE,wt=if(is.list(dat)) rep(1,length(dat)) else 1)
#saveRDS(nmf_output, file.path(output_dir, "intnmf_output.rds"))

nmf_output <- readRDS(file.path(output_dir, "nmf_best_fit.rds"))


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
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(output_dir, "intnmf_clusters.tsv"))
mm_clusters <- mm_clusters %>%
  inner_join(samples_map)

# combine IntNMF clusters with RNA-derived molecular subtypes
anno_file_rna <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# 3) generate balloon plot with at least 5 samples in a group
# IntNMF clusters vs RNA-derived molecular subtypes
dat <- anno_file_rna %>%
  filter(!is.na(molecular_subtype)) %>%
  group_by(molecular_subtype, intnmf_cluster)  %>%
  dplyr::summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = intnmf_cluster, value = n, fill = 0) %>%
  column_to_rownames("molecular_subtype")
pdf(file = file.path(plots_dir, "intnmf_clusters_vs_molsubtype_balloonplot.pdf"), width = 10)
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "IntNMF clusters vs Molecular subtypes",
            xlab = "", ylab = "",
            label = TRUE, 
            show.margins = FALSE)
dev.off()

# 4) survival curve
# format survival data
surv_data <- anno_file_rna
surv_data$OS_status <- ifelse(surv_data$OS_status == "DECEASED", 1, 0)

surv_data$intnmf_cluster <- factor(surv_data$intnmf_cluster, levels = sort(as.numeric(unique(surv_data$intnmf_cluster))))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ intnmf_cluster, data = surv_data)
pdf(file = file.path(plots_dir, "survival_IntNMF_clusters.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "Survival stratified by Newly derived clusters",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.coord = c(6000, 0.75),
                pval.method = TRUE,
                pval.method.coord = c(6000, 0.80),
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")
# guides(colour = guide_legend(ncol = 4)) 
# Modify the plot component to add guides
p$plot <- p$plot + guides(colour = guide_legend(ncol = 4))
print(p)
dev.off()

# 5) survival curves stratified by RNA-derived molecular subtype 
surv_data$molecular_subtype <- factor(surv_data$molecular_subtype, levels = sort(unique(surv_data$molecular_subtype)))
fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ molecular_subtype, data = surv_data)
pdf(file = file.path(plots_dir, "survival_molsubtype.pdf"), height = 8, width = 10, onefile = FALSE)
p <- ggsurvplot(fit,
                title = "Survival stratified by RNA-derived molecular subtypes",
                data = surv_data,
                font.x = c(12, face = "bold"),
                font.tickslab = c(10, face = "bold"),
                font.y = c(12, face = "bold"),
                font.legend = c(10, "bold"),
                pval = TRUE, 
                pval.coord = c(6000, 0.75),
                pval.method = TRUE,
                pval.method.coord = c(6000, 0.80),
                ggtheme = theme_minimal(),
                linetype = "strata",
                legend = "bottom")  

# Modify the plot component to add guides
p$plot <- p$plot + guides(colour = guide_legend(ncol = 2))
#guides(colour = guide_legend(ncol = 2)) 
print(p)
dev.off()

# 6) chi-square test of independence between various subtypes and IntNMF derived clusters
sink(file.path(output_dir, "chisq_intnmf_vs_subtypes.txt"), type = c("output"))
print("IntNMF clusters vs RNA-derived molecular_subtypes")
table(anno_file_rna$molecular_subtype, anno_file_rna$intnmf_cluster)
chisq.test(x = anno_file_rna$molecular_subtype, y = anno_file_rna$intnmf_cluster)
sink()

# 7) feature-level heatmaps using top 10 most representative features per modality
# read nmf output for best fit
nmf_output <- readRDS(file.path(output_dir, "nmf_best_fit.rds"))

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

# snv 
snv_selected <- nmf_output$H$H3
snv_selected <- snv_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  slice_head(n = 20) %>% 
  spread(key = "feature", value = "value", fill = 0) %>%
  column_to_rownames("cluster")
p3 <- pheatmap::pheatmap(mat = t(snv_selected), 
                         fontsize = 5, cellwidth = 8, cellheight = 5,
                         scale = "none", 
                         angle_col = 0,
                         silent = T,
                         main = paste0("Mutation Data\nTop 20 features per cluster"))

# cnv
cnv_selected <- nmf_output$H$H4
cnv_selected <- cnv_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  slice_head(n = 20) %>% 
  spread(key = "feature", value = "value", fill = 0) %>%
  column_to_rownames("cluster")
p4 <- pheatmap::pheatmap(mat = t(cnv_selected), 
                         fontsize = 5, cellwidth = 8, cellheight = 5,
                         scale = "row", 
                         angle_col = 0,
                         silent = T,
                         main = paste0("CNV Data\nTop 20 features per cluster"))

# splicing
splice_selected <- nmf_output$H$H5
splice_selected <- splice_selected %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  gather(key = "feature", value = "value", -c(cluster)) %>% 
  group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  slice_head(n = 20) %>% 
  spread(key = "feature", value = "value", fill = 0) %>%
  column_to_rownames("cluster")
p5 <- pheatmap::pheatmap(mat = t(splice_selected), 
                         fontsize = 5, cellwidth = 8, cellheight = 5,
                         scale = "row",
                         angle_col = 0,
                         silent = T,
                         main = paste0("Splice Data\nTop 20 features per cluster"))

# combine all in one
pdf(file.path(plots_dir, "feature_level_heatmaps.pdf"), width = 12, height = 23)
grid.arrange(arrangeGrob(grobs = list(p1$gtable, p2$gtable, p3$gtable, p4$gtable, p5$gtable), ncol = 3))
dev.off()

# 8) sample-level heatmaps
# heatmaps ordered by sample cluster, by molecular subtype, methylation derived subclass
# define distinct color palette for each annotation
set.seed(100)
cluster_palette <- distinctColorPalette(length(unique(anno_file_rna$intnmf_cluster)))
names(cluster_palette) <- sort(unique(anno_file_rna$intnmf_cluster))

set.seed(100)
mol_subtype_palette <- distinctColorPalette(length(unique(anno_file_rna$molecular_subtype)))
names(mol_subtype_palette) <- sort(unique(anno_file_rna$molecular_subtype))

# list of annotation and their colors
annots_colors <- list(cluster = cluster_palette, 
                      molecular_subtype = mol_subtype_palette)

# arrange
annots <- anno_file_rna %>% 
  column_to_rownames('sample_id') %>%
  dplyr::select(molecular_subtype, intnmf_cluster) %>%
  dplyr::arrange(intnmf_cluster, molecular_subtype)
annots$intnmf_cluster <- as.character(annots$intnmf_cluster)

# expression
pdf(file.path(plots_dir, "sample_level_heatmaps.pdf"), width = 11, height = 12, onefile = T)
count_data <- read_tsv(file.path(input_dir, "norm_counts.tsv")) %>% column_to_rownames()
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
methyl_data <- read_tsv(file.path(input_dir, "methyl_data.tsv")) %>% column_to_rownames()
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

