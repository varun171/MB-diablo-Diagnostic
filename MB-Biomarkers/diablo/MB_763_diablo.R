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
nmf_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/results/MB763")
output_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/results/MB763/intnmf")
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/plots/MB763/intnmf")
dir.create(plots_dir, showWarnings = F, recursive = T)


# read data
count_data <-
  read_tsv(file.path(data_dir, "norm_763counts.tsv")) %>% column_to_rownames()
methyl_data <-
  read_tsv(file.path(data_dir, "methyl_763data.tsv")) %>% column_to_rownames() 

samples_map <- read_tsv(file.path(data_dir, "samples_763map.tsv"))
nmf_clusters <- read_tsv(file.path(nmf_dir, "intnmf_clusters.tsv")) 
nmf_clusters <- nmf_clusters %>%
  inner_join(samples_map)

# Create the list of datasets
X <-
  list(
    mRNA = as.matrix(count_data),
    methyl = as.matrix(methyl_data))


# Outcome
nmf_clusters$intnmf_cluster <- paste("cluster", nmf_clusters$intnmf_cluster, sep = "")
#Y <- as.factor(nmf_clusters$intnmf_cluster)
Y <- as.factor(nmf_clusters$subtype)
# derive weights
res1.pls <- pls(X$mRNA, X$methyl, ncomp = 1)
cor(res1.pls$variates$X, res1.pls$variates$Y)

# create the design

design <- matrix(
  0.1,
  ncol = length(X),
  nrow = length(X),
  dimnames = list(names(X), names(X))
)

design["mRNA", "methyl"] <- 0.98
design["methyl", "mRNA"] <- 0.98

diag(design) <- 0


# Select the number of components

diablo.MB <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123)
perf.diablo.MB = perf(diablo.MB,
                      validation = 'Mfold',
                      folds = 10,
                      nrepeat = 10)

pdf(file.path(plots_dir, "perf_diablo_MB.pdf"))
plot(perf.diablo.MB)
dev.off()

ncomp = perf.diablo.MB$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo.MB$choice.ncomp$WeightedVote 

# select the number of variables (cab be optimized further with more variables and nrepeat)

set.seed(123)
test.keepX <-
  list (
    count_data = c(5:9, seq(10, 18, 2), seq(20, 30, 5)),
    methyl_data = c(5:9, seq(10, 18, 2), seq(20, 30, 5))
  )

tune.diablo.MB <- tune.block.splsda(
  X,
  Y,
  ncomp = ncomp,
  test.keepX = test.keepX,
  design = design,
  validation = 'Mfold',
  folds = 10,
  nrepeat = 10,
  BPPARAM = BiocParallel::SnowParam(workers = 8),
  dist = "centroids.dist"
)

list.keepX <- tune.diablo.MB$choice.keepX


# Final Model

diablo_MB <- block.splsda(X,
                          Y,
                          ncomp = ncomp,
                          keepX = list.keepX,
                          design = design)

saveRDS(diablo_MB, file = file.path(output_dir, "diablo_MB.rds"))
#diablo_MB <- readRDS(file = file.path(output_dir, "diablo_MB.rds"))

# QC Plots

pdf(file.path(plots_dir, "diagnostic_plot.pdf"))
plotDiablo(diablo_MB, ncomp = 1)
dev.off()

pdf(file.path(plots_dir, "sample_plot.pdf"))
plotIndiv(diablo_MB,
          ind.names = FALSE,
          legend = TRUE,
          title = 'MB, DIABLO comp 1 - 2')
dev.off()

pdf(file.path(plots_dir, "arrow_plot.pdf"))
plotArrow(diablo_MB,
          ind.names = FALSE,
          legend = TRUE,
          title = 'MB, DIABLO comp 1 - 2')
dev.off()

# variable Plots

pdf(file.path(plots_dir, "correlation_plot.pdf"))
plotVar(
  diablo_MB,
  var.names = FALSE,
  style = 'graphics',
  legend = TRUE,
  pch = c(16, 17),
  cex = c(2, 2),
  col = c('darkorchid', 'brown1'),
  title = 'MB, DIABLO comp 1 - 2'
)
dev.off()

pdf(file.path(plots_dir, "circos_plot.pdf"),
    width = 11,
    height = 8)
circosPlot(
  diablo_MB,
  cutoff = 0.7,
  line = TRUE,
  color.blocks = c('darkorchid', 'brown1'),
  color.cor = c("chocolate3", "grey20"),
  size.labels = 1.5
)
dev.off()

pdf(file.path(plots_dir, "circos_plot_com1.pdf"),
    width = 11,
    height = 8)
circosPlot(
  diablo_MB,
  comp = 1,
  cutoff = 0.85,
  line = TRUE,
  color.blocks = c('darkorchid', 'brown1'),
  color.cor = c("chocolate3", "grey20"),
  size.labels = 1.5
)
dev.off()


pdf(file.path(plots_dir, "network_plot_com.pdf"),
    width = 30,
    height = 20)
network(
  diablo_MB,
  blocks = c(1, 2),
  graph.scale = 0.3,
  size.node = 0.2,
  cutoff = 0.85,
  cex.node.name = 0.5,
  color.node = c('darkorchid', 'brown1')
)
dev.off()

pdf(file.path(plots_dir, "loadings_plot_1_test.pdf"),
    width = 14,
    height = 4)
plotLoadings(diablo_MB,
             comp = 1,
             contrib = 'max',
             method = 'median',
             legend.color = rainbow(nlevels(Y)))
dev.off()

pdf(file.path(plots_dir, "loadings_plot_2.pdf"),
    width = 14,
    height = 4)
plotLoadings(diablo_MB,
             comp = 2,
             contrib = 'max',
             method = 'median',
             legend.color = rainbow(12))
dev.off()

pdf(file.path(plots_dir, "loadings_plot_3.pdf"),
    width = 14,
    height = 4)
plotLoadings(diablo_MB,
             comp = 3,
             contrib = 'max',
             method = 'median',
             legend.color = rainbow(12))
dev.off()

pdf(file.path(plots_dir, "loadings_plot_4.pdf"),
    width = 14,
    height = 4)
plotLoadings(diablo_MB,
             comp = 4,
             contrib = 'max',
             method = 'median',
             legend.color = rainbow(12))
dev.off()

pdf(file.path(plots_dir, "loadings_plot_5.pdf"),
    width = 14,
    height = 4)
plotLoadings(diablo_MB,
             comp = 5,
             contrib = 'max',
             method = 'median',
             legend.color = rainbow(12))
dev.off()



pdf(file.path(plots_dir, "image_plot.pdf"),
    width = 18,
    height = 12)
cimDiablo(
  diablo_MB,
  color.blocks = c('darkorchid', 'brown1'),
  comp = 1,
  margin = c(8, 20),
  legend.position = "right"
)
dev.off()

# Model Performance

set.seed(123)
perf.diablo.MB = perf(
  diablo.MB,
  validation = 'Mfold',
  folds = 3,
  nrepeat = 10,
  dist = 'centroids.dist'
)

pdf(file.path(plots_dir, "auroc_mRNA_com1.pdf"),
    width = 18,
    height = 12)
auc.diablo.MB <-
  auroc(diablo.MB,
        roc.block = "mRNA",
        roc.comp = 1,
        print = FALSE)
dev.off()

pdf(
  file.path(plots_dir, "auroc_methyl_com1.pdf"),
  width = 18,
  height = 12
)
auc.diablo.MB <-
  auroc(diablo.MB,
        roc.block = "methyl",
        roc.comp = 1,
        print = FALSE)
dev.off()

