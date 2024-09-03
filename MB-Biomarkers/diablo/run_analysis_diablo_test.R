# diable analysis using multiomics data
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(mixOmics)
})

# parse command line options
option_list <- list(
    make_option(c("--cluster_file"), type = "character",
              help = "path to cluster annotation file"),
    make_option(c("--output_dir"), type = "character",
              help = "path to results directory"),
    make_option(c("--plots_dir"), type = "character",
                help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)


# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T) 

# Read clusterfile
mm_clusters <- read_tsv(opt$cluster_file)

# define directories
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
input_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/input")
#output_dir <- file.path("multimodal_clustering/diablo/results/intnmf")
#plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/diablo/plots/intnmf")
#mm_clusters <- read_tsv(file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/intNMF/results/intnmf_clusters.tsv"))

# read data
count_data <- read_tsv(file.path(input_dir, "norm_counts.tsv")) %>% column_to_rownames() 
methyl_data <- read_tsv(file.path(input_dir, "methyl_data1.tsv")) %>% column_to_rownames() 
#snv_data <- read_tsv(file.path(input_dir, "snv_data.tsv")) %>% column_to_rownames()
#cnv_data <- read_tsv(file.path(input_dir, "cnv_data.tsv")) %>% column_to_rownames()
splice_data <- read_tsv(file.path(input_dir, "splice_data1.tsv")) %>% column_to_rownames() 
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv")) 
# merge clusters
nmf_clusters <- mm_clusters %>% inner_join(samples_map)
nmf_anno <- file.path(input_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(nmf_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# Create the list of datasets
#X <- list(mRNA = as.matrix(count_data), methyl = as.matrix(methyl_data), mutation = as.matrix(snv_data), copy_number = as.matrix(cnv_data), splice_isoforms = as.matrix(splice_data))
X <- list(mRNA = as.matrix(count_data), methyl = as.matrix(methyl_data), splice_isoforms = as.matrix(splice_data))


# Outcome
#nmf_clusters$intnmf_cluster <- paste("cluster", nmf_clusters$intnmf_cluster, sep = "")
Y <- as.factor(nmf_anno$molecular_subtype)

# derive weights
res1.pls <- pls(X$mRNA, X$methyl, ncomp = 1)
cor(res1.pls$variates$X, res1.pls$variates$Y)

res2.pls <- pls(X$mRNA, X$splice_isoforms, ncomp = 1)
cor(res2.pls$variates$X, res2.pls$variates$Y)

res3.pls <- pls(X$methyl, X$splice_isoforms, ncomp = 1)
cor(res3.pls$variates$X, res3.pls$variates$Y)

# create the design

design <- matrix(0.1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))

design["mRNA", "methyl"] <- 0.9
design["methyl", "mRNA"] <- 0.9  

design["mRNA", "splice_isoforms"] <- 0.9
design["splice_isoforms", "mRNA"] <- 0.9  

design["methyl", "splice_isoforms"] <- 0.8
design["splice_isoforms", "methyl"] <- 0.8  

diag(design) <- 0


# Select the number of components

diablo.MB <- block.plsda(X, Y, ncomp = 5, design = design)

set.seed(123) 
perf.diablo.MB = perf(diablo.MB, validation = 'Mfold', folds = 3, nrepeat = 10)

ncomp <- perf.diablo.MB$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]

# select the number of variables (cab be optimized further with more variables and nrepeat)

set.seed(123) 
test.keepX <- list(count_data = c(seq(10, 25, 5)),
                   methyl_data = c(seq(10, 25, 5)),
                   splice_data = c(seq(10, 25, 5)))

tune.diablo.MB <- tune.block.splsda(X, Y, ncomp = 4, 
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = 3, nrepeat = 1, 
                                      BPPARAM = BiocParallel::SnowParam(workers = 2),
                                      dist = "mahalanobis.dist")

list.keepX <- tune.diablo.MB$choice.keepX


# Final Model

diablo_MB <- block.splsda(X, Y, ncomp = 4, 
                            keepX = list.keepX, design = design)

saveRDS(diablo_MB, file = file.path(plots_dir, "diablo_MB.rds"))

#selectVar(diablo_MB, block = 'mRNA', comp = 4)

# QC Plots

pdf(file.path(plots_dir, "diagnostic_plot.pdf"))
plotDiablo(diablo_MB, ncomp = 1) 
dev.off()

pdf(file.path(plots_dir, "sample_plot.pdf"))
plotIndiv(diablo_MB, ind.names = FALSE, legend = TRUE, 
          title = 'MB, DIABLO comp 1 - 2')
dev.off()

pdf(file.path(plots_dir, "arrow_plot.pdf"))
plotArrow(diablo_MB, ind.names = FALSE, legend = TRUE, 
          title = 'MB, DIABLO comp 1 - 2')
dev.off()

# variable Plots (splice data col_name needs to be reduced)

pdf(file.path(plots_dir, "correlation_plot.pdf"))
plotVar(diablo_MB, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2, 2, 2), 
        col = c('darkorchid', 'brown1', 'lightgreen'),
        title = 'MB, DIABLO comp 1 - 2')
dev.off()

pdf(file.path(plots_dir, "circos_plot.pdf"), width = 11, height = 8)
circosPlot(diablo_MB, cutoff = 0.7, line = TRUE, 
           color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
dev.off()

pdf(file.path(plots_dir, "network_plot.pdf"), width = 30, height = 20)
network(diablo_MB, blocks = c(1,2,3), 
        cutoff = 0.4,
        color.node = c('darkorchid', 'brown1', 'lightgreen'))
dev.off()

pdf(file.path(plots_dir, "loadings_plot_1.pdf"), width = 14, height = 8)
plotLoadings(diablo_MB, comp = 1, contrib = 'max', method = 'median')
dev.off()

pdf(file.path(plots_dir, "loadings_plot_2.pdf"), width = 14, height = 8)
plotLoadings(diablo_MB, comp = 2, contrib = 'max', method = 'median')
dev.off()

pdf(file.path(plots_dir, "loadings_plot_3.pdf"), width = 14, height = 8)
plotLoadings(diablo_MB, comp = 3, contrib = 'max', method = 'median')
dev.off()

pdf(file.path(plots_dir, "loadings_plot_4.pdf"), width = 14, height = 8)
plotLoadings(diablo_MB, comp = 4, contrib = 'max', method = 'median')
dev.off()


pdf(file.path(plots_dir, "image_plot.pdf"), width = 18, height = 12)
cimDiablo(diablo_MB, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8,20), legend.position = "right")
dev.off()
        
# Model Performance

set.seed(123) 
perf.diablo.MB = perf(diablo.MB, validation = 'Mfold', folds = 3, nrepeat = 10, dist = 'mahalanobis.dist')

pdf(file.path(plots_dir, "auroc_mRNA_com1.pdf"), width = 18, height = 12)
auc.diablo.MB <- auroc(diablo.MB, roc.block = "mRNA", roc.comp = 1, print = FALSE)
dev.off()

pdf(file.path(plots_dir, "auroc_methyl_com1.pdf"), width = 18, height = 12)
auc.diablo.MB <- auroc(diablo.MB, roc.block = "methyl", roc.comp = 1, print = FALSE)
dev.off()

pdf(file.path(plots_dir, "auroc_splice_com1.pdf"), width = 18, height = 12)
auc.diablo.MB <- auroc(diablo.MB, roc.block = "splice_isoforms", roc.comp = 1, print = FALSE)
dev.off()