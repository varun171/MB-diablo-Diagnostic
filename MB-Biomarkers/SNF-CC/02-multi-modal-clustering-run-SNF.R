# run multi-modal clustering
# Script adapted from @KomalRathi MB MM-clustering using IntNMF
suppressPackageStartupMessages({
  library(SNFtool)
  library(tidyverse)
  library(CancerSubtypes)
  library(fpc)
  library(dplyr)
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

# read data
count_data <- read_tsv(file.path(input_dir, "norm_counts1.tsv")) %>% column_to_rownames() %>% t()
methyl_data <- read_tsv(file.path(input_dir, "methyl_data1.tsv")) %>% column_to_rownames() %>% t()
snv_data <- read_tsv(file.path(input_dir, "snv_data.tsv")) %>% column_to_rownames() %>% t()
cnv_data <- read_tsv(file.path(input_dir, "cnv_data.tsv")) %>% column_to_rownames()%>% t()
splice_data <- read_tsv(file.path(input_dir, "splice_data1.tsv")) %>% column_to_rownames()%>% t()
samples_map <- read_tsv(file.path(input_dir, "samples_map.tsv")) 


# combine into a list 

dat <- list(as.matrix(count_data), as.matrix(methyl_data), as.matrix(snv_data), as.matrix(cnv_data), as.matrix(splice_data))

########################################################## did not work
##Error in silhouette.default(clustering, dmatrix = dmat) :
#'dmatrix' is not a dissimilarity matrix compatible to 'x'
#'
wt = if(is.list(dat)) rep(1,length(dat)) else 1
# Execute SNFCC
source(file.path(analysis_dir, "utils", "run_clusterstats_SNF.R"))

snf_output <- run_clusterstats(dat = dat, 
                               wt = wt, 
                               output_dir = output_dir, clusterNum = 20)

###############################################################################
#optimize both cluster number and K with nested loop
# Create an empty list to store the results for different alpha and clusterNum values
results_list <- list()

# Iterate through clusterNum values from 4 to 12 with an increment of 1
for (clusterNum in 4:20) {
  # Iterate through alpha values from 0.3 to 0.8 with an increment of 0.1
  for (alpha in seq(0.3, 0.8, by = 0.1)) {
    # Execute SNF with the current alpha and clusterNum values
    result_SNF <- ExecuteSNF.CC(dat, clusterNum = clusterNum, K = 20, alpha = alpha, t = 50, 
                                maxK = 20, pItem= 0.8, reps=500, title = "MB", 
                                plot = "png", finalLinkage ="average")
    
    # Store the result in the results_list using a unique key for each combination of alpha and clusterNum
    results_list[[paste0("alpha_", alpha, "_clusterNum_", clusterNum)]] <- result_SNF
  }
}

####### Average Sil width
# Iterate through each element of the results_list
for (key in names(results_list)) {
  result_SNF <- results_list[[key]]
  
  # Calculate silhouette similarity matrix
  sil <- silhouette_SimilarityMatrix(result_SNF$group, result_SNF$distanceMatrix)
  
  # Create a plot for the current result
  pdf(file = file.path(plots_dir, paste0("snf_silhouette_plot_", key, ".pdf")), 
      width = 10, height = 10)
  plot(sil, main = paste0("Silhouette Plot (", key, ")"))
  dev.off()
}
saveRDS(result_SNF, file = file.path(output_dir, "result_SNF.rds"))

# Run cluster again with optimum parameters and save the cluster groups
Final_SNF <- ExecuteSNF.CC(dat, clusterNum = 6, K = 20, alpha = 0.5, t = 50, 
                           maxK = 15, pItem= 0.8, reps=500, title = "Snfcc", 
                           plot = "pdf", finalLinkage ="average")

sil=silhouette_SimilarityMatrix(Final_SNF$group, Final_SNF$distanceMatrix)
summary(sil)

sil1=silhouette(Final_SNF$group, Final_SNF$distanceMatrix)
distance <- 1 - Final_SNF$distanceMatrix
sil1=silhouette(Final_SNF$group, distance)
summary(sil1)

pdf(file = file.path(plots_dir, "MB", "Final_snf_silhouette_plot.pdf"), 
    width = 10, height = 10)
plot(sil)
dev.off()

SNF_cluster <- Final_SNF$group

# Save the data frame as TSV
write.table(SNF_cluster, file = file.path(output_dir,"SNF_clusters.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


