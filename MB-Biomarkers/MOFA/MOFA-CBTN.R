# One sungroup or subtype for features and pathways identification 
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

# define directories
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
analysis_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MOFA")

# input directory
input_dir <- file.path(analysis_dir, "input")

# output directory
output_dir <- file.path(analysis_dir, "results")
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = F, recursive = T)

# read data
count_data <- read_tsv(file.path(data_dir, "norm_counts.tsv")) %>% column_to_rownames() %>% t()
methyl_data <- read_tsv(file.path(data_dir, "methyl_data.tsv")) %>% column_to_rownames() %>% t()

# combine into a list 
dat <- list(count_data = as.matrix(count_data), methyl_data = as.matrix(methyl_data))

# combine samples map with IntNMF derived clusters
samples_map <- read_tsv(file.path(data_dir, "samples_map.tsv"))
mm_clusters <- read_tsv(file.path(data_dir, "intnmf_clusters.tsv"))
mm_clusters <- mm_clusters %>%
  inner_join(samples_map)

# combine IntNMF clusters with RNA-derived molecular subtypes
anno_file_rna <- file.path(data_dir, "histologies_medullo.tsv") %>%
  read_tsv() %>%
  dplyr::select(Kids_First_Biospecimen_ID, OS_days, OS_status, molecular_subtype) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA") %>%
  dplyr::filter(intnmf_cluster == 2)

colnames(anno_file_rna)[colnames(anno_file_rna) == "sample_id"] <- "sample"

# Filter count_data
count_data <- count_data[, colnames(count_data) %in% anno_file_rna$sample]

# Filter methyl_data
methyl_data <- methyl_data[, colnames(methyl_data) %in% anno_file_rna$sample]

# Create MOFA object

MOFAobject <- create_mofa(dat)

#plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 5

model_opts

# train the model
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 40

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)

# Train the MOFA model
MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

# add metadata to the mofa object
samples_metadata(MOFAobject) <- anno_file_rna

# correlation between the factors
plot_factor_cor(MOFAobject)

# Plot variance decomposition
plot_variance_explained(MOFAobject, max_r2=5)

plot_variance_explained(MOFAobject, plot_total = T)[[2]]

# Characterisation of Factor 1
# Association analysis
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("OS_status","molecular_subtype","intnmf_cluster"), 
                                  plot="log_pval"
)

# Plot factor values
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)

# Plot feature weights
#MOFAobject@samples_metadata[["group"]] <- anno_file_rna$intnmf_cluster - did not work
plot_weights(MOFAobject,
             view = "methyl_data",
             factor = 1,
             nfeatures = 10,     
             scale = TRUE)

plot_data_scatter(MOFAobject, 
                  view = "count_data",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "intnmf_cluster"
) + labs(y="RNA expression")

plot_data_heatmap(MOFAobject, 
                  view = "count_data",
                  factor = 1,  
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


# Prediction of clinical subgroup

# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

# Train the model for IGHV
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
model.ighv <- randomForest(IGHV ~ ., data=df[!is.na(df$IGHV),], ntree=10)
df$IGHV <- NULL

# Do predictions
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)


# Train the model for Trisomy12
df$trisomy12 <- as.factor(MOFAobject@samples_metadata$trisomy12)
model.trisomy12 <- randomForest(trisomy12 ~ ., data=df[!is.na(df$trisomy12),], ntree=10)
df$trisomy12 <- NULL

MOFAobject@samples_metadata$trisomy12.pred <- stats::predict(model.trisomy12, df)

MOFAobject@samples_metadata$IGHV.pred_logical <- c("True","Predicted")[as.numeric(is.na(MOFAobject@samples_metadata$IGHV))+1]

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "IGHV.pred",
                  shape_by = "IGHV.pred_logical",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)


# GSEA 
utils::data(reactomeGS)

# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "count_data",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "count_data",
                               sign = "negative"
)

plot_enrichment_heatmap(res.positive)

plot_enrichment_heatmap(res.negative)

plot_enrichment(res.positive, factor = 5, max.pathways = 15)

plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 5, 
  max.pathways = 3
)

# Fit cox model

SurvObject <- Surv(MOFAobject@samples_metadata$TTT, MOFAobject@samples_metadata$treatedAfter)
Z <- get_factors(MOFAobject)[[1]]
fit <- coxph(SurvObject ~ Z) 
fit

# Plot hazard ratio

s <- summary(fit)
coef <- s[["coefficients"]]

df <- data.frame(
  factor = factor(rownames(coef), levels = rev(rownames(coef))),
  p      = coef[,"Pr(>|z|)"], 
  coef   = coef[,"exp(coef)"], 
  lower  = s[["conf.int"]][,"lower .95"], 
  higher = s[["conf.int"]][,"upper .95"]
)

ggplot(df, aes(x=factor, y=coef, ymin=lower, ymax=higher)) +
  geom_pointrange( col='#619CFF') + 
  coord_flip() +
  scale_x_discrete() + 
  labs(y="Hazard Ratio", x="") + 
  geom_hline(aes(yintercept=1), linetype="dotted") +
  theme_bw()

# Plot KM curve

df <- data.frame(
  time = SurvObject[,1], 
  event = SurvObject[,2], Z1 = Z[,1]
)
cut <- surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint
fit <- survfit(Surv(time, event) ~ FactorCluster, df)

ggsurvplot(fit, data = df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100,
           legend = "top", legend.labs = c(paste("low LF 1"), paste("high LF 1")),
           xlab = "Time to treatment", ylab="Survival probability (%)", title= "Factor 1"
)$plot
