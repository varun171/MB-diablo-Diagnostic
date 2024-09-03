# Load necessary libraries
library(dplyr)
library(randomForest)
library(e1071)
library(caret)
library(pROC)
library(preprocessCore)
library(readr)
library(tibble)
library(sva)
library(limma)
library(gbm)
library(class)
library(xgboost)
library(ggplot2)
library(doParallel)

# Directories
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
result_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/plots")

# Load the methyl data (cavalli etal)
microarray_data <- readRDS(file.path(data_dir, "methyl_763data.rds")) # all features

# Read methyl cbtn data
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds")) 
methyl_cbtn_data <- as.data.frame(t(methyl_cbtn_data))

# Select columns in the transposed methyl_cbtn_data that are present in the columns of microarray_data
methyl_cbtn_data <- methyl_cbtn_data[, colnames(methyl_cbtn_data) %in% colnames(microarray_data)]

# Select columns in microarray_data that are present in the filtered_transposed_methyl_cbtn_data
microarray_data <- microarray_data[, colnames(microarray_data) %in% colnames(methyl_cbtn_data)]

# Sample Data
samples_map <- read_tsv(file.path(data_dir, "samples_763map.tsv"))

# Select relevant columns from samples_map
samples_map <- samples_map %>% dplyr::select(sample_id, subtype)

# Merge the microarray data with the samples_map to align sample ids with subtypes
merged_data <- microarray_data %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  inner_join(samples_map, by = "sample_id")

# Prepare the features (X) and target (y)
X <- merged_data %>% dplyr::select(-sample_id, -subtype)
y <- merged_data$subtype

# Encode the target variable
y <- as.factor(y)

# Ensure factor levels are valid R variable names
levels(y) <- make.names(levels(y))

# Split the data into training and testing sets
set.seed(142)
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Define a control function for training
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary)

# Set up parallel processing
cl <- makeCluster(8) # Use 2 cores for parallel processing
registerDoParallel(cl)

# Define hyperparameter grids for Random Forest and Naive Bayes
rf_grid <- expand.grid(
  mtry = c(50, 100, 150)
)
rf_ntree <- c(100, 200, 300)
rf_nodesize <- c(1, 5, 10)
rf_maxnodes <- c(50, 100, 200)

nb_grid <- expand.grid(
  fL = c(0, 0.5, 1),
  usekernel = c(TRUE, FALSE),
  adjust = c(0.5, 1, 1.5)
)

# Train Random Forest with additional parameters
rf_results <- list()
for (mtry_val in rf_grid$mtry) {
  for (ntree_val in rf_ntree) {
    for (nodesize_val in rf_nodesize) {
      for (maxnodes_val in rf_maxnodes) {
        set.seed(42)
        rf_fit <- train(
          x = X_train, y = y_train, method = "rf",
          trControl = train_control,
          tuneGrid = data.frame(mtry = mtry_val),
          ntree = ntree_val,
          nodesize = nodesize_val,
          maxnodes = maxnodes_val
        )
        rf_results[[paste(mtry_val, ntree_val, nodesize_val, maxnodes_val, sep = "_")]] <- rf_fit
      }
    }
  }
}

# Train Naive Bayes with additional parameters
nb_results <- list()
for (i in 1:nrow(nb_grid)) {
  set.seed(42)
  nb_fit <- train(
    x = X_train, y = y_train, method = "nb",
    trControl = train_control,
    tuneGrid = nb_grid[i, , drop = FALSE]
  )
  nb_results[[paste(nb_grid[i,], collapse = "_")]] <- nb_fit
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Evaluate models on the test set
rf_evaluations <- lapply(rf_results, function(fit) {
  y_pred <- predict(fit, X_test)
  y_prob <- predict(fit, X_test, type = "prob")
  cm <- confusionMatrix(y_pred, y_test)
  auc <- multiclass.roc(y_test, as.matrix(y_prob))
  list(accuracy = cm$overall['Accuracy'], confusion_matrix = cm$table, class_report = cm$byClass, auc = auc)
})

nb_evaluations <- lapply(nb_results, function(fit) {
  y_pred <- predict(fit, X_test)
  y_prob <- predict(fit, X_test, type = "prob")
  cm <- confusionMatrix(y_pred, y_test)
  auc <- multiclass.roc(y_test, as.matrix(y_prob))
  list(accuracy = cm$overall['Accuracy'], confusion_matrix = cm$table, class_report = cm$byClass, auc = auc)
})

# Save the results
saveRDS(rf_evaluations, file.path(result_dir, "rf_model_results_with_tuning.rds"))
saveRDS(nb_evaluations, file.path(result_dir, "nb_model_results_with_tuning.rds"))

# Plot ROC curves for the best models
best_rf <- rf_results[[which.max(sapply(rf_evaluations, function(x) x$accuracy))]]
best_nb <- nb_results[[which.max(sapply(nb_evaluations, function(x) x$accuracy))]]

# Calculate ROC AUC for the best Random Forest model
roc_rf <- multiclass.roc(y_test, predict(best_rf, X_test, type = "prob"))
auc_rf <- auc(roc_rf)

# Calculate ROC AUC for the best Naive Bayes model
roc_nb <- multiclass.roc(y_test, predict(best_nb, X_test, type = "prob"))
auc_nb <- auc(roc_nb)

# Plot ROC curves
plot.roc(roc_rf[[1]], col = "blue", main = "ROC Curves for Best Models")
plot.roc(roc_nb[[1]], col = "red", add = TRUE)
legend("bottomright", legend = c("Random Forest", "Naive Bayes"),
       col = c("blue", "red"), lwd = 2)
ggsave(file.path(plots_dir, "best_model_roc_curves_with_tuning.pdf"))

# Predict subtypes using the best models
rf_predictions <- predict(best_rf, methyl_cbtn_data)
nb_predictions <- predict(best_nb, methyl_cbtn_data)

# Assuming methyl_cbtn_data row names are sample IDs and match with sample_cbtn_map
sample_cbtn_map <- read_tsv(file.path(data_dir, "histologies_medullo.tsv"))

# Ensure the order of samples matches between predictions and actual subtypes
sample_cbtn_map <- sample_cbtn_map[match(rownames(methyl_cbtn_data), sample_cbtn_map$Kids_First_Biospecimen_ID), ]

# Add predicted subtypes to the sample_cbtn_map
sample_cbtn_map$rf_predicted_subtype <- rf_predictions
sample_cbtn_map$nb_predicted_subtype <- nb_predictions

# Create balloon plots comparing predicted subtypes and actual molecular subtypes for Random Forest and Naive Bayes
create_balloon_plot <- function(data, model_name) {
  balloon_data <- as.data.frame(table(data$molecular_subtype, data[[paste0(model_name, "_predicted_subtype")]]))
  p <- ggplot(balloon_data, aes(Var1, Var2, size = Freq)) +
    geom_point(alpha = 0.7) +
    labs(title = paste("Balloon Plot:", model_name, "Predicted vs Actual Subtypes"),
         x = "Actual Molecular Subtype",
         y = paste("Predicted Subtype (", model_name, ")", sep = "")) +
    theme_minimal()
  ggsave(file.path(plots_dir, paste0("tune_balloon_plot_", model_name, ".pdf")), plot = p)
}

create_balloon_plot(sample_cbtn_map, "rf")
create_balloon_plot(sample_cbtn_map, "nb")
