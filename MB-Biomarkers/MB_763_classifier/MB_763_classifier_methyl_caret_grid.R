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
library(klaR)

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
set.seed(42)
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Define a control function for training
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary)

# Define hyperparameter grids for each model
tune_grids <- list(
  rf = expand.grid(mtry = c(50, 100)),
  nb = expand.grid(fL = c(0, 0.5, 1), usekernel = c(TRUE, FALSE), adjust = c(1)),
  svmLinear = expand.grid(C = c(0.1, 1)),
  gbm = expand.grid(interaction.depth = c(1, 3), n.trees = c(50, 100), shrinkage = c(0.01, 0.1), n.minobsinnode = c(10)),
  knn = expand.grid(k = c(3, 5)),
  xgbTree = expand.grid(nrounds = c(50, 100), max_depth = c(3, 6), eta = c(0.01, 0.1), gamma = c(0, 0.1), colsample_bytree = c(0.5), min_child_weight = c(1))
)

# Set up parallel processing
cl <- makeCluster(4) # Use 2 cores for parallel processing
registerDoParallel(cl)

# Train models using caret with hyperparameter tuning
results <- list()

for (model in names(tune_grids)) {
  set.seed(42)
  fit <- train(x = X_train, y = y_train, method = model, trControl = train_control, tuneGrid = tune_grids[[model]])
  results[[model]] <- fit
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Evaluate models on the test set
evaluations <- lapply(results, function(fit) {
  y_pred <- predict(fit, X_test)
  y_prob <- predict(fit, X_test, type = "prob")
  cm <- confusionMatrix(y_pred, y_test)
  auc <- multiclass.roc(y_test, as.matrix(y_prob))
  list(accuracy = cm$overall['Accuracy'], confusion_matrix = cm$table, class_report = cm$byClass, auc = auc)
})

# Save the results
saveRDS(evaluations, file.path(result_dir, "caret_model_results_with_tuning.rds"))

# Plot ROC curves
colors <- c("rf" = "blue", "nb" = "red", "svmLinear" = "green", "gbm" = "purple", "knn" = "orange", "xgbTree" = "pink")
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC Curves for Various Models")

for (model in names(evaluations)) {
  roc <- evaluations[[model]]$auc
  plot.roc(roc[[1]], col = colors[[model]], add = TRUE)
}

legend("bottomright", legend = names(colors), col = colors, lwd = 2)
ggsave(file.path(plots_dir, "caret_roc_curves_with_tuning.pdf"))

# Predict subtypes using the trained models
predictions <- lapply(results, function(fit) {
  predict(fit, methyl_cbtn_data)
})

# Revert factor levels to original names for evaluation
levels(y_train) <- make.names(levels(y_train), unique = TRUE)
levels(y_test) <- make.names(levels(y_test), unique = TRUE)

# Add predicted subtypes to the sample_cbtn_map
for (model in names(predictions)) {
  sample_cbtn_map[[paste0(model, "_predicted_subtype")]] <- predictions[[model]]
}

# Create balloon plots comparing predicted subtypes and actual molecular subtypes
create_balloon_plot <- function(data, model_name) {
  balloon_data <- as.data.frame(table(data$molecular_subtype, data[[paste0(model_name, "_predicted_subtype")]]))
  p <- ggplot(balloon_data, aes(Var1, Var2, size = Freq)) +
    geom_point(alpha = 0.7) +
    labs(title = paste("Balloon Plot:", model_name, "Predicted vs Actual Subtypes"),
         x = "Actual Molecular Subtype",
         y = paste("Predicted Subtype (", model_name, ")", sep = "")) +
    theme_minimal()
  ggsave(file.path(plots_dir, paste0("caret_balloon_plot_", model_name, ".pdf")), plot = p)
}

for (model in names(predictions)) {
  create_balloon_plot(sample_cbtn_map, model)
}
