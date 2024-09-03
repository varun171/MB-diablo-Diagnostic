# Load necessary libraries
library(caret)
library(pROC)
library(readr)
library(tibble)
library(ggplot2)

# Directories
data_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/data")
result_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/results")
plots_dir <- file.path("~/Documents/GitHub/MB-diablo-Diagnostic/MB-Biomarkers/MB_763_classifier/plots")

# Load the data
microarray_data <- readRDS(file.path(data_dir, "methyl_763data.rds"))
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds"))
methyl_cbtn_data <- as.data.frame(t(methyl_cbtn_data))

# Select columns in the transposed methyl_cbtn_data that are present in the columns of microarray_data
methyl_cbtn_data <- methyl_cbtn_data[, colnames(methyl_cbtn_data) %in% colnames(microarray_data)]
microarray_data <- microarray_data[, colnames(microarray_data) %in% colnames(methyl_cbtn_data)]

# Sample Data
samples_map <- read_tsv(file.path(data_dir, "samples_763map.tsv"))
samples_cbtn_map <- read_tsv(file.path(data_dir, "histologies_medullo.tsv"))

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

# Split the data into training and testing sets
set.seed(42)
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Define a control function for training
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = multiClassSummary)

# Train models using caret
models <- c("rf", "nb", "svmLinear", "gbm", "knn", "xgbTree")
results <- list()

for (model in models) {
  set.seed(42)
  fit <- train(x = X_train, y = y_train, method = model, trControl = train_control)
  results[[model]] <- fit
}

# Evaluate models on the test set
evaluations <- lapply(results, function(fit) {
  y_pred <- predict(fit, X_test)
  y_prob <- predict(fit, X_test, type = "prob")
  cm <- confusionMatrix(y_pred, y_test)
  auc <- multiclass.roc(y_test, y_prob)
  list(accuracy = cm$overall['Accuracy'], confusion_matrix = cm$table, class_report = cm$byClass, auc = auc)
})

# Save the results
saveRDS(evaluations, file.path(result_dir, "caret_model_results.rds"))

# Plot ROC curves
colors <- c("rf" = "blue", "nb" = "red", "svmLinear" = "green", "gbm" = "purple", "knn" = "orange", "xgbTree" = "pink")
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC Curves for Various Models")

for (model in names(evaluations)) {
  roc <- evaluations[[model]]$auc
  plot.roc(roc[[1]], col = colors[[model]], add = TRUE)
}

legend("bottomright", legend = names(colors), col = colors, lwd = 2)
ggsave(file.path(plots_dir, "caret_roc_curves.pdf"))

# Predict subtypes using the trained models
predictions <- lapply(results, function(fit) {
  predict(fit, methyl_cbtn_data)
})

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
