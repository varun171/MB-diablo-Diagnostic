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
library(shapr)
library(caret)

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

# Split the data into training and testing sets
set.seed(42)
trainIndex <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Train the Random Forest model
rf_model <- randomForest(x = X_train, y = y_train, ntree = 100)

# Train the Naive Bayes model
nb_model <- naiveBayes(x = X_train, y = y_train)

# Train the SVM model
svm_model <- svm(x = X_train, y = y_train, kernel = "linear", probability = TRUE)

# Train the GBM model
gbm_model <- gbm.fit(x = as.matrix(X_train), y = as.numeric(y_train) - 1, distribution = "multinomial", n.trees = 100, interaction.depth = 3)

# Train the k-NN model and evaluate on the test set
k <- 5  # Choose the number of neighbors

# Prepare the data for XGBoost
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = as.numeric(y_train) - 1)
dtest <- xgb.DMatrix(data = as.matrix(X_test), label = as.numeric(y_test) - 1)

# Train the XGBoost model
param <- list(objective = "multi:softprob", num_class = length(levels(y_train)), eval_metric = "mlogloss")
xgb_model <- xgb.train(params = param, data = dtrain, nrounds = 100)

# Evaluate the models on the test set for validation
# Random Forest
y_pred_rf_test <- predict(rf_model, X_test)
rf_probs <- predict(rf_model, X_test, type = "prob")
accuracy_rf <- confusionMatrix(y_pred_rf_test, y_test)$overall['Accuracy']
conf_matrix_rf <- confusionMatrix(y_pred_rf_test, y_test)$table
class_report_rf <- confusionMatrix(y_pred_rf_test, y_test)$byClass

# Naive Bayes
y_pred_nb_test <- predict(nb_model, X_test)
nb_probs <- predict(nb_model, X_test, type = "raw")
accuracy_nb <- confusionMatrix(y_pred_nb_test, y_test)$overall['Accuracy']
conf_matrix_nb <- confusionMatrix(y_pred_nb_test, y_test)$table
class_report_nb <- confusionMatrix(y_pred_nb_test, y_test)$byClass

# SVM
y_pred_svm_test <- predict(svm_model, X_test)
svm_probs <- attr(predict(svm_model, X_test, probability = TRUE), "probabilities")
accuracy_svm <- confusionMatrix(y_pred_svm_test, y_test)$overall['Accuracy']
conf_matrix_svm <- confusionMatrix(y_pred_svm_test, y_test)$table
class_report_svm <- confusionMatrix(y_pred_svm_test, y_test)$byClass

# GBM
y_pred_gbm_test <- predict(gbm_model, as.matrix(X_test), n.trees = 100, type = "response")
y_pred_gbm_prob <- as.data.frame(y_pred_gbm_test)
accuracy_gbm <- confusionMatrix(factor(apply(y_pred_gbm_prob, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train)), y_test)$overall['Accuracy']
conf_matrix_gbm <- confusionMatrix(factor(apply(y_pred_gbm_prob, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train)), y_test)$table
class_report_gbm <- confusionMatrix(factor(apply(y_pred_gbm_prob, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train)), y_test)$byClass

# k-NN
y_pred_knn_test <- knn(train = X_train, test = X_test, cl = y_train, k = k)
accuracy_knn <- confusionMatrix(y_pred_knn_test, y_test)$overall['Accuracy']
conf_matrix_knn <- confusionMatrix(y_pred_knn_test, y_test)$table
class_report_knn <- confusionMatrix(y_pred_knn_test, y_test)$byClass

# XGBoost
y_pred_xgb_test <- predict(xgb_model, dtest)
y_pred_xgb_prob <- matrix(y_pred_xgb_test, ncol = length(levels(y_train)), byrow = TRUE)
accuracy_xgb <- confusionMatrix(factor(apply(y_pred_xgb_prob, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train)), y_test)$overall['Accuracy']
conf_matrix_xgb <- confusionMatrix(factor(apply(y_pred_xgb_prob, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train)), y_test)$table
class_report_xgb <- confusionMatrix(factor(apply(y_pred_xgb_prob, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train)), y_test)$byClass

# Save the results
results <- list(
  Random_Forest = list(accuracy = accuracy_rf, confusion_matrix = conf_matrix_rf, class_report = class_report_rf),
  Naive_Bayes = list(accuracy = accuracy_nb, confusion_matrix = conf_matrix_nb, class_report = class_report_nb),
  SVM = list(accuracy = accuracy_svm, confusion_matrix = conf_matrix_svm, class_report = class_report_svm),
  GBM = list(accuracy = accuracy_gbm, confusion_matrix = conf_matrix_gbm, class_report = class_report_gbm),
  kNN = list(accuracy = accuracy_knn, confusion_matrix = conf_matrix_knn, class_report = class_report_knn),
  XGBoost = list(accuracy = accuracy_xgb, confusion_matrix = conf_matrix_xgb, class_report = class_report_xgb)
)

saveRDS(results, file.path(result_dir, "model_results.rds"))

# Models Performance
# Calculate ROC AUC for Random Forest
roc_rf <- multiclass.roc(y_test, rf_probs)
auc_rf <- auc(roc_rf)

# Calculate ROC AUC for Naive Bayes
roc_nb <- multiclass.roc(y_test, nb_probs)
auc_nb <- auc(roc_nb)

# Calculate ROC AUC for SVM
roc_svm <- multiclass.roc(y_test, svm_probs)
auc_svm <- auc(roc_svm)

# Plot ROC curves # not working need to fix
plot.roc(roc_rf[[1]], col = "blue", main = "ROC Curves for Various Models")
plot.roc(roc_nb[[1]], col = "red", add = TRUE)
plot.roc(roc_svm[[1]], col = "green", add = TRUE)
legend("bottomright", legend = c("Random Forest", "Naive Bayes", "SVM"),
       col = c("blue", "red", "green"), lwd = 2)

# Save the ROC plot
ggsave(file.path(plots_dir, "roc_curves.pdf"))

# Ensure the methyl_cbtn_data has the same columns as the training data
methyl_cbtn_data <- methyl_cbtn_data[, colnames(methyl_cbtn_data) %in% colnames(X_train)]

# Predict subtypes using the trained models
rf_predictions <- predict(rf_model, methyl_cbtn_data)
nb_predictions <- predict(nb_model, methyl_cbtn_data)
svm_predictions <- predict(svm_model, methyl_cbtn_data)
gbm_predictions <- predict(gbm_model, as.matrix(methyl_cbtn_data), n.trees = 100, type = "response")
gbm_predictions <- factor(apply(gbm_predictions, 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train))
knn_predictions <- knn(train = X_train, test = methyl_cbtn_data, cl = y_train, k = k)
#xgb_predictions <- predict(xgb_model, xgb.DMatrix(data = as.matrix(methyl_cbtn_data)))
#xgb_predictions <- factor(apply(matrix(xgb_predictions, ncol = length(levels(y_train)), byrow = TRUE), 1, which.max), levels = 1:length(levels(y_train)), labels = levels(y_train))

# Assuming methyl_cbtn_data row names are sample IDs and match with sample_cbtn_map
sample_cbtn_map <- read_tsv(file.path(data_dir, "histologies_medullo.tsv"))

# Ensure the order of samples matches between predictions and actual subtypes
sample_cbtn_map <- sample_cbtn_map[match(rownames(methyl_cbtn_data), sample_cbtn_map$Kids_First_Biospecimen_ID), ]

# Add predicted subtypes to the sample_cbtn_map # xgb not working
sample_cbtn_map$rf_predicted_subtype <- rf_predictions
sample_cbtn_map$nb_predicted_subtype <- nb_predictions
sample_cbtn_map$svm_predicted_subtype <- svm_predictions
sample_cbtn_map$gbm_predicted_subtype <- gbm_predictions
sample_cbtn_map$knn_predicted_subtype <- knn_predictions
#sample_cbtn_map$xgb_predicted_subtype <- xgb_predictions

# Create balloon plots comparing predicted subtypes and actual molecular subtypes # xgb not working
create_balloon_plot <- function(data, model_name) {
  balloon_data <- as.data.frame(table(data$molecular_subtype, data[[paste0(model_name, "_predicted_subtype")]]))
  p <- ggplot(balloon_data, aes(Var1, Var2, size = Freq)) +
    geom_point(alpha = 0.7) +
    labs(title = paste("Balloon Plot:", model_name, "Predicted vs Actual Subtypes"),
         x = "Actual Molecular Subtype",
         y = paste("Predicted Subtype (", model_name, ")", sep = "")) +
    theme_minimal()
  ggsave(file.path(plots_dir, paste0("balloon_plot_", model_name, ".pdf")), plot = p)
}

create_balloon_plot(sample_cbtn_map, "rf")
create_balloon_plot(sample_cbtn_map, "nb")
create_balloon_plot(sample_cbtn_map, "svm")
create_balloon_plot(sample_cbtn_map, "gbm")
create_balloon_plot(sample_cbtn_map, "knn")
create_balloon_plot(sample_cbtn_map, "xgb")


########################################################

# Create balloon plots comparing predicted subtypes and actual molecular subtypes
create_balloon_plot <- function(data, model_name) {
  balloon_data <- as.data.frame(table(data$molecular_subtype, data[[paste0(model_name, "_predicted_subtype")]]))
  p <- ggplot(balloon_data, aes(Var1, Var2, size = Freq)) +
    geom_point(alpha = 0.7) +
    labs(title = paste("Balloon Plot:", model_name, "Predicted vs Actual Subtypes"),
         x = "Actual Molecular Subtype",
         y = paste("Predicted Subtype (", model_name, ")", sep = "")) +
    theme_minimal()
  ggsave(file.path(plots_dir, paste0("balloon_plot_", model_name, ".pdf")), plot = p)
}

create_balloon_plot(sample_cbtn_map, "rf")
create_balloon_plot(sample_cbtn_map, "nb")
create_balloon_plot(sample_cbtn_map, "svm")
create_balloon_plot(sample_cbtn_map, "gbm")
create_balloon_plot(sample_cbtn_map, "knn")


# Select the specified columns
selected_columns <- c("sample_id", "molecular_subtype", "dkfz_v12_methylation_subclass",
                      "rf_predicted_subtype", "nb_predicted_subtype",
                      "svm_predicted_subtype", "gbm_predicted_subtype", "knn_predicted_subtype")

# Create a new data frame with only the selected columns
sample_cbtn_map_selected <- sample_cbtn_map[selected_columns]

# Save the data frame with only the selected columns as a TSV file
write.table(sample_cbtn_map_selected, file.path(result_dir, "sample_cbtn_map_selected.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Save the full data frame as a TSV file
write.table(sample_cbtn_map, file.path(result_dir, "sample_cbtn_map.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
######################################################

intnmf_clusters <- read_tsv(file.path(result_dir, "intnmf_clusters_Komal.tsv"))

# Assuming intnmf_clusters is another DataFrame with sample_id and intnmf_cluster_id
# Merge the selected columns with the intnmf_clusters data frame using sample_id as the common field
sample_cbtn_map_selected_intnmf <- merge(sample_cbtn_map_selected, intnmf_clusters, by = "sample_id")

# Save the data frame with only the selected columns and intnmf_cluster_id as a TSV file
write.table(sample_cbtn_map_selected_intnmf, file.path(result_dir, "sample_cbtn_map_selected_with_intnmf_cluster.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Create a contingency table for mm_cluster and rf_predicted_subtype
contingency_table <- table(sample_cbtn_map_selected_intnmf$mm_cluster, sample_cbtn_map_selected_intnmf$rf_predicted_subtype)

# Convert the table to a data frame for ggplot2
contingency_df <- as.data.frame(contingency_table)
names(contingency_df) <- c("mm_cluster", "rf_predicted_subtype", "count")

# Create the balloon plot
balloon_plot <- ggplot(contingency_df, aes(x = rf_predicted_subtype, y = mm_cluster)) +
  geom_point(aes(size = count), shape = 21, fill = "blue") +
  scale_size(range = c(3, 15)) +
  theme_minimal() +
  labs(title = "Balloon Plot between mm_cluster and rf_predicted_subtype",
       x = "rf_predicted_subtype",
       y = "mm_cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the balloon plot
plot_path <- file.path(plots_dir, "balloon_plot_mm_cluster_rf_predicted_subtype.pdf")
ggsave(plot_path, plot = balloon_plot, width = 12, height = 12)

# Print message indicating the plot has been saved
print(paste("Balloon plot saved at:", plot_path))

############################################################################################
# SHAP analysis for feature identifications - not for random forest

df1 <- shapr:::get_supported_models()

# Use the shapr package to calculate SHAP values
explainer <- shapr(data, rf_model)

# Define the sample data to explain
explainer_data <- data[1:5, , drop = FALSE]  # You can change this to the data you want to explain

# Calculate the SHAP values
shap_values <- explain(explainer_data, explainer, approach = "empirical", prediction_zero = mean(target))

# Plot the SHAP values for the first instance
plot_shap_values(shap_values$shap_values[1, ], explainer_data[1, ])

# Save SHAP values
saveRDS(shap_values, file.path(result_dir, "shap_values_rf.rds"))

# Print SHAP values summary for all instances
print(shap_values$shap_values)

######################################################################

# Load the methyl cbtn data
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds")) 

# Extract one column with all rows (assuming the first column for this example)
single_column <- methyl_cbtn_data[, 1:10, drop = FALSE]

# Save the single column as an RDS file
saveRDS(single_column, file.path(data_dir, "ten_sample_column.rds"))

# Output the column name for reference
cat("Column Name:", colnames(single_column), "\n")

########################################################################

# Load the methyl cbtn data
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds")) 

# Extract one column with all rows (assuming the first ten columns for this example)
single_column <- methyl_cbtn_data[, 1:2, drop = FALSE]

# Save the single column as a CSV file, including row names
write.csv(single_column, file.path(data_dir, "ten_sample_column.csv"), row.names = TRUE)

# Output the column names for reference
cat("Column Names:", colnames(single_column), "\n")

#################################################################
# Load the methyl cbtn data
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds")) 

# Extract one column with all rows (assuming the first ten columns for this example)
single_column <- methyl_cbtn_data[, 1:10, drop = FALSE]

# Save the single column as a text file, including row names
write.table(single_column, file.path(data_dir, "ten_sample_column.txt"), row.names = TRUE, col.names = TRUE, sep = "\t")

# Output the column names for reference
cat("Column Names:", colnames(single_column), "\n")


