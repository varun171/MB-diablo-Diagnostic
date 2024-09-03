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

# Load the microarray data
#microarray_data <- read_tsv(file.path(data_dir, "norm_763counts.tsv")) %>% column_to_rownames()
microarray_data <- read_tsv(file.path(data_dir, "methyl_763data.tsv")) %>% column_to_rownames() # 1000 features
microarray_data <- readRDS(file.path(data_dir, "methyl_763data.rds")) # all features

# read methyl cbtn data
methyl_cbtn_data <- readRDS(file.path(data_dir, "methyl-beta-values_medullo.rds")) 
transposed_methyl_cbtn_data <- as.data.frame(t(methyl_cbtn_data))

# Select columns in the transposed methyl_cbtn_data that are present in the columns of microarray_data
filtered_transposed_methyl_cbtn_data <- transposed_methyl_cbtn_data[, colnames(transposed_methyl_cbtn_data) %in% colnames(microarray_data)]

# Step 3: Select columns in microarray_data that are present in the filtered_transposed_methyl_cbtn_data
filtered_microarray_data <- microarray_data[, colnames(microarray_data) %in% colnames(filtered_transposed_methyl_cbtn_data)]
microarray_data <- filtered_microarray_data
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

# Quantile Normalization for Microarray Data
X_normalized <- normalize.quantiles(as.matrix(X))
X <- as.data.frame(X_normalized)
colnames(X) <- colnames(merged_data %>% dplyr::select(-sample_id, -subtype))
rownames(X) <- rownames(merged_data)

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

# Load RNA-seq data
rna_seq_data <- read_tsv(file.path(data_dir, "rna_seq_data.tsv")) %>% column_to_rownames()

# Apply Combat-Seq for RNA-seq data normalization
# Here, assuming rna_seq_data is count data
rna_seq_data_normalized <- ComBat_seq(as.matrix(rna_seq_data), batch = rep(1, nrow(rna_seq_data)))

# Perform batch correction using Combat
combined_data <- cbind(X_train, rna_seq_data_normalized)
batch <- c(rep("microarray", nrow(X_train)), rep("RNAseq", nrow(rna_seq_data_normalized)))
combat_data <- ComBat(dat = as.matrix(combined_data), batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)

# Split the combined Combat data back
X_train_combat <- combat_data[, 1:ncol(X_train)]
rna_seq_data_combat <- combat_data[, (ncol(X_train) + 1):ncol(combat_data)]

# Predict subtypes on RNA-seq data using the trained Random Forest model
y_pred_rf <- predict(rf_model, rna_seq_data_combat)

# Predict subtypes on RNA-seq data using the trained Naive Bayes model
y_pred_nb <- predict(nb_model, rna_seq_data_combat)

# Evaluate the model on the test set for validation
# Random Forest
y_pred_rf_test <- predict(rf_model, X_test)
accuracy_rf <- confusionMatrix(y_pred_rf_test, y_test)$overall['Accuracy']
conf_matrix_rf <- confusionMatrix(y_pred_rf_test, y_test)$table
class_report_rf <- confusionMatrix(y_pred_rf_test, y_test)$byClass

print("Random Forest Classifier")
print(paste("Accuracy: ", accuracy_rf))
print("Confusion Matrix:")
print(conf_matrix_rf)
print("Classification Report:")
print(class_report_rf)

# Naive Bayes
y_pred_nb_test <- predict(nb_model, X_test)
accuracy_nb <- confusionMatrix(y_pred_nb_test, y_test)$overall['Accuracy']
conf_matrix_nb <- confusionMatrix(y_pred_nb_test, y_test)$table
class_report_nb <- confusionMatrix(y_pred_nb_test, y_test)$byClass

print("Naive Bayes Classifier")
print(paste("Accuracy: ", accuracy_nb))
print("Confusion Matrix:")
print(conf_matrix_nb)
print("Classification Report:")
print(class_report_nb)

# Models Performance
# Calculate ROC AUC for Random Forest
rf_probs <- predict(rf_model, X_test, type = "prob")
roc_rf <- multiclass.roc(y_test, rf_probs)
auc_rf <- auc(roc_rf)

# Calculate ROC AUC for Naive Bayes
nb_probs <- predict(nb_model, X_test, type = "raw")
roc_nb <- multiclass.roc(y_test, nb_probs)
auc_nb <- auc(roc_nb)

# Plot ROC curves
plot.roc(roc_rf[[1]], col = "blue", main = "ROC Curves for Random Forest and Naive Bayes")
plot.roc(roc_nb[[1]], col = "red", add = TRUE)
legend("bottomright", legend = c("Random Forest", "Naive Bayes"), col = c("blue", "red"), lwd = 2)

print(paste("Random Forest AUC: ", auc_rf))
print(paste("Naive Bayes AUC: ", auc_nb))
