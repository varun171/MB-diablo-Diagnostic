# Load Required Libraries
library(dplyr)
library(readr)
library(tidyr)

# Read the expression matrix and sample descriptions
expr_matrix <- read_tsv("GSE124814_HW_expr_matrix.tsv")
sample_descriptions <- read_csv("GSE124814_sample_descriptions.csv")

# Update the characteristics: subgroup relabeled column for normal samples
sample_descriptions <- sample_descriptions %>%
  mutate(`characteristics: subgroup relabeled` = ifelse(`source name` == "Normal", "Normal", `characteristics: subgroup relabeled`))

# Filter the expression matrix to include only MYC and FZD1 genes
genes_of_interest <- c("MYC", "FZD1")
expr_matrix_filtered <- expr_matrix %>%
  filter(Gene_Symbol %in% genes_of_interest)

# Reshape the expression matrix to have samples as columns and genes as rows
expr_matrix_long <- expr_matrix_filtered %>%
  gather(key = "SampleID", value = "Expression", -Gene_Symbol)

# Merge the expression data with the sample descriptions
expr_with_descriptions <- expr_matrix_long %>%
  inner_join(sample_descriptions, by = c("SampleID" = "Sample name"))

# Calculate the mean expression for MYC and FZD1 in normal and G3 MB samples
mean_expression <- expr_with_descriptions %>%
  group_by(Gene_Symbol, `characteristics: subgroup relabeled`) %>%
  summarise(MeanExpression = mean(Expression, na.rm = TRUE))

# Take the absolute values of the mean expression
mean_expression <- mean_expression %>%
  mutate(MeanExpression = abs(MeanExpression))

# Calculate the absolute fold change between normal and G3 MB samples
fold_change_absolute <- mean_expression %>%
  spread(key = `characteristics: subgroup relabeled`, value = MeanExpression) %>%
  mutate(FoldChange_G3_vs_Normal = G3 / Normal)

# Calculate the log2 fold change between normal and G3 MB samples
fold_change_log2 <- mean_expression %>%
  spread(key = `characteristics: subgroup relabeled`, value = MeanExpression) %>%
  mutate(Log2_FoldChange_G3_vs_Normal = log2(G3 / Normal))

# Perform t-test for each gene
t_test_results <- expr_with_descriptions %>%
  group_by(Gene_Symbol) %>%
  summarise(p_value = t.test(Expression[`characteristics: subgroup relabeled` == "G3"],
                             Expression[`characteristics: subgroup relabeled` == "Normal"])$p.value)

# Merge fold change results with p-values
fold_change_results <- fold_change_log2 %>%
  left_join(t_test_results, by = "Gene_Symbol")

# Print the results
print("Fold Change and P-values:")
print(fold_change_results)

