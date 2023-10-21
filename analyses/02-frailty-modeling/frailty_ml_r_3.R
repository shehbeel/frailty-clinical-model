# Building ML Model for Frailty prediction 3
# Author: Shehbeel Arif
# NASA GeneLab Multi-Omics AWG

## LOAD LIBRARIES
# Library data manipulation
library(dplyr)
# Set seed
set.seed(1234)

################################################################################
## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "02-frailty-modeling")


# Set output directories
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Set input directory
data_dir <- file.path(root_dir, "data")

# Declare input file paths
metadata_file <- file.path(root_dir, "data", "GSE144304_meta.txt")
data_file <- file.path(root_dir, "analyses", "01-convert-ensembl-ids", "results", "GSE144304_raw_counts_with_genes.csv")
rfe_file <- file.path(root_dir, "analyses", "03-rfe-feature-selection", "results", "rfecv_selected_features_2.csv")

#######
# Import metadata and data 
metadata <- readr::read_delim(metadata_file)
expression_df <- readr::read_csv(data_file)
rfe_genes <- readr::read_csv(rfe_file)

#######
## PREPROCESS THE DATA
rfe_genes <- rfe_genes[rfe_genes$Rank == 1,]

filtered_expression_df <- expression_df %>%
  distinct(gene_symbol, .keep_all = TRUE) %>% # Remove duplicate genes
  tibble::column_to_rownames("gene_symbol") %>% 
  t() %>%
  as.data.frame() %>%
  dplyr::select(rfe_genes$Gene) %>%
  tibble::rownames_to_column(var = "sample_id")

# # Make the data in the order of the metadata
# expression_df <- expression_df %>%
#   dplyr::select(metadata$sample_id)
# 
# # Check if the expression and metadata are in the same order
# all.equal(colnames(expression_df), metadata$sample_id)


# Merge meta and gene matrix data
df <- merge(metadata, filtered_expression_df, by = "sample_id")
df <- df %>%
  dplyr::select(-c("sample_id", "condition"))
##df$gender <- unclass(df$gender)
df$gender<-c(male=1,female=0)[df$gender]
#df$frailty<-c(1=frail,female=normal)[df$frailty]

# Move condition column to as first column
# df <- df %>%
#   select(frailty, everything())

# Convert to matrix
#df <- df %>% as.matrix()
#df <- as.matrix(sapply(df, as.numeric))  

# look at the class of each column
sapply(df, class)

#df <- df %>% mutate_at(2:20299, as.numeric)
#df[,2:20]

################################################################################
# MRMR Feature Selection
library(mRMRe)
# Ignore gender column
f_data <- mRMR.data(data = df[,2:165])
# Run MRMR
results <- mRMR.classic(data = f_data, target_indices = 1,
                        feature_count = 10)
# results <- mRMR.ensemble(data = f_data, target_indices = 1,
#                          solution_count = 1, feature_count = 10)
solutions(results)
causality(results)
df[,c(10,37,113,116,125,134,138,149,153,164)]

#################
require(nnet)
# feature_count = 5
multinom_fit <- multinom(frailty ~ `HABP2` + `MYBPHL` + `TJP3` + `TEX56P` + `CCL14`, data = df, model = TRUE)
# feature_count = 10
multinom_fit <- multinom(frailty ~ `EFNA5` + `HABP2` + `MYBPHL` + `TJP3` + `TEX56P` + `ACRP1` + `LRTM1` + `CCL14` + `SEMA3E` + `PDF`, data = df, model = TRUE)
# Final model
multinom_fit <- multinom(frailty ~ `HABP2` + `TJP3` + `TEX56P` + `CCL14` + `PDF`, data = df, model = TRUE)

summary(multinom_fit)
round((1 - pnorm(abs(summary(multinom_fit)$coefficients/summary(multinom_fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(multinom_fit))

head(pp <- fitted(multinom_fit))
library(sjPlot)
tab_model(multinom_fit, digits = 4)

prediction = predict(multinom_fit, newdata = df, "class")
sum(prediction == df$frailty)/nrow(df)

# Load ROC-AUC library
library(pROC)
# create roc curve
roc_object <- roc(df$frailty, as.numeric(prediction), plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, xlab="False Positive (%)",
                  ylab="True Positive (%)", col="#377eb8", lwd=4, print.auc=TRUE,
                  auc.polygon=TRUE)

# calculate area under curve
auc(roc_object)

# Goodness-of-fit Test
# https://peopleanalytics-regression-book.org/multinomial-logistic-regression-for-nominal-category-outcomes.html
DescTools::PseudoR2(multinom_fit, 
                    which = c("McFadden", "CoxSnell", "Nagelkerke"))
# The most approachable method to assess model confidence is the Hosmer-Lemeshow test

################################################################################
# Gene Boxplots
library(ggpubr)
library(ggplot2)

# HABP2 Expression
HABP2_bp <- ggboxplot(df, 
                        # Specify x values
                        x = "frailty",
                        # Specify y values
                        y = "HABP2",
                        # Color in the box plot
                        fill = "frailty",
                        # Specify color palette
                        palette = "npg", 
                        # Add x-axis label
                        xlab="Frailty",
                        # Add y-axis label
                        ylab="HABP2 Gene Expression",
                        # Add points
                        add = "jitter") +
  # Add p-value
  stat_compare_means(method = "t.test")

# TJP3
TJP3_bp <- ggboxplot(df, 
                      # Specify x values
                      x = "frailty",
                      # Specify y values
                      y = "TJP3",
                      # Color in the box plot
                      fill = "frailty",
                      # Specify color palette
                      palette = "npg", 
                      # Add x-axis label
                      xlab="Frailty",
                      # Add y-axis label
                      ylab="TJP3 Gene Expression",
                      # Add points
                      add = "jitter") +
  # Add p-value
  stat_compare_means(method = "t.test")

# TEX56P
TEX56P_bp <- ggboxplot(df, 
                     # Specify x values
                     x = "frailty",
                     # Specify y values
                     y = "TEX56P",
                     # Color in the box plot
                     fill = "frailty",
                     # Specify color palette
                     palette = "npg", 
                     # Add x-axis label
                     xlab="Frailty",
                     # Add y-axis label
                     ylab="TEX56P Gene Expression",
                     # Add points
                     add = "jitter") +
  # Add p-value
  stat_compare_means(method = "t.test")

# CCL14
CCL14_bp <- ggboxplot(df, 
                     # Specify x values
                     x = "frailty",
                     # Specify y values
                     y = "CCL14",
                     # Color in the box plot
                     fill = "frailty",
                     # Specify color palette
                     palette = "npg", 
                     # Add x-axis label
                     xlab="Frailty",
                     # Add y-axis label
                     ylab="CCL14 Gene Expression",
                     # Add points
                     add = "jitter") +
  # Add p-value
  stat_compare_means(method = "t.test")

# PDF
PDF_bp <- ggboxplot(df, 
                     # Specify x values
                     x = "frailty",
                     # Specify y values
                     y = "PDF",
                     # Color in the box plot
                     fill = "frailty",
                     # Specify color palette
                     palette = "npg", 
                     # Add x-axis label
                     xlab="Frailty",
                     # Add y-axis label
                     ylab="PDF Gene Expression",
                     # Add points
                     add = "jitter") +
  # Add p-value
  stat_compare_means(method = "t.test")

# Combine into one plot
figure <- ggarrange(HABP2_bp, TJP3_bp, TEX56P_bp,
                    CCL14_bp, PDF_bp,
                    labels = c("A", "B", "C", "D", "E"),
                    ncol = 3, nrow = 2)

# View combined boxplot figure
print(figure)

# Save plot
boxplot_figure_png <- file.path(plots_dir, "boxplot_figure_panel.png")
# Export plot
ggsave(boxplot_figure_png, figure, width = 12, height = 8)





