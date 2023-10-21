# Building ML Model for Frailty prediction 2
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
  dplyr::select(-c("sample_id", "condition", "gender"))
##df$gender <- unclass(df$gender)
#df$gender<-c(male=1,female=0)[df$gender]
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
# Build MSGL model
library(msgl)


fit_cv <- msgl::cv(as.matrix(df[,2:164]), df$frailty, 
                   fold = 10,alpha = 0.5, lambda = 0.1, use_parallel = TRUE, 
                   standardize = FALSE)
fit_cv

fit <- msgl::fit(as.matrix(sapply(df[,2:164], as.numeric)), df$frailty, 
                 alpha = 0.5, lambda = 0.1, standardize = FALSE)
fit

#find the best model index from cv
features(fit)[[best_model(fit_cv)]]
parameters(fit)[[best_model(fit_cv)]]
coef(fit, best_model(fit_cv))
coef = coef(fit, best_model(fit_cv))
coef = as.matrix(coef)
coef[which(coef == 0)] = NA
coef = data.frame(coef)
coef = coef[,-1]
coef$condition = rownames(coef)

library(reshape2)
coef = melt(coef)
coef$variable = as.character(coef$variable)


################################################################################
# MRMR Feature Selection
library(mRMRe)

f_data <- mRMR.data(data = df)
results <- mRMR.classic("mRMRe.Filter", data = f_data, target_indices = 1,
                        feature_count = 20)
solutions(results)
causality(results)
df[,c(6,10,35,37,96,102,108,113,116,120,124,125,127,132,134,138,149,153,155,164)]

#################
require(nnet)
multinom_fit <- multinom(frailty ~ `HABP2` + `MYBPHL` + `TJP3` + `TEX56P` + `CCL14`, data = df, model = TRUE)
multinom_fit <- multinom(frailty ~ `EFNA5` + `HABP2` + `MYBPHL` + `TJP3` + `TEX56P` + `ACRP1` + `LRTM1` + `CCL14` + `SEMA3E` + `PDF`, data = df, model = TRUE)
#multinom_fit <- multinom(frailty ~ `NME9` + `EFNA5` + `LINC02502` + `HABP2` + `GRIA2` + `KRTAP19-3` + `KCNJ16` + `MYBPHL` + `TJP3` + `LINC02660` + `MRAP` + `TEX56P` + `ALDH3A1` + `SAA4` + `ACRP1` + `LRTM1` + `CCL14` + `SEMA3E` + `STAP1` + `PDF`, data = df, model = TRUE)
multinom_fit <- multinom(frailty ~ `HABP2` + `TJP3` + `TEX56P` + `CCL14` + `PDF`, data = df, model = TRUE)


summary(multinom_fit)
round((1 - pnorm(abs(summary(multinom_fit)$coefficients/summary(multinom_fit)$standard.errors), 0, 1)) * 2, digits = 3)
exp(coef(multinom_fit))

head(pp <- fitted(multinom_fit))
library(sjPlot)
tab_model(multinom_fit, digits = 6)

prediction = predict(multinom_fit, newdata = df, "class")
sum(prediction == df$frailty)/nrow(df)



