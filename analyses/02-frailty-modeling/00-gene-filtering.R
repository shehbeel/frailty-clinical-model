# Filtering out lowly expressed genes in Frailty data
# Author: Shehbeel Arif
# NASA GeneLab Multi-Omics AWG

## LOAD LIBRARIES
# Library for DE Analysis
library(DESeq2)
# Library data manipulation
#library(dplyr)
# Library for Plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)

# Package that contains MSigDB gene sets in tidy format
library(msigdbr)

# Homo sapiens annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)

# Set seed because jitter plot function involves randomness
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
data_dir <- file.path(root_dir, "analyses", "01-convert-ensembl-ids", "results")

# Declare input file paths
data_file <- file.path(data_dir, "GSE144304_raw_counts_with_genes.csv")

# Load data
expression_df <- readr::read_csv(data_file) 

# Drop duplicate genes
expression_df <- expression_df %>% dplyr::distinct(gene_symbol, .keep_all = TRUE)

# Make genes as column names
expression_df <- expression_df %>%
  tibble::column_to_rownames("gene_symbol")


# Remove lowly expressing genes
filtered_expression_df <- round(expression_df) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50)

# Save results as CSV
write.csv(filtered_expression_df, file.path(results_dir, "GSE144304_raw_counts_with_genes_filtered50.csv"))


