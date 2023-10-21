# Differential Expression Analysis of Frailty data using DESeq2
# Author: Shehbeel Arif
# NASA GeneLab Multi-Omics AWG
# Script adapted from ALSF's DE Analysis Tutorial (https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html)

## LOAD LIBRARIES
# Library for DE Analysis
library(DESeq2)
# Library data manipulation
library(dplyr)
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
analysis_dir <- file.path(root_dir, "analyses", "04-de-analysis")


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

#######
# Import metadata and data 
metadata <- readr::read_delim(metadata_file)
expression_df <- readr::read_csv(data_file)

#######
## PREPROCESS THE DATA
expression_df <- expression_df %>%
  distinct(gene_symbol, .keep_all = TRUE) %>% # Remove duplicate genes
  tibble::column_to_rownames("gene_symbol")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$sample_id)

# Check if the expression and metadata are in the same order
all.equal(colnames(expression_df), metadata$sample_id)

# metadata <- metadata %>%
#   tibble::column_to_rownames("sample_id")

# Make cluster a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    condition = factor(condition, levels = c("frail", "fit", "young"))
  )

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)


## Create DESeq2Dataset
# round all expression counts (if there are decimal values present in data)
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~condition
)

## Run Differential Expression Analysis
deseq_object <- DESeq(ddset)

# Export normalized read counts
normCounts <- counts(deseq_object, normalized=TRUE)
# Save results as CSV
write.csv(normCounts, file.path(results_dir, "mrna_normalized_counts.csv"))

# Extract results table
deseq_results <- results(deseq_object)
resultsNames(deseq_object)

# Use lfcShrink() function to obtain shrunken log fold change estimates based on 
# negative binomial distribution. This will add the estimates to your results table. 
# Using lfcShrink() can help decrease noise and preserve large differences between 
# groups (it requires that apeglm package be installed) (Zhu et al., Bioinformatics 2018).
# deseq_results <- lfcShrink(
#   deseq_object, # The original DESeq2 object after running DESeq()
#   coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
#   res = deseq_results # The original DESeq2 results table
# )

# Sort and filter DESeq2 results table and convert to dataframe

## Frailty vs Young
frailvyoung_results <- results(deseq_object, c("condition", "frail", "young"))
frailvyoung_df <- frailvyoung_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  frailvyoung_df,
  file.path(
    results_dir,
    "FrailvYoung_DEG_mrna.csv"
  )
)

## Frailty vs Fit
frailvfit_results <- results(deseq_object, c("condition", "frail", "fit"))
frailvfit_df <- frailvfit_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  frailvfit_df,
  file.path(
    results_dir,
    "FrailvFit_DEG_mrna.csv"
  )
)

## Fit vs Young
fitvyoung_results <- results(deseq_object, c("condition", "fit", "young"))
fitvyoung_df <- fitvyoung_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  fitvyoung_df,
  file.path(
    results_dir,
    "FitvYoung_DEG_mrna.csv"
  )
)

################################################################################
# Frailty vs All

# Make cluster a factor and set the levels appropriately
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    frailty = factor(frailty, levels = c("0", "1"))
  )

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)


## Create DESeq2Dataset
# round all expression counts (if there are decimal values present in data)
gene_matrix <- round(filtered_expression_df)

ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~frailty
)

## Run Differential Expression Analysis
deseq_object <- DESeq(ddset)

# # Export normalized read counts
# normCounts <- counts(deseq_object, normalized=TRUE)
# # Save results as CSV
# write.csv(normCounts, file.path(results_dir, "mrna_normalized_counts.csv"))

# Extract results table
deseq_results <- results(deseq_object)
resultsNames(deseq_object)


## Frailty vs All
frailvall_results <- results(deseq_object, c("frailty", "1", "0"))
frailvall_df <- frailvall_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # sort by statistic -- the highest values will be genes with higher expression
  dplyr::arrange(dplyr::desc(log2FoldChange))
# Save results as CSV
readr::write_csv(
  frailvall_df,
  file.path(
    results_dir,
    "FrailvAll_DEG_mrna.csv"
  )
)


################################################################################

## Create volcano plot
# Frail v All
frailvall_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  frailvall_df,
  lab = frailvall_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
frailvall_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = frailvall_volcano_plot,
  file.path(plots_dir, "FrailvAll_de_mrna_volcano_plot.tiff"),
  height = 8,
  width = 10
)

# Frail v Fit
frailvfit_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  frailvall_df,
  lab = frailvfit_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
frailvfit_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = frailvfit_volcano_plot,
  file.path(plots_dir, "FrailvFit_de_mrna_volcano_plot.tiff"),
  height = 8,
  width = 10
)

# Frail v Young
frailvyoung_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  frailvyoung_df,
  lab = frailvyoung_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
frailvyoung_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = frailvyoung_volcano_plot,
  file.path(plots_dir, "FrailvYoung_de_mrna_volcano_plot.tiff"),
  height = 8,
  width = 10
)

# Fit v Young
fitvyoung_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  fitvyoung_df,
  lab = fitvyoung_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
# Print out plot here
fitvyoung_volcano_plot
# Save volcano plot
ggplot2::ggsave(
  plot = fitvyoung_volcano_plot,
  file.path(plots_dir, "FitvYoung_de_mrna_volcano_plot.tiff"),
  height = 8,
  width = 10
)
