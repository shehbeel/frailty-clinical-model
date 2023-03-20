# Convert Emsembl IDs to Gene Symbols
# Author: Shehbeel Arif

## LOAD LIBRARIES
# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# We will need this so we can use the pipe: %>%
library(magrittr)


## SET DIRECTORIES
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "01-convert-ensembl-ids")
data_dir <- file.path(root_dir, "data")

# Set output directories
results_dir <- file.path(analysis_dir, "results")

# Make output directories if they don't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Declare input file path
counts_file <- file.path(data_dir, "GSE144304_raw_counts.csv")
# Declare output file path
output_file <- file.path(results_dir, "GSE144304_raw_counts_with_genes.csv")

# Read in the counts file
counts <- readr::read_csv(counts_file) 

## GENE IDENTIFIER CONVERSION
# Look at what types of IDs are available to us
keytypes(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
# [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GENETYPE"     "GO"          
# [13] "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL" 
# [19] "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"      
# [25] "UCSCKG"       "UNIPROT"  

# Remove the values after decimal in the Ensembl IDs (represents different versions of the ID)
counts$ENSEMBL_ID <- gsub("\\.[0-9]*$", "", counts$ENSEMBL_ID)

# Map Ensembl IDs to Gene symbol
# This returns a named vector which we can convert to a data frame, where
# the keys (Ensembl IDs) are the names
gene_symbol_vector <- mapIds(
  # Replace with annotation package for the organism relevant to your data
  org.Hs.eg.db,
  # The vector of gene identifiers we want to map
  keys = counts$ENSEMBL_ID,
  # Replace with the type of gene identifiers in your data
  keytype = "ENSEMBL",
  # Replace with the type of gene identifiers you would like to map to
  column = "SYMBOL",
  # In the case of 1:many mappings, return the first one. This is default behavior!
  multiVals = "first"
)

# Create a dataframe containing both Ensembl IDs and gene symbols
# We would like a data frame we can join to the counts data
gene_key_df <- data.frame(
  ensembl_id = names(gene_symbol_vector),
  gene_symbol = gene_symbol_vector,
  stringsAsFactors = FALSE
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol))

# Merge Gene symbols with counts data
counts <- gene_key_df %>%
  # Using a left join removes the rows without gene symbols because those rows
  # have already been removed in `gene_key_df`
  dplyr::left_join(counts,
                   # The name of the column that contains the Ensembl gene IDs
                   # in the left data frame and right data frame
                   by = c("ensembl_id" = "ENSEMBL_ID")
  ) %>%
  # Drop 'ensembl_id' column
  dplyr::select(-ensembl_id)

# Export Gene counts data
readr::write_csv(counts, output_file)


