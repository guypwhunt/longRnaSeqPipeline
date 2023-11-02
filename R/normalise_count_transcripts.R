# Install Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2"))

install.packages(c("tidyr", "magrittr", "data.table", "dplyr"), dependencies = TRUE)

# Load required libraries
library(DESeq2)
library(tidyr)
library(magrittr)
library(dplyr)
library(data.table)

# Define the input and output directory
directory <- "alignment_output/count_transcripts/"

# Define a function to set row names from a column
column_to_rownames <- function(df) {
  row.names(df) <- df[, 1]
  df <- df[, -1]
  return(df)
}

# Define a function to select specific columns
select_columns <- function(df, columns) {
  df <- df[, columns]
  return(df)
}

# Define a function to filter rows with ENST IDs
select_rows <- function(df, pattern) {
  df <- df[grep(pattern, row.names(df)), ]
  return(df)
}

# Read the gene count data and set row names
gene_counts <- fread(
  paste0(directory, "barcode_counts_transcript.txt"),
  skip = 1,
  header = TRUE
) %>%
  as.data.frame() %>%
  column_to_rownames()

# Select columns containing count information and filter out non-ENST IDs
gene_counts <- gene_counts %>%
  select_columns(grep("bar", colnames(gene_counts)))

# Write the raw count data to a CSV file
fwrite(gene_counts,
       paste0(directory, "raw_counts_transcript.csv"),
       row.names = TRUE)

# Filter rows with ENST IDs
gene_counts <- gene_counts %>%
  select_rows("ENST")

# Write the filtered raw count data to a CSV file
fwrite(gene_counts,
       paste0(directory, "filtered_raw_counts_transcript.csv"),
       row.names = TRUE)

# Create a phenotype data frame
number_of_columns <- ncol(gene_counts)
halve_number_of_columns <-
  round(number_of_columns / 2)
remainder_number_of_columns <-
  number_of_columns - halve_number_of_columns
phe <- data.frame(id = colnames(gene_counts),
                  status = as.factor(c(
                    rep(1, halve_number_of_columns),
                    rep(2, remainder_number_of_columns)
                  )))

# Define the formula for DESeq2 analysis
formula <- as.formula("~ status")

# Create a DESeqDataSet
dds_genes <-
  DESeqDataSetFromMatrix(countData = gene_counts,
                         colData = phe,
                         design = formula)

# Estimate size factors
dds_genes <- estimateSizeFactors(dds_genes)
normalized_gene_counts <-
  counts(dds_genes, normalized = TRUE)

# Write the normalized count data to a CSV file
fwrite(
  normalized_gene_counts,
  paste0(directory, "normalized_counts_transcript.csv"),
  row.names = TRUE
)

# Filter low-abundance genes
idx <- rowSums(normalized_gene_counts >= 5) >= 10
dds_genes <- dds_genes[idx, ]
normalized_gene_counts <-
  counts(dds_genes, normalized = TRUE)

# Write the filtered normalized count data to a CSV file
fwrite(
  normalized_gene_counts,
  paste0(directory, "filtered_normalized_counts_transcript.csv"),
  row.names = TRUE
)
