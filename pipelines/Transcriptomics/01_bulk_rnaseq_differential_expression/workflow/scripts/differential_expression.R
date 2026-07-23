###############################################################################
# Differential expression analysis with DESeq2
# Bulk RNA-seq Snakemake pipeline
# Author: Enrique Oropeza Martínez
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

# =============================================================================
# Snakemake inputs, outputs and parameters
# =============================================================================

counts_file <- snakemake@input[["counts"]]
samples_file <- snakemake@input[["samples"]]

results_file <- snakemake@output[["results"]]
normalized_counts_file <- snakemake@output[["normalized_counts"]]
vst_counts_file <- snakemake@output[["vst_counts"]]
summary_file <- snakemake@output[["summary"]]

design_variable <- snakemake@params[["design"]]
reference_condition <- snakemake@params[["reference"]]
contrast_condition <- snakemake@params[["contrast"]]
alpha_threshold <- as.numeric(snakemake@params[["alpha"]])
log2fc_threshold <- as.numeric(snakemake@params[["log2fc"]])

# =============================================================================
# Helper functions
# =============================================================================

clean_sample_name <- function(path) {
  sample_name <- basename(path)
  sample_name <- sub("\\.sorted\\.bam$", "", sample_name)
  return(sample_name)
}

write_analysis_summary <- function(
    output_file,
    total_genes,
    tested_genes,
    significant_genes,
    upregulated_genes,
    downregulated_genes
) {
  summary_lines <- c(
    "DESeq2 differential expression summary",
    "======================================",
    paste("Design variable:", design_variable),
    paste("Reference condition:", reference_condition),
    paste("Contrast condition:", contrast_condition),
    paste("Adjusted p-value threshold:", alpha_threshold),
    paste("Absolute log2 fold-change threshold:", log2fc_threshold),
    "",
    paste("Total genes in count matrix:", total_genes),
    paste("Genes tested by DESeq2:", tested_genes),
    paste("Significant genes:", significant_genes),
    paste("Upregulated genes:", upregulated_genes),
    paste("Downregulated genes:", downregulated_genes)
  )

  writeLines(summary_lines, output_file)
}

# =============================================================================
# Read featureCounts output
# =============================================================================

message("Reading featureCounts matrix: ", counts_file)

featurecounts_data <- read.delim(
  counts_file,
  header = TRUE,
  comment.char = "#",
  check.names = FALSE
)

required_annotation_columns <- c(
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(featurecounts_data)
)

if (length(missing_annotation_columns) > 0) {
  stop(
    "Missing expected featureCounts columns: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

count_matrix <- featurecounts_data %>%
  select(-all_of(required_annotation_columns)) %>%
  as.data.frame()

rownames(count_matrix) <- featurecounts_data$Geneid

colnames(count_matrix) <- vapply(
  colnames(count_matrix),
  clean_sample_name,
  character(1)
)

count_matrix <- round(as.matrix(count_matrix))

storage.mode(count_matrix) <- "integer"

# Remove genes with zero counts across all samples
count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]

# =============================================================================
# Read and validate sample metadata
# =============================================================================

message("Reading sample metadata: ", samples_file)

sample_metadata <- read_csv(
  samples_file,
  show_col_types = FALSE
) %>%
  as.data.frame()

if (!"sample" %in% colnames(sample_metadata)) {
  stop("The metadata file must contain a column named 'sample'.")
}

if (!design_variable %in% colnames(sample_metadata)) {
  stop(
    "The design variable '",
    design_variable,
    "' was not found in the metadata."
  )
}

if (anyDuplicated(sample_metadata$sample)) {
  stop("Duplicated sample names were found in the metadata.")
}

missing_metadata_samples <- setdiff(
  colnames(count_matrix),
  sample_metadata$sample
)

missing_count_samples <- setdiff(
  sample_metadata$sample,
  colnames(count_matrix)
)

if (length(missing_metadata_samples) > 0) {
  stop(
    "Samples in the count matrix missing from metadata: ",
    paste(missing_metadata_samples, collapse = ", ")
  )
}

if (length(missing_count_samples) > 0) {
  stop(
    "Samples in metadata missing from count matrix: ",
    paste(missing_count_samples, collapse = ", ")
  )
}

# Reorder metadata to match the count matrix
sample_metadata <- sample_metadata[
  match(colnames(count_matrix), sample_metadata$sample),
  ,
  drop = FALSE
]

rownames(sample_metadata) <- sample_metadata$sample

sample_metadata[[design_variable]] <- factor(
  sample_metadata[[design_variable]]
)

available_conditions <- levels(sample_metadata[[design_variable]])

if (!reference_condition %in% available_conditions) {
  stop(
    "Reference condition not found in metadata: ",
    reference_condition
  )
}

if (!contrast_condition %in% available_conditions) {
  stop(
    "Contrast condition not found in metadata: ",
    contrast_condition
  )
}

sample_metadata[[design_variable]] <- relevel(
  sample_metadata[[design_variable]],
  ref = reference_condition
)

# =============================================================================
# DESeq2 analysis
# =============================================================================

message(
  "Running DESeq2 contrast: ",
  contrast_condition,
  " versus ",
  reference_condition
)

design_formula <- as.formula(
  paste("~", design_variable)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = design_formula
)

# Keep genes with at least 10 reads across all samples
dds <- dds[rowSums(counts(dds)) >= 10, ]

if (nrow(dds) == 0) {
  stop("No genes passed the minimum count filter.")
}

dds <- DESeq(dds)

deseq_results <- results(
  dds,
  contrast = c(
    design_variable,
    contrast_condition,
    reference_condition
  ),
  alpha = alpha_threshold
)

# =============================================================================
# Prepare DESeq2 results
# =============================================================================

results_table <- as.data.frame(deseq_results) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    regulation = case_when(
      !is.na(padj) &
        padj < alpha_threshold &
        log2FoldChange >= log2fc_threshold ~ "Upregulated",

      !is.na(padj) &
        padj < alpha_threshold &
        log2FoldChange <= -log2fc_threshold ~ "Downregulated",

      TRUE ~ "Not_significant"
    )
  ) %>%
  arrange(padj)

# =============================================================================
# Normalized and VST-transformed counts
# =============================================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

vst_object <- vst(
  dds,
  blind = FALSE
)

vst_counts <- assay(vst_object) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

# =============================================================================
# Write outputs
# =============================================================================

dir.create(
  dirname(results_file),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(results_table, results_file)
write_csv(normalized_counts, normalized_counts_file)
write_csv(vst_counts, vst_counts_file)

total_genes <- nrow(count_matrix)
tested_genes <- nrow(results_table)

significant_genes <- sum(
  results_table$regulation != "Not_significant"
)

upregulated_genes <- sum(
  results_table$regulation == "Upregulated"
)

downregulated_genes <- sum(
  results_table$regulation == "Downregulated"
)

write_analysis_summary(
  output_file = summary_file,
  total_genes = total_genes,
  tested_genes = tested_genes,
  significant_genes = significant_genes,
  upregulated_genes = upregulated_genes,
  downregulated_genes = downregulated_genes
)

message("DESeq2 analysis completed successfully.")