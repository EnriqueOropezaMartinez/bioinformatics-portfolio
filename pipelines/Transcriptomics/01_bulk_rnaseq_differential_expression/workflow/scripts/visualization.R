###############################################################################
# Bulk RNA-seq visualization
# Author: Enrique Oropeza Martínez
###############################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

results_file <- snakemake@input[["results"]]
vst_file <- snakemake@input[["vst_counts"]]
samples_file <- snakemake@input[["samples"]]

pca_file <- snakemake@output[["pca"]]
volcano_file <- snakemake@output[["volcano"]]
heatmap_file <- snakemake@output[["heatmap"]]
ma_file <- snakemake@output[["ma_plot"]]

design_variable <- snakemake@params[["design"]]
alpha_threshold <- as.numeric(snakemake@params[["alpha"]])
log2fc_threshold <- as.numeric(snakemake@params[["log2fc"]])
top_genes_number <- as.integer(snakemake@params[["top_genes"]])

dir.create(
  dirname(pca_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# Read files
# =============================================================================

results_table <- read_csv(
  results_file,
  show_col_types = FALSE
)

vst_table <- read_csv(
  vst_file,
  show_col_types = FALSE
)

sample_metadata <- read_csv(
  samples_file,
  show_col_types = FALSE
) %>%
  as.data.frame()

if (!"sample" %in% colnames(sample_metadata)) {
  stop("Metadata must contain a column named 'sample'.")
}

if (!design_variable %in% colnames(sample_metadata)) {
  stop(
    "Design variable not found in metadata: ",
    design_variable
  )
}

# =============================================================================
# Prepare VST matrix
# =============================================================================

vst_matrix <- vst_table %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

missing_samples <- setdiff(
  colnames(vst_matrix),
  sample_metadata$sample
)

if (length(missing_samples) > 0) {
  stop(
    "Samples missing from metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

sample_metadata <- sample_metadata[
  match(colnames(vst_matrix), sample_metadata$sample),
  ,
  drop = FALSE
]

rownames(sample_metadata) <- sample_metadata$sample

# =============================================================================
# PCA
# =============================================================================

pca_object <- prcomp(
  t(vst_matrix),
  center = TRUE,
  scale. = FALSE
)

percent_variance <- (
  pca_object$sdev^2 /
    sum(pca_object$sdev^2)
) * 100

pca_data <- as.data.frame(
  pca_object$x[, 1:2, drop = FALSE]
) %>%
  rownames_to_column("sample") %>%
  left_join(
    sample_metadata,
    by = "sample"
  )

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = .data[[design_variable]],
    label = sample
  )
) +
  geom_point(size = 4) +
  geom_text_repel(show.legend = FALSE) +
  labs(
    title = "Principal component analysis",
    x = paste0(
      "PC1 (",
      round(percent_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(percent_variance[2], 1),
      "%)"
    ),
    color = design_variable
  ) +
  theme_classic()

ggsave(
  filename = pca_file,
  plot = pca_plot,
  width = 8,
  height = 6
)

# =============================================================================
# Volcano plot
# =============================================================================

volcano_data <- results_table %>%
  mutate(
    minus_log10_padj = -log10(padj),
    minus_log10_padj = ifelse(
      is.infinite(minus_log10_padj),
      NA,
      minus_log10_padj
    )
  )

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = regulation
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(
      -log2fc_threshold,
      log2fc_threshold
    ),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(alpha_threshold),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano plot",
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "Regulation"
  ) +
  theme_classic()

ggsave(
  filename = volcano_file,
  plot = volcano_plot,
  width = 8,
  height = 6
)

# =============================================================================
# Heatmap
# =============================================================================

top_genes <- results_table %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = top_genes_number) %>%
  pull(gene_id)

top_genes <- intersect(
  top_genes,
  rownames(vst_matrix)
)

if (length(top_genes) < 2) {
  stop("Not enough genes available to generate the heatmap.")
}

heatmap_matrix <- vst_matrix[
  top_genes,
  ,
  drop = FALSE
]

annotation_columns <- sample_metadata[
  ,
  design_variable,
  drop = FALSE
]

pdf(
  heatmap_file,
  width = 9,
  height = 10
)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_columns,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  main = paste(
    "Top",
    length(top_genes),
    "differentially expressed genes"
  )
)

dev.off()

# =============================================================================
# MA plot
# =============================================================================

ma_plot <- ggplot(
  results_table,
  aes(
    x = baseMean,
    y = log2FoldChange,
    color = regulation
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.3,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = c(
      -log2fc_threshold,
      log2fc_threshold
    ),
    linetype = "dashed"
  ) +
  scale_x_log10() +
  labs(
    title = "MA plot",
    x = "Mean normalized expression",
    y = "log2 fold change",
    color = "Regulation"
  ) +
  theme_classic()

ggsave(
  filename = ma_file,
  plot = ma_plot,
  width = 8,
  height = 6
)

message("Visualization completed successfully.")