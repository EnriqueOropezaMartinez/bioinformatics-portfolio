###############################################################################
# Functional enrichment analysis
# Bulk RNA-seq Snakemake pipeline
# Author: Enrique Oropeza Martínez
###############################################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# =============================================================================
# Snakemake inputs, outputs and parameters
# =============================================================================

results_file <- snakemake@input[["results"]]

go_up_file <- snakemake@output[["go_up"]]
go_down_file <- snakemake@output[["go_down"]]
kegg_file <- snakemake@output[["kegg"]]

go_up_plot_file <- snakemake@output[["go_up_plot"]]
go_down_plot_file <- snakemake@output[["go_down_plot"]]
kegg_plot_file <- snakemake@output[["kegg_plot"]]

summary_file <- snakemake@output[["summary"]]

alpha_threshold <- as.numeric(
  snakemake@params[["alpha"]]
)

log2fc_threshold <- as.numeric(
  snakemake@params[["log2fc"]]
)

organism_code <- snakemake@params[["organism"]]
input_id_type <- snakemake@params[["id_type"]]

show_categories <- as.integer(
  snakemake@params[["show_categories"]]
)

# =============================================================================
# Helper functions
# =============================================================================

write_empty_table <- function(path) {
  write_csv(
    tibble(
      ID = character(),
      Description = character(),
      GeneRatio = character(),
      BgRatio = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric(),
      geneID = character(),
      Count = integer()
    ),
    path
  )
}

save_empty_plot <- function(path, title_text, message_text) {
  empty_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = message_text,
      size = 5
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title_text) +
    theme_void()

  ggsave(
    filename = path,
    plot = empty_plot,
    width = 9,
    height = 7
  )
}

run_go_enrichment <- function(entrez_ids) {
  if (length(entrez_ids) == 0) {
    return(NULL)
  }

  enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = alpha_threshold,
    qvalueCutoff = alpha_threshold,
    readable = TRUE
  )
}

# =============================================================================
# Read differential-expression results
# =============================================================================

message(
  "Reading DESeq2 results: ",
  results_file
)

results_table <- read_csv(
  results_file,
  show_col_types = FALSE
)

required_columns <- c(
  "gene_id",
  "log2FoldChange",
  "padj",
  "regulation"
)

missing_columns <- setdiff(
  required_columns,
  colnames(results_table)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in DESeq2 results: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Remove Ensembl version suffixes when present
results_table <- results_table %>%
  mutate(
    clean_gene_id = sub(
      "\\.[0-9]+$",
      "",
      gene_id
    )
  )

# =============================================================================
# Select significant genes
# =============================================================================

significant_results <- results_table %>%
  filter(
    !is.na(padj),
    padj < alpha_threshold,
    abs(log2FoldChange) >= log2fc_threshold
  )

upregulated_ids <- significant_results %>%
  filter(log2FoldChange >= log2fc_threshold) %>%
  pull(clean_gene_id) %>%
  unique()

downregulated_ids <- significant_results %>%
  filter(log2FoldChange <= -log2fc_threshold) %>%
  pull(clean_gene_id) %>%
  unique()

all_significant_ids <- significant_results %>%
  pull(clean_gene_id) %>%
  unique()

# =============================================================================
# Convert identifiers to Entrez
# =============================================================================

convert_to_entrez <- function(gene_ids) {
  if (length(gene_ids) == 0) {
    return(character())
  }

  converted <- bitr(
    gene_ids,
    fromType = input_id_type,
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  unique(converted$ENTREZID)
}

upregulated_entrez <- convert_to_entrez(
  upregulated_ids
)

downregulated_entrez <- convert_to_entrez(
  downregulated_ids
)

all_significant_entrez <- convert_to_entrez(
  all_significant_ids
)

# =============================================================================
# GO Biological Process
# =============================================================================

go_up <- run_go_enrichment(
  upregulated_entrez
)

go_down <- run_go_enrichment(
  downregulated_entrez
)

# =============================================================================
# KEGG
# =============================================================================

if (length(all_significant_entrez) > 0) {
  kegg_results <- enrichKEGG(
    gene = all_significant_entrez,
    organism = organism_code,
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = alpha_threshold,
    qvalueCutoff = alpha_threshold
  )
} else {
  kegg_results <- NULL
}

# =============================================================================
# Create output directories
# =============================================================================

dir.create(
  dirname(go_up_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(go_up_plot_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# Write GO outputs and plots
# =============================================================================

if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
  write_csv(
    as.data.frame(go_up),
    go_up_file
  )

  go_up_plot <- dotplot(
    go_up,
    showCategory = show_categories
  ) +
    labs(
      title = "GO biological processes: upregulated genes"
    )

  ggsave(
    filename = go_up_plot_file,
    plot = go_up_plot,
    width = 9,
    height = 7
  )
} else {
  write_empty_table(go_up_file)

  save_empty_plot(
    go_up_plot_file,
    "GO biological processes: upregulated genes",
    "No significantly enriched terms"
  )
}

if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
  write_csv(
    as.data.frame(go_down),
    go_down_file
  )

  go_down_plot <- dotplot(
    go_down,
    showCategory = show_categories
  ) +
    labs(
      title = "GO biological processes: downregulated genes"
    )

  ggsave(
    filename = go_down_plot_file,
    plot = go_down_plot,
    width = 9,
    height = 7
  )
} else {
  write_empty_table(go_down_file)

  save_empty_plot(
    go_down_plot_file,
    "GO biological processes: downregulated genes",
    "No significantly enriched terms"
  )
}

# =============================================================================
# Write KEGG output and plot
# =============================================================================

if (
  !is.null(kegg_results) &&
  nrow(as.data.frame(kegg_results)) > 0
) {
  write_csv(
    as.data.frame(kegg_results),
    kegg_file
  )

  kegg_plot <- dotplot(
    kegg_results,
    showCategory = show_categories
  ) +
    labs(
      title = "KEGG pathways: differentially expressed genes"
    )

  ggsave(
    filename = kegg_plot_file,
    plot = kegg_plot,
    width = 9,
    height = 7
  )
} else {
  write_empty_table(kegg_file)

  save_empty_plot(
    kegg_plot_file,
    "KEGG pathways: differentially expressed genes",
    "No significantly enriched pathways"
  )
}

# =============================================================================
# Write analysis summary
# =============================================================================

go_up_terms <- if (
  !is.null(go_up)
) {
  nrow(as.data.frame(go_up))
} else {
  0
}

go_down_terms <- if (
  !is.null(go_down)
) {
  nrow(as.data.frame(go_down))
} else {
  0
}

kegg_terms <- if (
  !is.null(kegg_results)
) {
  nrow(as.data.frame(kegg_results))
} else {
  0
}

summary_lines <- c(
  "Functional enrichment summary",
  "=============================",
  paste("Input identifier type:", input_id_type),
  paste("Organism code:", organism_code),
  paste("Adjusted p-value threshold:", alpha_threshold),
  paste("Absolute log2 fold-change threshold:", log2fc_threshold),
  "",
  paste("Significant genes:", length(all_significant_ids)),
  paste("Upregulated genes:", length(upregulated_ids)),
  paste("Downregulated genes:", length(downregulated_ids)),
  "",
  paste(
    "Mapped significant Entrez IDs:",
    length(all_significant_entrez)
  ),
  paste(
    "Mapped upregulated Entrez IDs:",
    length(upregulated_entrez)
  ),
  paste(
    "Mapped downregulated Entrez IDs:",
    length(downregulated_entrez)
  ),
  "",
  paste("Significant GO terms, upregulated:", go_up_terms),
  paste("Significant GO terms, downregulated:", go_down_terms),
  paste("Significant KEGG pathways:", kegg_terms)
)

writeLines(
  summary_lines,
  summary_file
)

message(
  "Functional enrichment completed successfully."
)
