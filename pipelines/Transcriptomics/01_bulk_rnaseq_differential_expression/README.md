# Bulk RNA-seq Differential Expression Pipeline

A reproducible and modular Snakemake workflow for bulk RNA-seq differential expression analysis from raw FASTQ files to biological interpretation.

---

## Overview

This workflow performs an end-to-end RNA sequencing analysis including:

- Quality assessment
- Adapter trimming
- Genome alignment
- Gene quantification
- Differential expression analysis
- Exploratory visualization
- Functional enrichment analysis

The pipeline was designed following reproducible bioinformatics practices using Snakemake and Conda.

---

## Workflow

```
FASTQ
   │
   ▼
FastQC
   │
   ▼
Trim Galore
   │
   ▼
FastQC
   │
   ▼
MultiQC
   │
   ▼
STAR
   │
   ▼
Sorted BAM
   │
   ▼
featureCounts
   │
   ▼
DESeq2
   │
   ├── Normalized counts
   ├── VST counts
   ├── Differential expression
   └── Summary
   │
   ▼
Visualization
   ├── PCA
   ├── Volcano
   ├── Heatmap
   └── MA Plot
   │
   ▼
GO Biological Process
KEGG Pathways
```

---

# Features

✔ Quality control with FastQC

✔ Adapter trimming using Trim Galore

✔ STAR genome alignment

✔ BAM sorting and indexing

✔ Alignment statistics

✔ Gene quantification with featureCounts

✔ Differential expression analysis using DESeq2

✔ Variance stabilizing transformation (VST)

✔ PCA

✔ Volcano plots

✔ Heatmaps

✔ MA plots

✔ GO enrichment

✔ KEGG enrichment

✔ Independent Conda environments

✔ Fully reproducible Snakemake workflow

---

# Directory structure

```
project/
│
├── config/
├── data/
├── workflow/
│   ├── rules/
│   ├── envs/
│   └── scripts/
├── results/
├── Snakefile
└── README.md
```

---

# Input

Raw paired-end FASTQ files

```
Control_1_R1.fastq.gz
Control_1_R2.fastq.gz
Tumor_1_R1.fastq.gz
Tumor_1_R2.fastq.gz
```

Sample sheet

```
sample,condition
Control_1,Control
Control_2,Control
Tumor_1,Tumor
Tumor_2,Tumor
```

---

# Software

| Tool | Version |
|-------|----------|
| Snakemake | |
| FastQC | |
| Trim Galore | |
| STAR | |
| SAMtools | |
| featureCounts | |
| DESeq2 | |
| clusterProfiler | |

---

# Installation

```bash
git clone ...
cd ...
```

Create Conda environments automatically

```
snakemake --use-conda
```

---

# Run

Dry run

```bash
snakemake --dry-run --use-conda
```

Execute

```bash
snakemake \
--use-conda \
--cores 16
```

---

# Outputs

Quality Control

```
results/qc/
```

Alignment

```
results/alignment/
```

Gene Counts

```
results/counts/
```

Differential Expression

```
results/differential_expression/
```

Figures

```
results/figures/
```

Functional Enrichment

```
results/enrichment/
```

---

# Differential Expression

Genes are considered significant when:

- Adjusted p-value < 0.05
- |log2FoldChange| ≥ 1

---

# Functional Enrichment

GO Biological Process

KEGG

clusterProfiler

---

# Validation

Current status

✅ Pipeline implemented

✅ Syntax validated

⬜ Tested on Linux

⬜ Tested using public RNA-seq datasets

⬜ Version 1.0

---

# Future Improvements

- Docker support
- Nextflow implementation
- Salmon pseudoalignment
- edgeR support
- limma-voom support
- GSEA
- Reactome analysis
- Automatic report generation

---

# Citation

If you use this workflow please cite ...

---

# Author

Enrique Oropeza Martínez

M.Sc. Candidate in Genomic Sciences

Bioinformatics • Cancer Genomics • Transcriptomics
