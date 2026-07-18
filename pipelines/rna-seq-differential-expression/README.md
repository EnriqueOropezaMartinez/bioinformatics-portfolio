# RNA-seq Differential Expression Pipeline

## Overview

This repository contains a reproducible RNA-seq analysis workflow for identifying differentially expressed genes between biological conditions.

The pipeline processes sequencing data from raw reads to biological interpretation.

---

## Workflow

Raw FASTQ files
↓
Quality Control
(FastQC + MultiQC)
↓
Read preprocessing
(Trim Galore)
↓
Genome alignment
(STAR)
↓
Gene quantification
(featureCounts)
↓
Differential expression analysis
(DESeq2)
↓
Visualization
(PCA, Heatmap, Volcano Plot)
↓
Functional enrichment
(GO analysis)

---

## Input Data

Required files:

- FASTQ sequencing reads
- Reference genome
- Gene annotation file (GTF)

---

## Output

The pipeline generates:

- Quality control reports
- Gene count matrices
- Differential expression tables
- PCA plots
- Heatmaps
- Volcano plots
- Functional enrichment results

---

## Tools

| Tool | Purpose |
|---|---|
| FastQC | Sequencing quality control |
| MultiQC | QC report aggregation |
| Trim Galore | Read preprocessing |
| STAR | Genome alignment |
| featureCounts | Gene quantification |
| DESeq2 | Differential expression |
| R | Statistical analysis |
| Snakemake | Workflow management |

---

## Status

Pipeline under active development.