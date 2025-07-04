# Bulk RNA-seq Differential Expression and GSEA Analysis

## 📌 Overview  
This repository contains an R script for performing bulk RNA-seq differential expression analysis and gene set enrichment analysis (GSEA) using **DESeq2** and **fgsea**. The workflow includes normalization, filtering, differential expression, visualization (heatmaps, volcano plots, PCA), and GSEA based on Hallmark gene sets.

## 🧬 Requirements  
Make sure you have the following R packages installed:
```r
library(DESeq2)
library(fgsea)
library(pheatmap)
library(ggplot2)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(EnhancedVolcano)

## 📂 Input Data

RNAseq.csv — A count matrix file where rows represent genes (with an underscore prefix) and columns represent samples.
h.all.v7.1.symbols.gmt — Hallmark gene sets file for GSEA (downloaded from MSigDB).
## 📝 Steps in the Analysis

### 1️⃣ Data Loading and Preprocessing
Reads a CSV file with raw counts.
Converts counts to matrix form and updates gene identifiers.
Defines sample conditions and constructs metadata.
Filters out lowly expressed genes (row sums ≤ 20).
### 2️⃣ Differential Expression (DE) Analysis
Performs DE analysis between:
SFM-T vs SFM
CM-T vs CM
Identifies significantly upregulated and downregulated genes (|log2FC| ≥ 1 and FDR < 0.05).
### 3️⃣ Visualization
Heatmaps of the top 50 differentially expressed genes for each contrast.
Volcano plots highlighting genes with |log2FC| > 2 and FDR < 0.05.
PCA plot to explore sample clustering.
Sample distance heatmap to assess sample similarity.
### 4️⃣ Gene Set Enrichment Analysis (GSEA)
Ranks genes by log2 fold change.
Runs GSEA using Hallmark pathways.
Plots normalized enrichment scores for significant pathways (FDR < 0.05).

## 📊 Example Plots
Volcano plots of DE genes
Heatmaps of top DE genes
PCA plot of samples
Bar plot of enriched pathways

## 🚀 How to Run

1️⃣ Update file paths to your local data files (e.g., RNAseq.csv, h.all.v7.1.symbols.gmt).
2️⃣ Run the R script in your preferred R environment (RStudio recommended).

## ⚠️ Notes

Ensure your sample conditions align with the condition vector.
Adjust filtering thresholds or plotting parameters as needed for your dataset.
The script uses 1000 permutations for GSEA; you can increase this for more robust results.

## Credits

This pipeline integrates popular R packages (DESeq2, fgsea, EnhancedVolcano, etc.) and standard best practices for bulk RNA-seq analysis.
