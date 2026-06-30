# Bulk RNA-seq Analysis of TRAIL- and Tumor-Conditioned Media-Induced Neutrophil Transcriptomic Changes

This repository contains an R-based bulk RNA-seq analysis workflow for evaluating transcriptional changes in neutrophils across TRAIL-associated and tumor-conditioned media experimental conditions.

The analysis is based on data associated with the following paper:

**TRAIL induces cytokine production via the NFkB2 pathway promoting neutrophil chemotaxis and neutrophil-mediated immune-suppression in triple negative breast cancer cells**
*Cancer Letters*, Volume 620, 2025, 217692
DOI: https://doi.org/10.1016/j.canlet.2025.217692

This project focuses specifically on the bulk RNA-seq differential expression and pathway enrichment portion of the study. It is not a full reproduction of all experiments in the paper.

---

## Project Overview

The goal of this project is to perform a reproducible bulk RNA-seq differential expression analysis using a gene count matrix and sample metadata.

The workflow includes:

* Loading and formatting a bulk RNA-seq count matrix
* Creating sample metadata for four experimental groups
* Running differential expression analysis with DESeq2
* Identifying significantly differentially expressed genes
* Generating volcano plots
* Generating top-50 differentially expressed gene heatmaps
* Generating a PCA plot
* Generating a sample-to-sample distance heatmap
* Running Hallmark gene set enrichment analysis using `fgsea`

---

## Experimental Conditions

The dataset contains 20 samples across four experimental conditions, with five samples per condition.

| Condition | Description                                      |
| --------- | ------------------------------------------------ |
| `SFM`     | Serum-free media                                 |
| `SFM-T`   | Serum-free media with TRAIL treatment            |
| `CM`      | Tumor-conditioned media                          |
| `CM-T`    | Conditioned media from TRAIL-treated tumor cells |

Note: In the paper, the `CM-T` condition may also be referred to as `T-CM`. This repository keeps the condition naming used in the local metadata and analysis script.

---

## Key Comparisons

This analysis focuses on two differential expression contrasts:

| Contrast       | Biological Question                                                                                                                                            |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SFM-T vs SFM` | What transcriptional changes occur in neutrophils exposed to TRAIL in serum-free media?                                                                        |
| `CM-T vs CM`   | What transcriptional changes occur in neutrophils exposed to conditioned media from TRAIL-treated tumor cells compared with untreated tumor-conditioned media? |

---

## Repository Structure

```text
.
├── README.md
├── scripts/
│   └── bulk_RNA.R
├── inputs/
│   ├── RNAseq.csv
│   ├── metadata.csv
│   └── h.all.v7.1.symbols.gmt
├── results/
│   ├── tables/
│   │   ├── results_sfm.csv
│   │   ├── results_cm.csv
│   │   ├── top_50_sfm.csv
│   │   └── top_50_cm.csv
│   └── figures/
│       ├── SFM_vs_SFMT_volcano_plot.png
│       ├── CM_vs_CMT_volcano_plot.png
│       ├── top50_DE_genes_SFM_vs_SFMT_heatmap.png
│       ├── top50_DE_genes_CM_vs_CMT_heatmap.png
│       ├── PCA_plot.png
│       ├── sample_distance_heatmap.png
│       └── fgsea_pathway_enrichment_SFM_vs_SFMT.png
└── README.md
```

---

## Input Files

### `RNAseq.csv`

Bulk RNA-seq count matrix used as input for DESeq2.

Expected format:

* Rows: genes
* Columns: samples
* Values: raw gene counts

The script converts the count data to a rounded count matrix before creating the DESeq2 object.

### `metadata.csv`

Sample metadata file describing the experimental condition for each sample.


### `h.all.v7.1.symbols.gmt`

Hallmark gene set file from MSigDB used for pathway enrichment analysis with `fgsea`.

---

## Methods

### 1. Data Loading and Formatting

The count matrix is loaded from `RNAseq.csv`, rounded, and converted into a matrix suitable for DESeq2 analysis.

Gene identifiers are formatted by extracting the gene symbol from row names when applicable.

### 2. Metadata Creation

Samples are assigned to one of four experimental conditions:

* `SFM`
* `SFM-T`
* `CM`
* `CM-T`

The metadata table is used as the sample annotation table for DESeq2.

### 3. Differential Expression Analysis

Differential expression analysis is performed using the `DESeq2` package.

The DESeq2 design formula used in this analysis is:

```r
design = ~ condition
```

Two contrasts are tested:

```r
SFM-T vs SFM
CM-T vs CM
```

### 4. Low-Expression Gene Filtering

Lowly expressed genes are filtered before differential expression analysis.

Filtering criterion:

```r
rowSums(counts(dds)) > 20
```

### 5. Significance Thresholds

Genes are considered significantly differentially expressed using the following criteria:

```text
adjusted p-value < 0.05
absolute log2 fold change >= 1
```

This threshold is used for this repository’s analysis and may differ from thresholds used in the original publication.

### 6. Visualization

The workflow generates the following visualizations:

* Volcano plots for each contrast
* Heatmaps of the top 50 differentially expressed genes
* PCA plot across all samples
* Sample-to-sample distance heatmap
* Hallmark pathway enrichment bar plot

### 7. Gene Set Enrichment Analysis

Gene set enrichment analysis is performed using `fgsea`.

The ranking metric used for enrichment analysis is:

```text
log2 fold change
```

Hallmark pathways are loaded from:

```text
h.all.v7.1.symbols.gmt
```

Significant pathways are filtered using:

```text
adjusted p-value < 0.05
```

---

## Output Files

### Differential Expression Results

| File              | Description                                              |
| ----------------- | -------------------------------------------------------- |
| `results_sfm.csv` | DESeq2 results for `SFM-T vs SFM`                        |
| `results_cm.csv`  | DESeq2 results for `CM-T vs CM`                          |
| `top_50_sfm.csv`  | Top 50 differentially expressed genes for `SFM-T vs SFM` |
| `top_50_cm.csv`   | Top 50 differentially expressed genes for `CM-T vs CM`   |


---

## Reproducibility Notes

* The count matrix should contain raw gene-level counts.
* The sample metadata should match the count matrix columns.
* The analysis uses DESeq2 normalization and statistical testing.
* Low-expression genes are filtered before running DESeq2.
* The volcano plots, heatmaps, PCA plot, distance heatmap, and GSEA plot are generated from the processed DESeq2 results.
* Results may vary slightly depending on package versions.

---

## Skills Demonstrated

This project demonstrates:

* Bulk RNA-seq analysis
* Differential expression analysis using DESeq2
* RNA-seq quality-control visualization
* PCA and sample distance analysis
* Heatmap generation
* Volcano plot generation
* Gene set enrichment analysis using fgsea
* R scripting for reproducible bioinformatics workflows
* Scientific documentation and GitHub project organization

---

## Citation

If using this dataset, please cite the original publication:

Huang et al. **TRAIL induces cytokine production via the NFkB2 pathway promoting neutrophil chemotaxis and neutrophil-mediated immune-suppression in triple negative breast cancer cells.** *Cancer Letters*. 2025;620:217692.
DOI: https://doi.org/10.1016/j.canlet.2025.217692

---

## Disclaimer

This repository is intended as a reproducible educational and portfolio analysis project based on publicly associated bulk RNA-seq data from the cited study. It does not reproduce all experiments, figures, or conclusions from the original publication.

