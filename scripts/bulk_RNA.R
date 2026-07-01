##Bulk RNAseq analysis###

#Load required libraries
library(DESeq2)
library(fgsea)
library(pheatmap)
library(ggplot2)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(EnhancedVolcano)

#Load Data: Count matrix

#read in the csv file
countdata <- read.csv("../inputs/RNAseq.csv", row.names = 1, header = TRUE)

# Convert count data to a matrix
countdata <- as.matrix(round(countdata))
colnames(countdata)

#Define experimental conditions
condition <- factor(c(rep("SFM", 5), rep("SFM-T", 5), rep("CM", 5), rep("CM-T", 5)))

metadata <- read.csv("../inputs/metadata.csv", row.names = 1)
metadata$condition <- factor(metadata$condition)

countdata <- countdata[, rownames(metadata)]


# Update rownames to keep only the gene symbol after the underscore
rownames(countdata) <- sub(".*_", "", rownames(countdata))
head(rownames(countdata))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~ condition)

# Filter out lowly expressed genes (row sums > 20)
dds <- dds[rowSums(counts(dds)) > 20, ]

# Normalization
dds <- estimateSizeFactors(dds)
dds

# Run DESeq2
dds <- DESeq(dds)

#Extract results

res_sfm <- results(dds, contrast = c("condition", "SFM-T", "SFM")) #SFM-T compared to SFM
res_cm <- results(dds, contrast = c("condition", "CM-T", "CM"))  #CM-T compared to CM


# Identify upregulated and downregulated genes
# Remove rows with NA values in log2FoldChange or padj columns
res_cleaned_sfm <- res_sfm[!is.na(res_sfm$log2FoldChange) & !is.na(res_sfm$padj), ]
res_cleaned_cm <- res_cm[!is.na(res_cm$log2FoldChange) & !is.na(res_cm$padj), ]

dir.create("../results/tables", recursive = TRUE, showWarnings = FALSE)
# Now apply the filtering conditions to sfm and sfm-t
res_filtered_sfm <- res_cleaned_sfm[abs(res_cleaned_sfm$log2FoldChange) >= 1 & res_cleaned_sfm$padj < 0.05, ]
write.csv(res_filtered_sfm, "../results/tables/deseq2_results_SFM-T_vs_SFM.csv", row.names = TRUE)

upregulated_genes_sfm <- res_filtered_sfm[res_filtered_sfm$log2FoldChange > 0,]
downregulated_genes_sfm <- res_filtered_sfm[res_filtered_sfm$log2FoldChange < 0,]
# Number of upregulated and downregulated genes
cat("Upregulated genes SFM: ", nrow(upregulated_genes_sfm))
cat("Downregulated genes SFM: ", nrow(downregulated_genes_sfm))

# Now apply the filtering conditions to cm and cm-t
res_filtered_cm <- res_cleaned_cm[abs(res_cleaned_cm$log2FoldChange) >= 1 & res_cleaned_cm$padj < 0.05, ]
write.csv(res_filtered_cm, "../results/tables/deseq2_results_CM-T_vs_CM.csv", row.names = TRUE)

upregulated_genes_cm <- res_filtered_cm[res_filtered_cm$log2FoldChange > 0,]
downregulated_genes_cm <- res_filtered_cm[res_filtered_cm$log2FoldChange < 0,]
# Number of upregulated and downregulated genes
cat("Upregulated genes CM: ", nrow(upregulated_genes_cm))
cat("Downregulated genes CM: ", nrow(downregulated_genes_cm))

#Heatmap

# Extract top 50 differentially expressed genes in SFM and SFM-T samples
top_genes_sfm <- rownames(res_sfm)[order(res_sfm$padj)][1:50]
write.csv(top_genes_sfm, "../results/tables/top50_DE_genes_SFM-T_vs_SFM.csv", row.names = TRUE)

# Normalize counts for visualization
rlog_counts_sfm <- rlog(dds, blind = FALSE)

#Get SFM and SFM-T samples 
sfm_samples <- rownames(metadata)[metadata$condition %in% c("SFM", "SFM-T")]

# Extract normalized counts for top genes
top_gene_counts_sfm <- assay(rlog_counts_sfm)[top_genes_sfm, sfm_samples ]

dir.create("../results/figures", recursive = TRUE, showWarnings = FALSE)
# Generate heatmap
pheatmap(top_gene_counts_sfm,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE, fontsize = 6,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Heatmap of Top 50 Differentially Expressed Genes in SFM and SFM-T",
         filename = "../results/figures/heatmap_top50_SFM-T_vs_SFM.png")

# Extract top 50 differentially expressed genes in CM and CM-T samples
top_genes_cm <- rownames(res_cm)[order(res_cm$padj)][1:50]
write.csv(top_genes_cm, "../results/tables/top50_DE_genes_CM-T_vs_CM.csv", row.names = TRUE)

# Normalize counts for visualization
rlog_counts_cm <- rlog(dds, blind = FALSE)

#Get SFM and SFM-T samples 
cm_samples <- rownames(metadata)[metadata$condition %in% c("CM", "CM-T")]

# Extract normalized counts for top genes
top_gene_counts_cm <- assay(rlog_counts_cm)[top_genes_cm, cm_samples ]

# Generate heatmap
pheatmap(top_gene_counts_cm,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE, fontsize = 6,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Heatmap of Top 50 Differentially Expressed Genes in CM and CM-T",
         filename = "../results/figures/heatmap_top50_CM-T_vs_CM.png")

# Volcano Plot
volcano_sfm <- EnhancedVolcano(
  res_sfm,
  lab = rownames(res_sfm),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 2,
  pointSize = 2.0,
  labSize = 3.0,
  col = c("lightskyblue4", "steelblue1", "rosybrown1", "hotpink1"),
  title = "Volcano Plot",
  subtitle = "FDR < 0.05, |log2FC| > 2"
)

ggsave(
  filename = "../results/figures/volcano_SFM-T_vs_SFM.png",
  plot = volcano_sfm,
  width = 8,
  height = 6,
  dpi = 300
)

volcano_cm <- EnhancedVolcano(res_cm,
                              lab = rownames(res_cm), # Gene labels
                              x = 'log2FoldChange', # X-axis: log2 fold change
                              y = 'padj', # Y-axis: adjusted p-value (FDR)
                              pCutoff = 0.05, # FDR cutoff of 0.05
                              FCcutoff = 2, # Fold-change cutoff of 4 (log2FC = 2)
                              pointSize = 2.0, # Size of points
                              labSize = 3.0, # Size of labels
                              col = c('lightskyblue4', 'steelblue1', 'rosybrown1', 'hotpink1'), # Custom colors
                              title = "Volcano Plot", # Plot title
                              subtitle = "FDR < 0.05, |log2FC| > 2" # Plot subtitle
)

ggsave(
  filename = "../results/figures/volcano_CM-T_vs_CM.png",
  plot = volcano_cm,
  width = 8,
  height = 6,
  dpi = 300
)

# PCA Plot

# Perform PCA on rlog-transformed data
rlog_counts <- rlog(dds, blind = FALSE)
pca_data <- prcomp(t(assay(rlog_counts)))

# Prepare data for plotting
pca_df <- as.data.frame(pca_data$x)
pca_df$condition <- metadata$condition
pca_df$sample_id <- rownames(metadata)

# Generate PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  scale_color_brewer(palette = "Set1")

ggsave(
  filename = "../results/figures/PCA_all_samples.png",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)


#Sample relationship

#construct distance matrix
sampleDists <- dist(t(assay(rlog_counts)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata$condition, rownames(metadata), sep = "_")
colnames(sampleDistMatrix) <- paste(metadata$condition, rownames(metadata), sep = "_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  filename = "../results/figures/sample_distance_heatmap.png",
  width = 8,
  height = 6
)


# Prepare data for GSEA
res1_sfm <- as.data.frame(res_sfm)
res1_sfm <- res1_sfm %>%
  tibble::rownames_to_column(var = "GeneSymbol")

res2_sfm <- res1_sfm %>% 
  dplyr::select(GeneSymbol, log2FoldChange)

ranks_sfm <- tibble::deframe(res2_sfm)
ranks_sfm <- sort(ranks_sfm, decreasing = TRUE)

# Load Hallmark gene sets
pathways.hallmark <- gmtPathways("../inputs/h.all.v7.1.symbols.gmt")

# Run GSEA
fgseaRes_sfm <- fgsea(pathways = pathways.hallmark, stats = ranks_sfm)

# Tidy GSEA results
fgseaResTidy_sfm <- fgseaRes_sfm %>%
  as_tibble() %>%
  arrange(desc(NES))

# Filter significant pathways
fgseaResTidy1_sfm <- filter(fgseaResTidy_sfm, fgseaResTidy_sfm$padj < 0.05)
fgseaResTidy1_sfm <- fgseaResTidy1_sfm %>% 
  mutate(variables = case_when(
    NES > 0 ~ "up", 
    NES < 0  ~ "down"))

fgseaResTidy1_sfm
# Plot GSEA results

fgsea_plot1 <- ggplot(fgseaResTidy1_sfm, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = variables)) +
  coord_flip() +
  scale_fill_manual(values = c("#4d9221", "#c51b7d")) +
  labs(
    x = "Pathway",
    y = "Normalized Enrichment Score",
    title = "FGSEA Pathway Enrichment SFM-T vs SFM"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = "../results/figures/gsea_hallmark_SFM-T_vs_SFM.png",
  plot = fgsea_plot1,
  width = 9,
  height = 7,
  dpi = 300
)

res1_cm <- as.data.frame(res_cm)
res1_cm <- res1_cm %>%
  tibble::rownames_to_column(var = "GeneSymbol")

res2_cm <- res1_cm %>% 
  dplyr::select(GeneSymbol, log2FoldChange)

ranks_cm <- tibble::deframe(res2_cm)
ranks_cm <- sort(ranks_cm, decreasing = TRUE)

# Run GSEA
fgseaRes_cm <- fgsea(pathways = pathways.hallmark, stats = ranks_cm)

# Tidy GSEA results
fgseaResTidy_cm <- fgseaRes_cm %>%
  as_tibble() %>%
  arrange(desc(NES))

# Filter significant pathways
fgseaResTidy1_cm <- filter(fgseaResTidy_cm, fgseaResTidy_cm$padj < 0.05)
fgseaResTidy1_cm <- fgseaResTidy1_cm %>% 
  mutate(variables = case_when(
    NES > 0 ~ "up", 
    NES < 0  ~ "down"))

fgseaResTidy1_cm
# Plot GSEA results

fgsea_plot2 <- ggplot(fgseaResTidy1_cm, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = variables)) +
  coord_flip() +
  scale_fill_manual(values = c("#4d9221", "#c51b7d")) +
  labs(
    x = "Pathway",
    y = "Normalized Enrichment Score",
    title = "FGSEA Pathway Enrichment CM-T vs CM"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = "../results/figures/gsea_hallmark_CM-T_vs_CM.png",
  plot = fgsea_plot2,
  width = 9,
  height = 7,
  dpi = 300
)
