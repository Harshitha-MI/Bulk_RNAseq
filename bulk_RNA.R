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
countdata <- read.csv("./RNAseq.csv", row.names = 1, header = TRUE)
View(countdata)

# Convert count data to a matrix
countdata <- as.matrix(round(countdata))
colnames(countdata)

#Define experimental conditions
condition <- factor(c(rep("SFM", 5), rep("SFM-T", 5), rep("CM", 5), rep("CM-T", 5)))

# Create a metadata table
metadata <- data.frame(row.names = colnames(countdata), condition)

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
head(res_sfm)
head(res_cm)

# Identify upregulated and downregulated genes
# Remove rows with NA values in log2FoldChange or padj columns
res_cleaned_sfm <- res_sfm[!is.na(res_sfm$log2FoldChange) & !is.na(res_sfm$padj), ]
res_cleaned_cm <- res_cm[!is.na(res_cm$log2FoldChange) & !is.na(res_cm$padj), ]

# Now apply the filtering conditions to sfm and sfm-t
res_filtered_sfm <- res_cleaned_sfm[abs(res_cleaned_sfm$log2FoldChange) >= 1 & res_cleaned_sfm$padj < 0.05, ]

upregulated_genes_sfm <- res_filtered_sfm[res_filtered_sfm$log2FoldChange > 0,]
downregulated_genes_sfm <- res_filtered_sfm[res_filtered_sfm$log2FoldChange < 0,]
# Number of upregulated and downregulated genes
cat("Upregulated genes SFM: ", nrow(upregulated_genes_sfm))
cat("Downregulated genes SFM: ", nrow(downregulated_genes_sfm))

# Now apply the filtering conditions to cm and cm-t
res_filtered_cm <- res_cleaned_cm[abs(res_cleaned_cm$log2FoldChange) >= 1 & res_cleaned_cm$padj < 0.05, ]

upregulated_genes_cm <- res_filtered_cm[res_filtered_cm$log2FoldChange > 0,]
downregulated_genes_cm <- res_filtered_cm[res_filtered_cm$log2FoldChange < 0,]
# Number of upregulated and downregulated genes
cat("Upregulated genes CM: ", nrow(upregulated_genes_cm))
cat("Downregulated genes CM: ", nrow(downregulated_genes_cm))

#Heatmap

# Extract top 50 differentially expressed genes in SFM and SFM-T samples
top_genes_sfm <- rownames(res_sfm)[order(res_sfm$padj)][1:50]

# Normalize counts for visualization
rlog_counts_sfm <- rlog(dds, blind = FALSE)

#Get SFM and SFM-T samples 
sfm_samples <- rownames(metadata)[metadata$condition %in% c("SFM", "SFM-T")]

# Extract normalized counts for top genes
top_gene_counts_sfm <- assay(rlog_counts_sfm)[top_genes_sfm, sfm_samples ]

# Generate heatmap
pheatmap(top_gene_counts_sfm,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE, fontsize = 6,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Heatmap of Top 50 Differentially Expressed Genes in SFM and SFM-T")

# Extract top 50 differentially expressed genes in CM and CM-T samples
top_genes_cm <- rownames(res_cm)[order(res_cm$padj)][1:50]

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
         main = "Heatmap of Top 50 Differentially Expressed Genes in CM and CM-T")

# Volcano Plot
EnhancedVolcano(res_sfm,
                lab = rownames(res_sfm), # Gene labels
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
EnhancedVolcano(res_cm,
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

# PCA Plot

# Perform PCA on rlog-transformed data
rlog_counts <- rlog(dds, blind = FALSE)
pca_data <- prcomp(t(assay(rlog_counts)))

# Prepare data for plotting
pca_df <- as.data.frame(pca_data$x)

# Generate PCA plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  scale_color_brewer(palette = "Set1")

#Sample relationship

#construct distance matrix
sampleDists <- dist(t(assay(rlog_counts)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata$condition, rownames(metadata), sep = "_")
colnames(sampleDistMatrix) <- paste(metadata$condition, rownames(metadata), sep = "_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)


# Prepare data for GSEA
res_sfm <- as.data.frame(res_sfm)
res_sfm <- res_sfm %>%
  tibble::rownames_to_column(var = "GeneSymbol")

res2 <- res_sfm %>% 
  dplyr::select(GeneSymbol, log2FoldChange)

ranks <- tibble::deframe(res2)

# Load Hallmark gene sets
pathways.hallmark <- gmtPathways("./h.all.v7.1.symbols.gmt")

# Run GSEA
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks, nperm = 1000)

# Tidy GSEA results
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Filter significant pathways
fgseaResTidy1 <- filter(fgseaResTidy, fgseaResTidy$padj < 0.05)
fgseaResTidy1 <- fgseaResTidy1 %>% 
  mutate(variables = case_when(
    NES > 0 ~ "up", 
    NES < 0  ~ "down"))

fgseaResTidy1
# Plot GSEA results
ggplot(fgseaResTidy1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = variables)) +
  coord_flip() + scale_fill_manual(values = c("#4d9221", "#c51b7d")) +
  labs(x = "Pathway", y = "Normalized Enrichment Score") +
  theme_bw(base_size = 12)
