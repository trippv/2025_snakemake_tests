#!/usr/bin/env Rscript

library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(RColorBrewer)

# Cargar datos
samples <- read.table(snakemake@input[["metadata"]], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sample

counts <- read.csv(snakemake@input[["counts"]],, header = TRUE, row.names = 1)
counts <- counts[, samples$sample, drop = FALSE]  # Asegura orden y correspondencia

# Filtrar genes con poca representación
keep <- rowSums(counts >= 10) >= 2
counts_filtered <- counts[keep, ]

# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = samples, design = ~ group)
vsd <- varianceStabilizingTransformation(dds)

# ====== Análisis exploratorio ======

# Matriz de distancias
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samples$sample
colnames(sampleDistMatrix) <- samples$sample
sampleDistMatrix <- rownames_to_column(as.data.frame(sampleDistMatrix), "sample")


#write.csv(counts, snakemake@output[["gene_matrix"]], row.names = FALSE, quote = FALSE)
# Guardar matriz de distancias
write.table(sampleDistMatrix, snakemake@output[["dist_matrix"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ========== Paleta de colores por grupo ==========
groups <- unique(samples$group)
n_groups <- length(groups)
colors <- setNames(brewer.pal(min(max(n_groups, 3), 8), "Set1")[1:n_groups], groups)


# Análisis de PCA
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pcaData$color <- colors[as.character(pcaData$group)]

pcaData_df <- data.frame(
  sample = pcaData$name,
  PC1 = pcaData$PC1,
  PC2 = pcaData$PC2,
  group = pcaData$group,
  color = pcaData$color
)

write.table(pcaData_df, snakemake@output[["pca"]],
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====== Expresión diferencial ======

# ==== Importante: A partir de aqui solo se ilustra los resultados de los dos primeros grupos. Para mas grupos
# se debera realizar el análisis de forma manual

dds <- DESeq(dds)
group_names <- unique(samples$group)
res <- results(dds, contrast = c("group", group_names[1], group_names[2]))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Etiquetar genes significativos
res_df$significant <- "No significativo"
res_df$significant[which(res_df$padj < 0.01 & abs(res_df$log2FoldChange) >= 1)] <- "Significativo"

# Preparar tabla para volcano plot
res_df_volcano <- data.frame(
  gene = rownames(res_df),
  logFC = res_df$log2FoldChange,
  pval = -log(res_df$padj),
  significant = res_df$significant
)

res_df_volcano$color <- ifelse(
  res_df$significant != "Significativo", "#b3aaaa",
  ifelse(res_df$log2FoldChange > 1, "red", "blue")
)

write.table(res_df_volcano, snakemake@output[["volcano"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====== Matriz para heatmap ======

# Solo genes con resultados
counts_de <- counts[rownames(res_df), ]
counts_de_mx <- as.matrix(counts_de)

# Escalar por fila
counts_de_mx <- t(scale(t(counts_de_mx), center = TRUE, scale = TRUE))
counts_de_df <- rownames_to_column(as.data.frame(counts_de_mx), "gene")

write.table(counts_de_df, snakemake@output[["abundance"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====== Top DEGs ======

deg_df <- data.frame(
  gene = res_df$gene,
  baseMean = res_df$baseMean,
  log2FC = res_df$log2FoldChange,
  pval = res_df$pvalue,
  adjustPval = res_df$padj
)

deg_df <- deg_df[order(deg_df$adjustPval), ]
write.table(deg_df[1:50, ], snakemake@output[["deg"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
