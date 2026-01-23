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

# Leer umbrales desde el config de Snakemake
# Se usa un valor por defecto (default) por si no están en el yaml
fdr_limit <- snakemake@config[["diff_exp"]][["fdr_threshold"]]
lfc_limit <- snakemake@config[["diff_exp"]][["lfc_threshold"]]

# ==== Importante: A partir de aqui solo se ilustra los resultados de los dos primeros grupos. Para mas grupos
# se debera realizar el análisis de forma manual

dds <- DESeq(dds)
group_names <- unique(samples$group)
res <- results(dds, contrast = c("group", group_names[1], group_names[2]))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)


# Etiquetar genes significativos usando las variables del config
res_df$significant <- "No significativo"
res_df$significant[which(res_df$padj < fdr_limit & abs(res_df$log2FoldChange) >= lfc_limit)] <- "Significativo"

#retener solo significativos
res_df_significant <- res_df[res_df$significant == "Significativo", ]

# si res_df__significant tiene mas de 200 filas, quedarnos con las primeras 200 ordenadas por pvalor
if (nrow(res_df_significant) > 200) {
  res_df_significant <- res_df_significant[order(res_df_significant$padj), ][1:200, ]
}

# Preparar tabla para volcano plot
res_df_volcano <- data.frame(
  gene = rownames(res_df),
  logFC = res_df$log2FoldChange,
  pval = -log(res_df$padj),
  significant = res_df$significant
)

# El color del volcano ahora respeta el lfc_limit del config
res_df_volcano$color <- ifelse(
  res_df$significant != "Significativo", "#b3aaaa",
  ifelse(res_df$log2FoldChange >= lfc_limit, "red", "blue")
)

write.table(res_df_volcano, snakemake@output[["volcano"]], 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====== Matriz para heatmap ======

# Solo genes con resultados
counts_de <- counts[rownames(res_df_significant), ]
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
