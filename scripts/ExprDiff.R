# Cargar librerías
library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(plotly)
library(here)

# Definir rutas
main_path <- here()
metadata_path <- file.path(main_path, "data/metadata.tsv")
counts_path <- file.path(main_path, "results/quant/transcript_count_matrix.csv")

# Cargar datos
samples <- read.table(metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sample

counts <- read.csv(counts_path, header = TRUE, row.names = 1)
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

# Guardar matriz de distancias
write.table(sampleDistMatrix, "results/summary_qc/distance_matrix.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# Análisis de PCA
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pcaData_df <- data.frame(
  name = pcaData$name,
  x = pcaData$PC1,
  y = pcaData$PC2,
  color = pcaData$group
)
write.table(pcaData_df, "results/summary_qc/pca.txt", 
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

write.table(res_df_volcano, "results/summary_qc/volcano.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====== Matriz para heatmap ======

# Solo genes con resultados
counts_de <- counts[rownames(res_df), ]
counts_de_mx <- as.matrix(counts_de)

# Escalar por fila
counts_de_mx <- t(scale(t(counts_de_mx), center = TRUE, scale = TRUE))
counts_de_df <- rownames_to_column(as.data.frame(counts_de_mx), "gene")

write.table(counts_de_df, "results/summary_qc/abundance_de_matrix.txt", 
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
write.table(deg_df[1:50, ], "results/summary_qc/top50_DEGs.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")
