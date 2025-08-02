
library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(RColorBrewer)


#here()

# Definir rutas absolutas basadas en el parámetro 'main'
main_path <- here()
metadata_path <- file.path(main_path, "data/metadata.tsv")
counts_path <- file.path(main_path, "results/quant/transcript_count_matrix.csv")

# Cargar metadatos y matriz de cuentas
samples <- read.table(metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sample

counts <- read.csv(counts_path, header = TRUE, row.names = 1)
counts <- counts[, samples$sample, drop = FALSE]  # Asegura que el orden y las muestras coincidan


# filtrar genes con poca repreesntacion
keep <- rowSums(counts >= 10) >= 2
counts_filtered <- counts[keep, ]
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = samples, design = ~ group)
vsd <- varianceStabilizingTransformation(dds)


### matriz de distancia entre muestras
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <-vsd$sample
colnames(sampleDistMatrix) <-vsd$sample
sampleDistMatrix <- as.data.frame(sampleDistMatrix)

sampleDistMatrix <- rownames_to_column(sampleDistMatrix, "sample")



# exportar matriz de distancia
#write.csv(sampleDists, snakemake@output[["gene_matrix"]], row.names = FALSE, quote = FALSE)
write.table(sampleDistMatrix, "results/summary_qc/distance_matrix.txt" , row.names = FALSE, quote = FALSE, sep = "\t")

# Crear objeto DESeq y transformar con VST
#dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ group)
#vsd <- vst(dds)


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


write.table(pcaData_df, "results/summary_qc/pca.txt" , row.names = FALSE, quote = FALSE, sep = "\t")


# Expresion diferencial

# Realizar análisis de expresión diferencial
dds <- DESeq(dds)

# construir tabla con los resultados de expresion difenrcial para un contraste específico
group_names <- unique(samples$group)
# Contraste: cambia los nombres según tus grupos reales
# Por ejemplo: groupA = "Control", groupB = "Tratado"
res <- results(dds, contrast = c("group", group_names[1], group_names[2]))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Etiquetas para el volcano
res_df$significant <- "No significativo"
res_df$significant[which(res_df$padj < 0.01 & abs(res_df$log2FoldChange) >= 1)] <- "Significativo"

# Volcano plot con el numero total del número de genes
res_df_volcano <- data.frame(gene = rownames(res_df), 
      logFC = res_df$log2FoldChange, 
      pval = -log(res_df$padj),
      significant = res_df$significant)

res_df_volcano$color <- ifelse(res_df_volcano$significant != "Significativo", '#b3aaaa', 
        ifelse(res_df_volcano$logFC > 1 , 'red', 'blue'))

# guardar tabla para volcano plot
write.table(res_df_volcano, "results/summary_qc/volcano.txt" , row.names = FALSE, quote = FALSE, sep = "\t")



##########################################################
# Filtrar los primeros 150 genes con menor pvalue
res_df_150 <- res_df[!is.na(res_df$padj), ]  # Eliminar filas con NA en 'padj'
res_df_150 <- res_df_150[res_df_150$significant == "Significativo", ]  # Filtrar solo los significativos
res_df_150 <- res_df_150[order(res_df_150$padj, decreasing = FALSE), ]
# Tomar los primeros 150 genes
res_df_150 <- head(res_df_150, 1500)


# Con estos resultados contruir el volcano plot

res_df_volcano_150 <- data.frame(gene = rownames(res_df_150), 
      logFC = res_df_150$log2FoldChange, 
      pval = -log(res_df_150$padj),
      significant = res_df_150$significant)

res_df_volcano_150$color <- ifelse(res_df_volcano_150$significant != "Significativo", '#b3aaaa', 
        ifelse(res_df_volcano$logFC > 1 , 'red', 'blue'))

# guardar tabla para volcano plot
write.table(res_df_volcano_150, "results/summary_qc/volcano_150.txt" , row.names = FALSE, quote = FALSE, sep = "\t")


## Heatmap con los valores de abundancia de los genes expresados diferencialmente 
counts_de <- counts[rownames(res_df_150), ] 

counts_de_mx <- as.matrix(counts_de)

# Escalar manualmente por fila
counts_de_mx <- t(scale(t(counts_de_mx), center = TRUE, scale = TRUE))

counts_de_df <- as.data.frame(counts_de_mx)
counts_de_df <- rownames_to_column(counts_de_df, "gene")

# Exportar matriz de abundancia para heatmap
write.table(counts_de_df, "results/summary_qc/abundance_de_matrix_150.txt" , row.names = FALSE, quote = FALSE, sep = "\t")


# Exportar lista de genes expresados diferencialmente (top 50)
deg_df <- data.frame(gene = res_df$gene,
                    baseMean = res_df$baseMean,
                    log2FC = res_df$log2FoldChange,
                    pval = res_df$pvalue,
                    adjustPval = res_df$padj)


deg_df <- deg_df[order(deg_df$adjustPval, decreasing = FALSE), ]


