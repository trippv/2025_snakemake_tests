library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(RColorBrewer)
library(here)


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
# NUEVA SECCIÓN PARA EL VOLCANO PLOT PERSONALIZADO
##########################################################

# Limpiar NA de los resultados
res_df_clean <- res_df[!is.na(res_df$padj), ]

# Identificar genes diferencialmente expresados (DEGs)
degs <- res_df_clean[res_df_clean$significant == "Significativo", ]

# Identificar genes no diferencialmente expresados (non-DEGs)
non_degs <- res_df_clean[res_df_clean$significant == "No significativo", ]

# 1. Seleccionar los 1000 genes más significativos (DEGs)
# Ordenar DEGs por p-valor ajustado de forma ascendente
degs <- degs[order(degs$padj, decreasing = FALSE), ]
# Tomar hasta 1000 genes, si hay menos, toma todos
top_1000_degs <- head(degs, 1000)

# 2. Seleccionar hasta 1000 genes no significativos (non-DEGs) al azar
# Si hay más de 1000, seleccionar una muestra aleatoria
if (nrow(non_degs) > 1000) {
  sampled_non_degs <- non_degs[sample(1:nrow(non_degs), 1000), ]
} else {
  # Si hay 1000 o menos, tomar todos
  sampled_non_degs <- non_degs
}

# 3. Combinar los dos conjuntos de genes
volcano_data_custom <- rbind(top_1000_degs, sampled_non_degs)

# 4. Preparar los datos para el volcano plot
volcano_data_custom$volcano_color <- "#b3aaaa"
volcano_data_custom$volcano_color[which(volcano_data_custom$significant == "Significativo" & volcano_data_custom$log2FoldChange > 1)] <- "red"
volcano_data_custom$volcano_color[which(volcano_data_custom$significant == "Significativo" & volcano_data_custom$log2FoldChange < -1)] <- "blue"
volcano_data_custom$padj = -log10(volcano_data_custom$padj)

# 5. Crear el volcano plot personalizado
#volcano_plot_custom <- ggplot(volcano_data_custom, aes(x = log2FoldChange, y = -log10(padj))) +
#  geom_point(aes(color = volcano_color), alpha = 0.6, size = 1.5) +
#  scale_color_manual(values = c("Upregulated" = "red", 
#                                "Downregulated" = "blue", 
#                                "Non-significant" = "#b3aaaa")) +
#  labs(
#    title = "Volcano Plot (Top DEGs + Sampled non-DEGs)",
#    x = expression(paste("log"[2], "Fold Change")),
#    y = expression(paste("-log"[10], "Adjusted p-value")),
#    color = "Gene Type"
#  ) +
#  theme_minimal() +
#  theme(legend.position = "bottom")

# Imprimir el plot para visualizarlo
#print(volcano_plot_custom)

# Opcionalmente, puedes guardar este plot en un archivo
# ggsave("results/summary_qc/volcano_custom.png", plot = volcano_plot_custom)

# Guardar tabla de volcano plot personalizado
write.table(volcano_data_custom, "results/summary_qc/volcano_custom.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# La sección de heatmap se mantiene, pero ahora se usarán los top DEGs
counts_de <- counts[rownames(top_1000_degs), ] 

counts_de_mx <- as.matrix(counts_de)

# Escalar manualmente por fila
counts_de_mx <- t(scale(t(counts_de_mx), center = TRUE, scale = TRUE))

counts_de_df <- as.data.frame(counts_de_mx)
counts_de_df <- rownames_to_column(counts_de_df, "gene")

# Exportar matriz de abundancia para heatmap
write.table(counts_de_df, "results/summary_qc/abundance_de_matrix_1000.txt" , row.names = FALSE, quote = FALSE, sep = "\t")