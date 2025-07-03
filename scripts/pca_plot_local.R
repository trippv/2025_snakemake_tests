
library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(plotly)
library(here)



here()

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
rownames(sampleDistMatrix) <-vst$sample
colnames(sampleDistMatrix) <-vst$sample
sampleDistMatrix <- as.data.frame(sampleDistMatrix)

sampleDistMatrix <- rownames_to_column(sampleDistMatrix, "sample")



# exportar matriz de distancia
#write.csv(sampleDists, snakemake@output[["gene_matrix"]], row.names = FALSE, quote = FALSE)
write.table(sampleDistMatrix, "results/summary_qc/distance_matrix.txt" , row.names = FALSE, quote = FALSE, sep = "\t")

# Crear objeto DESeq y transformar con VST
#dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ group)
#vsd <- vst(dds)

# Extraer datos de la PCA
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Graficar PCA con ggplot2 + plotly
gg <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  theme_minimal() 

ggplotly(gg)
```


# Expresion diferencial

# Realizar análisis de expresión diferencial
dds <- DESeq(dds)

group_names <- unique(samples$group)
# Contraste: cambia los nombres según tus grupos reales
# Por ejemplo: groupA = "Control", groupB = "Tratado"
res <- results(dds, contrast = c("group", group_names[1], group_names[2]))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Etiquetas para el volcano
res_df$significant <- "No significativo"
res_df$significant[which(res_df$padj < 0.01 & abs(res_df$log2FoldChange) >= 1)] <- "Significativo"



# Volcano Plot interactivo
gg_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), text = gene)) +
  geom_point(aes(color = significant), alpha = 0.7) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano plot", 
  x = "Log2 Fold Change", 
  y = "-Log10(padj)")

ggplotly(gg_volcano, tooltip = "text")
