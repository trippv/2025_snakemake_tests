
library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(plotly)
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

# convert
pcaData_df <- data.frame(name = pcaData$name, 
x = pcaData$PC1, 
y = pcaData$PC2, 
color = pcaData$group) # falta modificar este parametro para asignar un color a cada grupo

write.table(pcaData_df, "results/summary_qc/pca.txt" , row.names = FALSE, quote = FALSE, sep = "\t")


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

res_df_volcano <- data.frame(gene = rownames(res_df), 
      logFC = res_df$log2FoldChange, 
      pval = -log(res_df$padj),
      significant = res_df$significant)


res_df_volcano$color <- ifelse(res_df$significant != "Significativo", '#b3aaaa', 
        ifelse(res_df$log2FoldChange > 1 , 'red', 'blue'))

# guardar tabla para volcano plot
write.table(res_df_volcano, "results/summary_qc/volcano.txt" , row.names = FALSE, quote = FALSE, sep = "\t")
