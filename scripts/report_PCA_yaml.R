# Cargar librería necesaria
library(yaml)

# Leer el archivo PCA
pca <- read.delim("results/summary_qc/pca_df.txt", header = TRUE, sep = "\t")

# Crear la lista de puntos
pca_points <- lapply(seq_len(nrow(pca)), function(i) {
  list(
    x = round(pca$x[i], 4),
    y = round(pca$y[i], 4),
    color = pca$color[i]
  )
})
names(pca_points) <- pca$sample

# Crear estructura YAML para MultiQC
pca_yaml <- list(
  id = "mqc_pca",
  section_name = "Análisis PCA",
  description = "Este gráfico muestra los resultados del PCA con colores por grupo",
  plot_type = "scatter",
  pconfig = list(
    id = "mqc_pca_scatter_plot",
    title = "PCA",
    ylab = "Componente 2",
    xlab = "Componente 1"
  ),
  data = pca_points
)

# Escribir el archivo YAML
writeLines(as.yaml(pca_yaml), "results/summary_qc/pac_mqc.yaml")
