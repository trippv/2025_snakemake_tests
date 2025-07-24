# Cargar la librería yaml
library(yaml)

# Leer el archivo volcano.txt
volcano <- read.delim("results/summary_qc/volcano.txt", header = TRUE, sep = "\t")

# Crear lista de puntos con coordenadas y color
volcano_points <- lapply(seq_len(nrow(volcano)), function(i) {
  list(
    x = round(volcano$logFC[i], 4),
    y = round(volcano$pval[i], 4),
    color = volcano$color[i]
  )
})
names(volcano_points) <- volcano$gene

# Crear el objeto completo
volcano_yaml <- list(
  id = "mqc_volcano",
  section_name = "Gráfico Volcano",
  description = "Este gráfico muestra los resultados del análisis de expresión diferencial",
  plot_type = "scatter",
  pconfig = list(
    id = "mqc_volcano_plot",
    title = "Volcano plot",
    xlab = "log2 Fold Change",
    ylab = "-log10(p-valor)"  # Cambia esto si estás usando -log10(p)
  ),
  data = volcano_points
)

# Escribir archivo YAML
writeLines(as.yaml(volcano_yaml), "results/summary_qc/volcano_mqc.yaml")
