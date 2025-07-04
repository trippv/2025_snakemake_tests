# Cargar librería
library(yaml)

# Leer el archivo de abundancias
abund <- read.delim("results/summary_qc/abundance_de_matrix.txt", header = TRUE, sep = "\t")

# Convertir los datos a una lista anidada
abund_data <- lapply(seq_len(nrow(abund)), function(i) {
  gene_values <- abund[i, -1]
  as.list(setNames(as.numeric(gene_values), colnames(gene_values)))
})
names(abund_data) <- abund$gene

# Construir objeto YAML
abund_yaml <- list(
  id = "mqc_abundance",
  section_name = "Mapa de calor de abundancias",
  description = "Este heatmap muestra los valores de abundancia normalizados por gen",
  plot_type = "heatmap",
  pconfig = list(
    id = "mqc_abundance_heatmap",
    title = "Abundancias normalizadas",
    xlab = "Muestras",
    ylab = "Genes"
  ),
  data = abund_data
)

# Guardar el archivo YAML
writeLines(as.yaml(abund_yaml), "results/summary_qc/abundance_heatmap_mqc.yaml")
