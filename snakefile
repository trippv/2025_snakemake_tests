import pandas as pd

# Cargar archivo de configuración
configfile: "config/config.yaml"

# Leer tabla de muestras
SAMPLES_TABLE = pd.read_csv(config["samples_file"], sep="\t")

# Filtrar muestras incluidas
INCLUDED_SAMPLES = SAMPLES_TABLE[SAMPLES_TABLE["include"] == 1]

# Diccionario con datos de muestras
SAMPLE_DICT = {
    row["sample_id"]: {
        "group": row["group"],
        "fastq1": row["fastq1"],
        "fastq2": row["fastq2"]
    }
    for _, row in INCLUDED_SAMPLES.iterrows()
}

# Lista de IDs de muestras incluidas
SAMPLES = list(SAMPLE_DICT.keys())

# Funciones para obtener rutas de fastq
def get_fastq1(wildcards):
    return SAMPLE_DICT[wildcards.sample]["fastq1"]

def get_fastq2(wildcards):
    return SAMPLE_DICT[wildcards.sample]["fastq2"]

# Regla principal para ejecutar todo el flujo
rule all:
    input:
        expand("results/quant/{sample}/{sample}.gtf", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}.annotated.gtf", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}.stats", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}.loci", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}.tracking", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        expand("results/summary_qc/{sample}.hisat2.log", sample=SAMPLES),
        expand("results/summary_qc/samtools/{sample}.samtools.stats", sample=SAMPLES),
        "results/quant/samples_table.txt",
        "results/quant/stringtie_merged.gtf",
        "results/quant/samples_quant_table.txt",
        "results/quant/gene_count_matrix.csv",
        "results/quant/transcript_count_matrix.csv",
        "results/summary_qc/pca_plot_mqc.html",
        "results/summary_qc/multiqc_report.html"


# Incluir reglas desde archivos
include: "rules/fastp.smk"
include: "rules/hisat2.smk"
include: "rules/samtools.smk"
include: "rules/stringtie.smk"
include: "rules/gffcompare.smk"
include: "rules/prepde.smk"
include: "rules/metadata.smk"
include: "rules/pca.smk"
include: "rules/multiqc.smk"
