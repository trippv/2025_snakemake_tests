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
        "fastq2": row["fastq2"],
        "extension": row["extension"]
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
# obtiene la extension del archivo
def get_extension(wildcards):
    return SAMPLE_DICT[wildcards.sample]["extension"]

# Regla principal para ejecutar todo el flujo
rule all:
    input:
        expand("results/quant/{sample}/{sample}.gtf", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}_gff.annotated.gtf", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}_gff.stats", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}_gff.loci", sample=SAMPLES),
        expand("results/summary_qc/gffcompare/{sample}_gff.tracking", sample=SAMPLES),
        # incluir archivos temporales de fastp
        expand("results/fastp/{sample}_R1.clean.fastq.gz", sample=SAMPLES),
        expand("results/fastp/{sample}_R2.clean.fastq.gz", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        expand("results/summary_qc/{sample}.hisat2.log", sample=SAMPLES),
        expand("results/summary_qc/samtools/{sample}.samtools.stats", sample=SAMPLES),
        # refseQC results
        expand("results/summary_qc/rseqc/{sample}.read_distribution.txt", sample=SAMPLES),
        expand("results/summary_qc/rseqc/{sample}.inner_distance.txt", sample=SAMPLES),
        expand("results/summary_qc/rseqc/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("results/summary_qc/rseqc/{sample}.seq.DupRate.xls", sample=SAMPLES),

        "results/abundance/samples_table.txt",
        "results/abundance/stringtie_merged.gtf",
        "results/abundance/samples_quant_table.txt",
        "results/abundance/gene_count_matrix.csv",
        "results/abundance/transcript_count_matrix.csv",
        #"results/summary_qc/pca_plot_mqc.html",
        #resultados para multiqc
        "results/summary_qc/abundance_heatmap_mqc.yaml",
        "results/summary_qc/pca_mqc.yaml",
        "results/summary_qc/volcano_mqc.yaml",
        "results/summary_qc/multiqc_report.html"


# Incluir reglas desde archivos
include: "rules/fastp.smk"
include: "rules/hisat2.smk"
include: "rules/samtools.smk"
include: "rules/rseqc.smk"
include: "rules/stringtie.smk"
include: "rules/gffcompare.smk"
include: "rules/prepde.smk"
include: "rules/metadata.smk"
include: "rules/expr_diff.smk"
include: "rules/mqc_reformat.smk"
include: "rules/multiqc.smk"
