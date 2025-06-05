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
        expand("results/summary_qc/{sample}.annotated.gtf", sample=SAMPLES),
        expand("results/summary_qc/{sample}.stats", sample=SAMPLES),
        expand("results/summary_qc/{sample}.loci", sample=SAMPLES),
        expand("results/summary_qc/{sample}.tracking", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.html", sample=SAMPLES),
        expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        expand("results/summary_qc/{sample}.hisat2.log", sample=SAMPLES),
        expand("results/summary_qc/{sample}.samtools.stats", sample=SAMPLES),
        "results/quant/samples_table.txt",
        "results/quant/stringtie_merged.gtf",
        "results/quant/samples_quant_table.txt",
        "results/quant/gene_count_matrix.csv",
        "results/quant/transcript_count_matrix.csv",
        "results/summary_qc/pca_plot_mqc.html"

# Regla para control de calidad y eliminación de adaptadores con fastp
rule fastp:
    input:
        r1 = get_fastq1,
        r2 = get_fastq2
    output:
        r1_clean="results/fastp/{sample}_R1.clean.fastq.gz",
        r2_clean="results/fastp/{sample}_R2.clean.fastq.gz",
        json="results/summary_qc/{sample}_fastp.json",
        html="results/summary_qc/{sample}_fastp.html"
    conda:
        "envs/fastp.yaml"
    threads: 4
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1_clean} --out2 {output.r2_clean} \
              --json {output.json} --html {output.html} \
              --thread {threads} \
              --detect_adapter_for_pe
        """

# Regla para extraer sitios de empalme para HISAT2
rule hisat2_extract_splicesites:
    input:
        gtf=config["gtf"] 
    output:
        "results/hisat2/splice_sites.txt"
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        hisat2_extract_splice_sites.py {input.gtf} > {output}
        """

# Regla para indexar el genoma con HISAT2
rule hisat2_index:
    input:
        fasta=config["genome_fasta"],
        splice_sites="results/hisat2/splice_sites.txt"
    output:
        directory("results/hisat2_index")
    conda:
        "envs/hisat2.yaml"
    threads: 8
    shell:
        """
        mkdir -p {output}
        hisat2-build -p {threads} --ss {input.splice_sites} {input.fasta} {output}/index
        """

# Regla para alinear lecturas con HISAT2
rule hisat2_align:
    input:
        r1="results/fastp/{sample}_R1.clean.fastq.gz",
        r2="results/fastp/{sample}_R2.clean.fastq.gz",
        index="results/hisat2_index",
        splice_sites="results/hisat2/splice_sites.txt"
    output:
        bam="results/align/{sample}.bam",
        bai="results/align/{sample}.bam.bai",
        log="results/summary_qc/{sample}.hisat2.log"
    conda:
        "envs/hisat2.yaml"
    threads: 8
    shell:
        """
        hisat2 -p {threads} \
               --dta \
               --known-splicesite-infile {input.splice_sites} \
               -x {input.index}/index \
               -1 {input.r1} -2 {input.r2} \
               2> {output.log} | \
               samtools view -bS - | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# Regla para generar estadísticas de alineamiento para MultiQC
rule samtools_stats:
    input:
        bam="results/align/{sample}.bam"
    output:
        stats="results/summary_qc/{sample}.samtools.stats"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        """

#Cuantificación inicial con ensamblaje
rule stringtie_quant:
    input:
        bam="results/align/{sample}.bam",
        gtf=config["gtf"]
    output:
        gtf="results/quant/{sample}/{sample}.gtf"
    conda:
        "envs/stringtie.yaml"
    threads: 8
    shell:
        """
        mkdir -p results/quant/{wildcards.sample}
        stringtie -p {threads} -G {input.gtf} -o {output.gtf} -l {wildcards.sample} {input.bam}
        """

# Crear tabla de ensamblajes
rule build_samples_table:
    input:
        gtfs=expand("results/quant/{sample}/{sample}.gtf", sample=SAMPLES)
    output:
        table="results/quant/samples_table.txt"
    run:
        with open(output.table, "w") as f:
            for sample in SAMPLES:
                f.write(f"results/quant/{sample}/{sample}.gtf\n")

# Unir ensamblajes
rule stringtie_merge:
    input:
        table="results/quant/samples_table.txt",
        gtf=config["gtf"]
    output:
        merged="results/quant/stringtie_merged.gtf"
    conda:
        "envs/stringtie.yaml"
    threads: 8
    shell:
        """
        stringtie --merge -p {threads} -G {input.gtf} -o {output.merged} {input.table}
        """

# Comparar con anotación de referencia
rule gffcompare:
    input:
        gtf="results/quant/{sample}/{sample}.gtf",
        ref=config["gtf"]

    output:
        stat="results/summary_qc/{sample}.stats",
        loci="results/summary_qc/{sample}.loci",
        annotated="results/summary_qc/{sample}.annotated.gtf",
        track="results/summary_qc/{sample}.tracking"
    conda:
        "envs/gffcompare.yaml"

    shell:
        """
        gffcompare -r {input.ref} -o results/summary_qc/{wildcards.sample} {input.gtf}
        """

# Cuantificación final usando ensamblaje unificado
rule stringtie_final_quant:
    input:
        bam="results/align/{sample}.bam",
        merged_gtf="results/quant/stringtie_merged.gtf"
    output:
        gtf="results/quant_final/{sample}/{sample}.gtf"
    conda:
        "envs/stringtie.yaml"
    threads: 8
    shell:
        """
        mkdir -p results/quant_final/{wildcards.sample}
        stringtie -e -B -p {threads} -G {input.merged_gtf} -o {output.gtf} {input.bam}
        """

# regla construir samples files para cuantificacion final
rule build_quant_sample_table:
    input:
        gtfs=expand("results/quant_final/{sample}/{sample}.gtf", sample=SAMPLES)
    output:
        table="results/quant/samples_quant_table.txt"
    run:
        with open(output.table, 'w') as f:
            for sample in SAMPLES:
                gtf_path = f"results/quant_final/{sample}/{sample}.gtf"
                f.write(f"{sample}\t{gtf_path}\n")



# Generar matrices con prepDE.py
rule prepDE:
    input:
        table="results/quant/samples_quant_table.txt"
    output:
        genes="results/quant/gene_count_matrix.csv",
        transcripts="results/quant/transcript_count_matrix.csv"
    conda:
        "envs/prepde.yaml"
    shell:
        """
        prepDE.py -i {input.table} -g {output.genes} -t {output.transcripts}
        """


rule generate_metadata:
    input:
        sample_file = "config/samples.tsv"
    output:
        metadata = "data/metadata.tsv"
    run:
        import os
        import csv

        def load_samples(sample_file):
            samples = []
            with open(sample_file, newline='') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    if row["include"].strip() == "1":
                        samples.append({
                            "sample_id": row["sample_id"],
                            "group": row["group"]
                        })
            return samples

        samples = load_samples(input.sample_file)
        os.makedirs(os.path.dirname(output.metadata), exist_ok=True)
        with open(output.metadata, "w") as f:
            f.write("sample\tgroup\n")
            for sample in samples:
                f.write(f"{sample['sample_id']}\t{sample['group']}\n")




rule pca_html:
    input:
        counts = "results/quant/gene_count_matrix.csv",
        metadata = "data/metadata.tsv",
        rmd = "scripts/pca_plot.Rmd"
    output:
        html = "results/summary_qc/pca_plot_mqc.html"
    conda:
        "envs/r.yaml"
    shell:
        """
        MAIN=$(pwd)
        Rscript -e "rmarkdown::render('{input.rmd}', params = list(main='$MAIN', counts='{input.counts}', metadata='{input.metadata}'))"
        mv scripts/pca_plot.html {output.html}

        """
