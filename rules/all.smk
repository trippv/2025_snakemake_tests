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
