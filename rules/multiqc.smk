###############################################################################
# MultiQC – reporte global de control de calidad
###############################################################################
rule multiqc:
    input:
        # Todos los archivos que MultiQC debe escanear
        fastp_json      = expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        fastp_html      = expand("results/summary_qc/{sample}_fastp.html", sample=SAMPLES),
        hisat2_logs     = expand("results/summary_qc/{sample}.hisat2.log",   sample=SAMPLES),
        samtools_stats  = expand("results/summary_qc/samtools/{sample}.samtools.stats", sample=SAMPLES),
        gff_stats       = expand("results/summary_qc/gffcompare/{sample}.stats",        sample=SAMPLES),
        gff_tracking    = expand("results/summary_qc/gffcompare/{sample}.tracking",     sample=SAMPLES),
        gff_loci        = expand("results/summary_qc/gffcompare/{sample}.loci",         sample=SAMPLES),
        pca_html        = "results/summary_qc/pca_plot_mqc.html"
    output:
        report = "results/summary_qc/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    threads: 4
    shell:
        """
        multiqc results/summary_qc -o results/summary_qc -f
        """
