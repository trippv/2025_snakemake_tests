###############################################################################
# MultiQC – reporte global de control de calidad
###############################################################################
rule multiqc:
    input:
        # Todos los archivos que MultiQC debe escanear
        sample_names = "data/samples_mqc.txt",
        fastp_json      = expand("results/summary_qc/{sample}_fastp.json", sample=SAMPLES),
        fastp_html      = expand("results/summary_qc/{sample}_fastp.html", sample=SAMPLES),
        hisat2_logs     = expand("results/summary_qc/{sample}.hisat2.log",   sample=SAMPLES),
        samtools_stats  = expand("results/summary_qc/samtools/{sample}.samtools.stats", sample=SAMPLES),
        gff_stats       = expand("results/summary_qc/gffcompare/{sample}_gff.stats",        sample=SAMPLES),
        gff_tracking    = expand("results/summary_qc/gffcompare/{sample}_gff.tracking",     sample=SAMPLES),
        gff_loci        = expand("results/summary_qc/gffcompare/{sample}_gff.loci",         sample=SAMPLES),
        heatmap = "results/summary_qc/abundance_heatmap_mqc.yaml",
        pca = "results/summary_qc/pca_mqc.yaml",
        volcano = "results/summary_qc/volcano_mqc.yaml",
    output:
        report = "results/summary_qc/multiqc_report.html"
    params:
        comment = config["multiqc_comment"]
    conda:
        "../envs/multiqc.yaml"
    threads: 4
    shell:
        """
        multiqc results/summary_qc -o results/summary_qc -f --config config/multiqc_config.yaml --template mi_multiqc --comment "{params.comment}"
        """
