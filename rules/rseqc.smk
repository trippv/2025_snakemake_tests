
rule rseqc_gtf_to_bed:
    input:
        gtf=config["gtf"]
    output:
        bed="results/summary_qc/rseqc/annotation.bed"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        perl scripts/gtf2bed.pl < {input.gtf} > {output.bed}
        """

# Read distribution para evaluar donde cae las lecturas en el genoma
rule rseqc_read_distribution:
    input:
        bam="results/align/{sample}.bam",
        bed="results/summary_qc/rseqc/annotation.bed"
    output:
        txt="results/summary_qc/rseqc/{sample}.read_distribution.txt"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} > {output.txt}

        """

# Inner distance para evaluar la distancia entre pares de lecturas
rule rseqc_inner_distance:
    input:
        bam="results/align/{sample}.bam",
        bed="results/summary_qc/rseqc/annotation.bed"
    output:
        txt="results/summary_qc/rseqc/{sample}.inner_distance.txt",
        pdf="results/summary_qc/rseqc/{sample}.inner_distance_plot.pdf"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        inner_distance.py \
            -i {input.bam} \
            -r {input.bed} \
            -o results/summary_qc/rseqc/{wildcards.sample}
        """


# Gene body coverage para evaluar la cobertura a lo largo del cuerpo del gen
rule rseqc_gene_body_coverage:
    input:
        bam="results/align/{sample}.bam",
        bed="results/summary_qc/rseqc/annotation.bed"
    output:
        txt="results/summary_qc/rseqc/{sample}.geneBodyCoverage.txt",
        pdf="results/summary_qc/rseqc/{sample}.geneBodyCoverage.curves.pdf"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        geneBody_coverage.py \
            -i {input.bam} \
            -r {input.bed} \
            -o results/summary_qc/rseqc/{wildcards.sample}
        """


# Read duplication para evaluar la duplicacion de lecturas
rule rseqc_read_duplication:
    input:
        bam="results/align/{sample}.bam"
    output:
        pos="results/summary_qc/rseqc/{sample}.pos.DupRate.xls",
        seq="results/summary_qc/rseqc/{sample}.seq.DupRate.xls"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_duplication.py \
            -i {input.bam} \
            -o results/summary_qc/rseqc/{wildcards.sample}
        """
