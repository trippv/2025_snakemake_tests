rule samtools_stats:
    input:
        bam="results/align/{sample}.bam"
    output:
        stats="results/summary_qc/{sample}.samtools.stats"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        """
