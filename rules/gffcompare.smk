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
        "../envs/gffcompare.yaml"
    shell:
        """
        
        gffcompare -r {input.ref} -o results/summary_qc/{wildcards.sample} {input.gtf}
        """
