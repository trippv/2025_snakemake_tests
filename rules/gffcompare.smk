rule gffcompare:
    input:
        gtf="results/quant/{sample}/{sample}.gtf",
        ref=config["gtf_quant"]
    output:
        stat="results/summary_qc/gffcompare/{sample}_gff.stats",
        loci="results/summary_qc/gffcompare/{sample}_gff.loci",
        annotated="results/summary_qc/gffcompare/{sample}_gff.annotated.gtf",
        track="results/summary_qc/gffcompare/{sample}_gff.tracking"
    conda:
        "../envs/gffcompare.yaml"
    shell:
        """
        
        gffcompare -r {input.ref} -o results/summary_qc/gffcompare/{wildcards.sample}_gff {input.gtf}
        """
