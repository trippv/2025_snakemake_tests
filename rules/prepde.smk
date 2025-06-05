rule prepDE:
    input:
        table="results/quant/samples_quant_table.txt"
    output:
        genes="results/quant/gene_count_matrix.csv",
        transcripts="results/quant/transcript_count_matrix.csv"
    conda:
        "../envs/prepde.yaml"
    shell:
        """
        prepDE.py -i {input.table} -g {output.genes} -t {output.transcripts}
        """
