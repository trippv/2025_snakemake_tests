rule prepDE:
    input:
        table="results/abundance/samples_quant_table.txt"
    output:
        genes="results/abundance/gene_count_matrix.csv",
        transcripts="results/abundance/transcript_count_matrix.csv"
    conda:
        "../envs/prepde.yaml"
    shell:
        """
        prepDE.py -i {input.table} -g {output.genes} -t {output.transcripts}
        """
