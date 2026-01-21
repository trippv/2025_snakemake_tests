rule fix_gtf_format:
    input:
        gtf = config["gtf"]
    output:
        fixed_gtf = "results/genome/fixed_genome.gtf"
    shell:
        """
        # Verificamos si el archivo tiene el error 'transcript_id ""'
        if grep -q 'transcript_id ""' {input.gtf}; then
            echo "Detectado formato NCBI incompatible. Limpiando..."
            awk '$3 != "gene" && $0 !~ /transcript_id ""/' {input.gtf} > {output.fixed_gtf}
        else
            echo "El GTF parece correcto. Procediendo sin cambios."
            cp {input.gtf} {output.fixed_gtf}
        fi
        """


#Cuantificación inicial con ensamblaje
rule stringtie_quant:
    input:
        bam="results/align/{sample}.bam",
        gtf="results/genome/fixed_genome.gtf"
    output:
        gtf="results/quant/{sample}/{sample}.gtf"
    conda:
        "../envs/stringtie.yaml"
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
        gtf="results/genome/fixed_genome.gtf"
    output:
        merged="results/quant/stringtie_merged.gtf"
    conda:
        "../envs/stringtie.yaml"
    threads: 8
    params:
        memory=16
    shell:
        """
        stringtie --merge -p {threads} -G {input.gtf} -o {output.merged} {input.table}
        """

# Cuantificación final usando ensamblaje unificado
rule stringtie_final_quant:
    input:
        bam="results/align/{sample}.bam",
        merged_gtf="results/quant/stringtie_merged.gtf"
    output:
        gtf="results/quant_final/{sample}/{sample}.gtf"
    conda:
        "../envs/stringtie.yaml"
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

