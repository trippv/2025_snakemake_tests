# reformatear el gtf si es necesario. Stringtie no acepta ciertos formatos de gtf
# En caso de que el formato sea incorrecto, se genera un nuevo gtf corregido

rule fix_gtf_format:
    input:
        gtf = config["gtf"]
    output:
        fixed_gtf = "results/genome/fixed_genome.gtf"
    log:
        "results/logs/fix_gtf_format.log"
    shell:
        """
        # Redirigimos la salida de los mensajes al archivo log definido en la regla
        exec > {log} 2>&1

        if grep -q 'transcript_id ""' {input.gtf}; then
            echo "[$(date)] Detectado formato NCBI incompatible con stringtie. El GTF corregido se guardará en {log}."
            awk '$3 != "gene" && $0 !~ /transcript_id ""/' {input.gtf} > {output.fixed_gtf}
        else
            echo "[$(date)] El GTF parece correcto. Procediendo sin cambios."
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
    threads: config["stringtie_threads"]
    resources:
        mem_mb=config["stringtie_mem"]
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
    threads: config["stringtie_merge_threads"]
    resources:
        mem_mb=config["stringtie_merge_mem"]
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
    threads: config["stringtie_threads"]
    resources:
        mem_mb=config["stringtie_mem"]
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

