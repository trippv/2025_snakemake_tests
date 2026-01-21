# regla para crear softlinks temporales con los nombres modificados de las lecturas crudas
rule create_soft_links:
    input:
        r1 = get_fastq1,
        r2 = get_fastq2
    output:
        r1_link = temp("results/fastp/raw/{sample}_R1_raw.{ext}"),
        r2_link = temp("results/fastp/raw/{sample}_R2_raw.{ext}")
    params:
        ext = get_extension,
        r1_abs = lambda wildcards: os.path.abspath(get_fastq1(wildcards)),
        r2_abs = lambda wildcards: os.path.abspath(get_fastq2(wildcards))
    shell:
        """
        ln -sf {params.r1_abs} {output.r1_link}
        ln -sf {params.r2_abs} {output.r2_link}
        """


rule fastp:
    input:
        r1 = lambda wildcards: f"results/fastp/raw/{wildcards.sample}_R1_raw.{SAMPLE_DICT[wildcards.sample]['extension']}",
        r2 = lambda wildcards: f"results/fastp/raw/{wildcards.sample}_R2_raw.{SAMPLE_DICT[wildcards.sample]['extension']}"
    output:
        r1_clean = "results/fastp/{sample}_R1.clean.fastq.gz",
        r2_clean = "results/fastp/{sample}_R2.clean.fastq.gz",
        json = "results/summary_qc/{sample}_fastp.json",
        html = "results/summary_qc/{sample}_fastp.html"
    conda:
        "../envs/fastp.yaml"
    threads: config["fastp_threads"]
    resources:
        mem_mb=config["fastp_mem"]
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1_clean} --out2 {output.r2_clean} \
              --json {output.json} --html {output.html} \
              --thread {threads} \
              --detect_adapter_for_pe
        """
