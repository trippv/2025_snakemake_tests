# Regla para extraer sitios de empalme para HISAT2
rule hisat2_extract_splicesites:
    input:
        gtf=config["gtf"] 
    output:
        "results/hisat2/splice_sites.txt"
    conda:
        "../envs/hisat2.yaml"
    shell:
        """
        hisat2_extract_splice_sites.py {input.gtf} > {output}
        """

# Regla para indexar el genoma con HISAT2
rule hisat2_index:
    input:
        fasta=config["genome_fasta"],
        splice_sites="results/hisat2/splice_sites.txt"
    output:
        directory("results/hisat2_index")
    conda:
        "../envs/hisat2.yaml"
    threads: 8
    shell:
        """
        mkdir -p {output}
        hisat2-build -p {threads} --ss {input.splice_sites} {input.fasta} {output}/index
        """

# Regla para alinear lecturas con HISAT2
rule hisat2_align:
    input:
        r1="results/fastp/{sample}_R1.clean.fastq.gz",
        r2="results/fastp/{sample}_R2.clean.fastq.gz",
        index="results/hisat2_index",
        splice_sites="results/hisat2/splice_sites.txt"
    output:
        bam="results/align/{sample}.bam",
        bai="results/align/{sample}.bam.bai",
        log="results/summary_qc/{sample}.hisat2.log"
    conda:
        "../envs/hisat2.yaml"
    threads: 8
    shell:
        """
        hisat2 -p {threads} \
               --dta \
               --known-splicesite-infile {input.splice_sites} \
               -x {input.index}/index \
               -1 {input.r1} -2 {input.r2} \
               2> {output.log} | \
               samtools view -bS - | samtools sort -o {output.bam}
        samtools index {output.bam}
        """