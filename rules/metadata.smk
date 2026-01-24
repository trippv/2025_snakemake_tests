# genera tabla con metadatos a partir del archivo de muestras para el análisis de DESeq2 y otros análisis downstream
rule generate_metadata:
    input:
        sample_file = config["samples_file"]
    output:
        metadata = "data/metadata.tsv"
    shell:
        """
        mkdir -p $(dirname {output.metadata})
        awk -F '\t' 'NR==1 {{print "sample\tgroup"}} NR>1 && $6==1 {{print $1"\t"$2}}' {input.sample_file} > {output.metadata}
        """



rule generate_sample_names_for_multiqc:
    input:
        sample_file = config["samples_file"]
    output:
        names = "data/samples_mqc.txt"
    run:
        import os
        import csv

        def get_basename(path):
            return os.path.splitext(os.path.basename(path))[0]

        os.makedirs(os.path.dirname(output.names), exist_ok=True)

        with open(input.sample_file, newline='') as f_in, open(output.names, "w") as f_out:
            reader = csv.DictReader(f_in, delimiter="\t")
            for row in reader:
                if row["include"].strip() == "1":
                    sample_id = row["sample_id"]
                    raw1 = get_basename(row["fastq1"])
                    raw2 = get_basename(row["fastq2"])
                    f_out.write(f"{raw1}\t{sample_id}\n")
                    f_out.write(f"{raw2}\t{sample_id}\n")


