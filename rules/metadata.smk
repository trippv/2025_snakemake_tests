rule generate_metadata:
    input:
        sample_file = config["samples_file"]
    output:
        metadata = "data/metadata.tsv"
    run:
        import os
        import csv

        def load_samples(sample_file):
            samples = []
            with open(sample_file, newline='') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    if row["include"].strip() == "1":
                        samples.append({
                            "sample_id": row["sample_id"],
                            "group": row["group"]
                        })
            return samples

        samples = load_samples(input.sample_file)
        os.makedirs(os.path.dirname(output.metadata), exist_ok=True)
        with open(output.metadata, "w") as f:
            f.write("sample\tgroup\n")
            for sample in samples:
                f.write(f"{sample['sample_id']}\t{sample['group']}\n")


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


