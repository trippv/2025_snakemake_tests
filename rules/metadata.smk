rule generate_metadata:
    input:
        sample_file = "config/samples.tsv"
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
