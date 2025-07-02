# rules/rename_raw.smk

import os

def get_renamed_raw_output(wildcards):
    original_r1 = SAMPLE_DICT[wildcards.sample]["fastq1"] # Access the original path from SAMPLE_DICT
    # Determine if the original file is gzipped based on its extension
    is_gzipped = original_r1.endswith(".gz")

    # Construct the new filename, conditionally adding .gz
    base_name_r1 = f"results/raw_renamed/{wildcards.sample}_R1.raw.fastq"
    base_name_r2 = f"results/raw_renamed/{wildcards.sample}_R2.raw.fastq"

    output_r1 = base_name_r1 + (".gz" if is_gzipped else "")
    output_r2 = base_name_r2 + (".gz" if is_gzipped else "")

    return temp(output_r1), temp(output_r2)

rule rename_raw:
    input:
        r1_orig = get_fastq1,
        r2_orig = get_fastq2
    output:
        # Use the function to get dynamic output names
        # Snakemake automatically unpacks the tuple returned by the function
        # into the r1_renamed and r2_renamed outputs.
        r1_renamed=lambda wildcards: get_renamed_raw_output(wildcards)[0],
        r2_renamed=lambda wildcards: get_renamed_raw_output(wildcards)[1]
    shell:
        """
        mkdir -p $(dirname {output.r1_renamed})

        # Create symbolic links
        ln -s {input.r1_orig} {output.r1_renamed}
        ln -s {input.r2_orig} {output.r2_renamed}
        """