rule expr_diff:
    input:
        counts = "results/quant/transcript_count_matrix.csv",
        metadata = "data/metadata.tsv"
    output:
        dist_matrix="results/summary_qc/distance_matrix.txt",
        pca="results/summary_qc/pca.txt",
        abundance="results/summary_qc/abundance.txt",
        volcano="results/summary_qc/volcano.txt",
        deg="results/summary_qc/deg.txt"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/ExprDiff.R"
