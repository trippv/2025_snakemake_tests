rule pca_html:
    input:
        counts = "results/quant/gene_count_matrix.csv",
        metadata = "data/metadata.tsv",
        rmd = "scripts/pca_plot.Rmd"
    output:
        html = "results/summary_qc/pca_plot_mqc.html"
    conda:
        "../envs/r.yaml"
    shell:
        """
        MAIN=$(pwd)
        Rscript -e "rmarkdown::render('{input.rmd}', params = list(main='$MAIN', counts='{input.counts}', metadata='{input.metadata}'))"
        mv scripts/pca_plot.html {output.html}
        """
