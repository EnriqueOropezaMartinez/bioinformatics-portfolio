rule visualization:
    input:
        results="results/differential_expression/deseq2_results.csv",
        vst_counts="results/differential_expression/vst_counts.csv",
        samples=config["samples"]

    output:
        pca="results/figures/pca_plot.pdf",
        volcano="results/figures/volcano_plot.pdf",
        heatmap="results/figures/top_genes_heatmap.pdf",
        ma_plot="results/figures/ma_plot.pdf"

    conda:
        "../envs/visualization.yaml"

    params:
        design=config["deseq2"]["design"],
        alpha=config["deseq2"]["alpha"],
        log2fc=config["deseq2"]["log2fc_threshold"],
        top_genes=config["visualization"].get("top_genes", 50)

    script:
        "../scripts/visualization.R"