rule functional_enrichment:
    input:
        results="results/differential_expression/deseq2_results.csv"

    output:
        go_up="results/enrichment/go_upregulated.csv",
        go_down="results/enrichment/go_downregulated.csv",
        kegg="results/enrichment/kegg_differentially_expressed.csv",
        go_up_plot="results/figures/go_upregulated_dotplot.pdf",
        go_down_plot="results/figures/go_downregulated_dotplot.pdf",
        kegg_plot="results/figures/kegg_dotplot.pdf",
        summary="results/enrichment/enrichment_summary.txt"

    conda:
        "../envs/enrichment.yaml"

    params:
        alpha=config["deseq2"]["alpha"],
        log2fc=config["deseq2"]["log2fc_threshold"],
        organism=config["enrichment"].get("organism", "hsa"),
        id_type=config["enrichment"].get("id_type", "ENSEMBL"),
        show_categories=config["enrichment"].get("show_categories", 15)

    script:
        "../scripts/functional_enrichment.R"