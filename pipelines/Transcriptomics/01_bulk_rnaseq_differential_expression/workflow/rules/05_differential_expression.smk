rule differential_expression:
    input:
        counts="results/counts/gene_counts.txt",
        samples=config["samples"]

    output:
        results="results/differential_expression/deseq2_results.csv",
        normalized_counts="results/differential_expression/normalized_counts.csv",
        vst_counts="results/differential_expression/vst_counts.csv",
        summary="results/differential_expression/deseq2_summary.txt"

    conda:
        "../envs/deseq2.yaml"

    params:
        design=config["deseq2"]["design"],
        reference=config["deseq2"]["reference_condition"],
        contrast=config["deseq2"]["contrast_condition"],
        alpha=config["deseq2"]["alpha"],
        log2fc=config["deseq2"]["log2fc_threshold"]

    script:
        "../scripts/differential_expression.R"