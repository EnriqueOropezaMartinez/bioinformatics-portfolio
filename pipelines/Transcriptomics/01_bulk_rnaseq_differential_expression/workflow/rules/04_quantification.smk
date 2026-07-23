rule featurecounts:
    input:
        annotation=config["reference"]["annotation"],
        bams=expand(
            "results/alignment/{sample}/{sample}.sorted.bam",
            sample=SAMPLES
        )

    output:
        counts="results/counts/gene_counts.txt",
        summary="results/counts/gene_counts.txt.summary"

    threads:
        config["threads"]["featurecounts"]

    conda:
        "../envs/quantification.yaml"

    params:
        feature_type=config["featurecounts"].get("feature_type", "exon"),
        attribute=config["featurecounts"].get("attribute", "gene_id"),
        strand=config["featurecounts"].get("strand", 0)

    shell:
        """
        mkdir -p results/counts

        featureCounts \
            -T {threads} \
            -p \
            --countReadPairs \
            -s {params.strand} \
            -t {params.feature_type} \
            -g {params.attribute} \
            -a {input.annotation} \
            -o {output.counts} \
            {input.bams}
        """