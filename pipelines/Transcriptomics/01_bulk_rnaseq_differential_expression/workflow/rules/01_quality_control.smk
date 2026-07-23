rule fastqc_raw:
    input:
        fastq="data/raw/{sample}_{read}.fastq.gz"

    output:
        html="results/qc/raw/fastqc/{sample}_{read}_fastqc.html",
        zip="results/qc/raw/fastqc/{sample}_{read}_fastqc.zip"

    threads:
        config["threads"]["fastqc"]

    conda:
        "../envs/qc.yaml"

    shell:
        """
        mkdir -p results/qc/raw/fastqc

        fastqc \
            --threads {threads} \
            --outdir results/qc/raw/fastqc \
            {input.fastq}
        """


rule multiqc_raw:
    input:
        html=expand(
            "results/qc/raw/fastqc/{sample}_{read}_fastqc.html",
            sample=SAMPLES,
            read=READS
        ),
        zip=expand(
            "results/qc/raw/fastqc/{sample}_{read}_fastqc.zip",
            sample=SAMPLES,
            read=READS
        )

    output:
        report="results/qc/raw/multiqc_report.html"

    conda:
        "../envs/qc.yaml"

    shell:
        """
        mkdir -p results/qc/raw

        multiqc \
            results/qc/raw/fastqc \
            --outdir results/qc/raw \
            --force
        """