rule trim_galore:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz"

    output:
        r1="results/trimmed/{sample}_R1_val_1.fq.gz",
        r2="results/trimmed/{sample}_R2_val_2.fq.gz"

    threads:
        config["threads"]["trimming"]

    conda:
        "../envs/trimming.yaml"

    shell:
        """
        mkdir -p results/trimmed

        trim_galore \
            --paired \
            --gzip \
            --cores {threads} \
            --output_dir results/trimmed \
            {input.r1} \
            {input.r2}
        """


rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_R1_val_1.fq.gz",
        r2="results/trimmed/{sample}_R2_val_2.fq.gz"

    output:
        r1_html="results/qc/trimmed/fastqc/{sample}_R1_val_1_fastqc.html",
        r1_zip="results/qc/trimmed/fastqc/{sample}_R1_val_1_fastqc.zip",
        r2_html="results/qc/trimmed/fastqc/{sample}_R2_val_2_fastqc.html",
        r2_zip="results/qc/trimmed/fastqc/{sample}_R2_val_2_fastqc.zip"

    threads:
        config["threads"]["fastqc"]

    conda:
        "../envs/qc.yaml"

    shell:
        """
        mkdir -p results/qc/trimmed/fastqc

        fastqc \
            --threads {threads} \
            --outdir results/qc/trimmed/fastqc \
            {input.r1} \
            {input.r2}
        """


rule multiqc_trimmed:
    input:
        r1_html=expand(
            "results/qc/trimmed/fastqc/{sample}_R1_val_1_fastqc.html",
            sample=SAMPLES
        ),
        r1_zip=expand(
            "results/qc/trimmed/fastqc/{sample}_R1_val_1_fastqc.zip",
            sample=SAMPLES
        ),
        r2_html=expand(
            "results/qc/trimmed/fastqc/{sample}_R2_val_2_fastqc.html",
            sample=SAMPLES
        ),
        r2_zip=expand(
            "results/qc/trimmed/fastqc/{sample}_R2_val_2_fastqc.zip",
            sample=SAMPLES
        )

    output:
        report="results/qc/trimmed/multiqc_report.html"

    conda:
        "../envs/qc.yaml"

    shell:
        """
        mkdir -p results/qc/trimmed

        multiqc \
            results/qc/trimmed/fastqc \
            --outdir results/qc/trimmed \
            --force
        """