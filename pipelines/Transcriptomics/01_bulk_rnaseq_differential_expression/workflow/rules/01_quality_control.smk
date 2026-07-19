rule fastqc:

    input:
        "data/raw/{sample}.fastq.gz"

    output:
        html="results/qc/fastqc/{sample}_fastqc.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip"

    shell:
        """
        mkdir -p results/qc/fastqc

        fastqc \
        {input} \
        --outdir results/qc/fastqc
        """