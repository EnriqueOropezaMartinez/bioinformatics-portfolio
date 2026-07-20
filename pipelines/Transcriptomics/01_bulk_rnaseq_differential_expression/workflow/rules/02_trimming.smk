rule trim_galore:

    input:
        fastq="data/raw/{sample}.fastq.gz"

    output:
        trimmed="results/trimmed/{sample}_trimmed.fq.gz"

    threads: 4

    shell:
        """
        mkdir -p results/trimmed

        trim_galore \
        --gzip \
        --fastqc \
        --cores {threads} \
        --output_dir results/trimmed \
        {input.fastq}
        """