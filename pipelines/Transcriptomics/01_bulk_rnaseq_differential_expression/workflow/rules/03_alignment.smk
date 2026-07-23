rule star_index:
    input:
        genome=config["reference"]["genome"],
        annotation=config["reference"]["annotation"]

    output:
        directory(config["reference"]["star_index"])

    threads:
        config["threads"]["star_index"]

    conda:
        "../envs/alignment.yaml"

    params:
        sjdb_overhang=config["reference"].get("sjdb_overhang", 149)

    shell:
        """
        mkdir -p {output}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.annotation} \
            --sjdbOverhang {params.sjdb_overhang}
        """


rule star_alignment:
    input:
        index=config["reference"]["star_index"],
        r1="results/trimmed/{sample}_R1_val_1.fq.gz",
        r2="results/trimmed/{sample}_R2_val_2.fq.gz"

    output:
        bam="results/alignment/{sample}/{sample}.sorted.bam",
        log_final="results/alignment/{sample}/{sample}.Log.final.out",
        log_out="results/alignment/{sample}/{sample}.Log.out",
        log_progress="results/alignment/{sample}/{sample}.Log.progress.out",
        splice_junctions="results/alignment/{sample}/{sample}.SJ.out.tab"

    threads:
        config["threads"]["star_alignment"]

    conda:
        "../envs/alignment.yaml"

    params:
        output_prefix="results/alignment/{sample}/{sample}."

    shell:
        """
        mkdir -p results/alignment/{wildcards.sample}

        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM SortedByCoordinate

        mv \
            {params.output_prefix}Aligned.sortedByCoord.out.bam \
            {output.bam}
        """


rule samtools_index:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam"

    output:
        bai="results/alignment/{sample}/{sample}.sorted.bam.bai"

    threads:
        config["threads"]["samtools"]

    conda:
        "../envs/alignment.yaml"

    shell:
        """
        samtools index \
            -@ {threads} \
            {input.bam} \
            {output.bai}
        """


rule samtools_flagstat:
    input:
        bam="results/alignment/{sample}/{sample}.sorted.bam"

    output:
        report="results/alignment/{sample}/{sample}.flagstat.txt"

    threads:
        config["threads"]["samtools"]

    conda:
        "../envs/alignment.yaml"

    shell:
        """
        samtools flagstat \
            -@ {threads} \
            {input.bam} \
            > {output.report}
        """