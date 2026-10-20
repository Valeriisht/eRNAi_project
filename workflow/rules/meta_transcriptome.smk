# Snakefile
configfile: "config/config.yaml"

# config 
SRA_IDS = config["metatranscriptome"]["sra"]["sra_id"]
INP_DIR = config["metatranscriptome"]["input_dir"]
OUT_DIR = config["metatranscriptome"]["output_dir"]

SRA_DIR = "metatranscriptome/sra"

rule all:
    input:
        expand(
            OUT_DIR + "/fastp/{SRA_ID}_1.fastq.gz",
            SRA_ID=SRA_IDS
        ),
        expand(
            OUT_DIR + "/fastp/{SRA_ID}_2.fastq.gz",
            SRA_ID=SRA_IDS
        ),
        expand(
            OUT_DIR + "/fastp/{SRA_ID}_fastp.html",
            SRA_ID=SRA_IDS
        )

        

# pre-download
rule prefetch_metatranscriptome:
    output:
        sra_file = SRA_DIR + "/{SRA_ID}.sra"  # проверить 
    params:
        sra_id = "{SRA_ID}"
    log:
        OUT_DIR + "/logs/{SRA_ID}_prefetch.log"
    shell:
        """
        set -euo pipefail

        mkdir -p {SRA_DIR}
        mkdir -p {OUT_DIR}/logs

        prefetch {params.sra_id} \
            --output-file {output.sra_file} \
            > {log} 2>&1
        """

rule download_metatranscriptome:
    input:
        sra_file = SRA_DIR + "/{SRA_ID}.sra"
    output: 
        r1 = INP_DIR + "/{SRA_ID}_1.fastq",
        r2 = INP_DIR + "/{SRA_ID}_2.fastq"
    # выдавало ошибку надо исправить !!!! пока захэшила чтоб не мешало
    threads:
        4
    params:
        outdir = INP_DIR
    
    log:
        OUT_DIR + "/logs/{SRA_ID}_download.log"
    shell: 
        """
        set -euo pipefail

        mkdir -p {params.outdir}
        mkdir -p {OUT_DIR}/logs

        fasterq-dump {input.sra_file} \
            --outdir {params.outdir} \
            --split-files \
            --threads {threads} \
            > {log} 2>&1

        test -f {output.r1}
        test -f {output.r2}
        """

rule fastp_metatranscriptome:
    input:
        r1 = INP_DIR + "/{SRA_ID}_1.fastq",
        r2 = INP_DIR + "/{SRA_ID}_2.fastq"
    output:
        r1 = OUT_DIR + "/fastp/{SRA_ID}_1.fastq.gz",
        r2 = OUT_DIR + "/fastp/{SRA_ID}_2.fastq.gz",
        html = OUT_DIR + "/fastp/{SRA_ID}_fastp.html",
        json = OUT_DIR + "/fastp/{SRA_ID}_fastp.json"
    threads:
        config["fastp"]["thread"]
    params:
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"],
        min_length = config["fastp"]["min_length"]
    log:
        OUT_DIR + "/logs/{SRA_ID}_fastp.log"
    shell:
        """
        set -euo pipefail

        mkdir -p {OUT_DIR}/fastp
        mkdir -p {OUT_DIR}/logs

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --thread {threads} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --length_required {params.min_length} \
            --html {output.html} \
            --json {output.json} \
            > {log} 2>&1
        """