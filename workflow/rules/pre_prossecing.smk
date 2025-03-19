####### Configuration #######
configfile: "config/config.yaml"

# для выходных файлов
OUTPUT_DIR = config["output_dir"]
SRA_ID = config["sra"]["sra_id"]



# Создаем директорию для выходных файлов, если она не существует
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)

# предварительная загрузка 
rule prefetch_data:
    output:
        sra_file = "results/sra/{sra_id}.sra"
    params:
        sra_id = SRA_ID
    log:
        OUTPUT_DIR + "/logs/{sra_id}_prefetch.log"
    shell:
        """
        prefetch {params.sra_id} --output-file {output.sra_file} > {log} 2>&1
        """

ruleorder: process_paired_data > download_data

rule download_data:
    input:
        sra_file = "results/sra/{sra_id}.sra"
    output: 
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq",
        r2 = OUTPUT_DIR + "/{sra_id}_2.fastq" if config["sra"].get("paired", False) else temp(OUTPUT_DIR + "/{sra_id}_2.fastq")
    params:
        sra_id = SRA_ID,
        threads = config["sra"]["thread"],
        paired = config["sra"].get("paired", False),
        detect_adapters="--detect_adapter_for_pe"
    log:
        OUTPUT_DIR + "/logs/{sra_id}_download.log"
    shell: 
        """
        if {params.paired}; then
            fasterq-dump {params.sra_id} \
            --outdir {OUTPUT_DIR} \
            --split-files \
            --threads {params.threads} > {log} 2>&1
        else
            fasterq-dump {params.sra_id} \
            --outdir {OUTPUT_DIR} \
            --threads {params.threads} > {log} 2>&1
        fi
        """

rule process_paired_data:
    input:
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq",
        r2 = OUTPUT_DIR + "/{sra_id}_2.fastq"
    output: 
        filtered_r1 = OUTPUT_DIR + "/{sra_id}_filtered_1.fastq",
        filtered_r2 = OUTPUT_DIR + "/{sra_id}_filtered_2.fastq",
        report_json = OUTPUT_DIR + "/{sra_id}_fastp_report.json"
    params:
        threads = config["fastp"]["threads"],
        quality_threshold = config["fastp"]["qualified_quality_phred"],
        min_length = config["fastp"]["min_length"]
    log:
        OUTPUT_DIR + "/logs/{sra_id}_fastp_paired.log"
    shell: 
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.filtered_r1} -O {output.filtered_r2} \ 
        --threads {params.threads} \
        {params.detect_adapters} \
        -q {params.quality_threshold} \
        --length_required {params.min_length} \
        --json {output.report_json} > {log} 2>&1
        """

rule clean_temp_files:
    input: 
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq",
        r2 = OUTPUT_DIR + "/{sra_id}_2.fastq" if config["sra"].get("paired", False) else []
    params:
        sra_id = SRA_ID
    shell: 
        """

        rm -f {input.r1} {input.r2} 

        """

