###### Configuration #######
configfile: "config/config.yaml"
 
# для выходных файлов
OUTPUT_DIR = config["output_dir"]
SRA_ID = "SRR8265535"
 
# Создаем директорию для выходных файлов, если она не существует
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)
 
# предварительная загрузка 
rule prefetch_data:
    output:
        sra_file = "results/sra/{SRA_ID}.sra"
    params:
        sra_id = "{SRA_ID}"
    log:
        OUTPUT_DIR + "/logs/{SRA_ID}_prefetch.log"
    shell:
        """
        prefetch {params.sra_id} --output-file {output.sra_file} > {log} 2>&1
        """

rule download_data:
    input:
        sra_file = "results/sra/{SRA_ID}.sra"
    output: 
        f1 = OUTPUT_DIR + "/{SRA_ID}_1.fastq",
        r1 = OUTPUT_DIR + "/{SRA_ID}_2.fastq" 
    params:
        sra_id = "{SRA_ID}",
        threads = config["sra"]["thread"],
        paired = config["sra"]["paired"]
    log:
        OUTPUT_DIR + "/logs/{SRA_ID}_download.log"
    shell: 
        """
        set -euo pipefail
        fasterq-dump {input.sra_file} \
        --outdir {OUTPUT_DIR} \
        --split-files  \
        --threads {params.threads} > {log} 2>&1
        """
 
rule process_paired_data:
    input:
        f1 = OUTPUT_DIR + "/{SRA_ID}_1.fastq",
        r1 = OUTPUT_DIR + "/{SRA_ID}_2.fastq" 
    output: 
        filtered_r1 = OUTPUT_DIR + "/{SRA_ID}_filtered_1.fastq",
        filtered_r2 = OUTPUT_DIR + "/{SRA_ID}_filtered_2.fastq" if config["sra"].get("paired", False) else [], 
        report_json = OUTPUT_DIR + "/{SRA_ID}_fastp_report.json"
    params:
        threads = config["fastp"]["thread"],
        quality_threshold = config["fastp"]["qualified_quality_phred"],
        min_length = config["fastp"]["min_length"],
        detect_adapters = config["fastp"]["detect_adapters"]
    log:
        OUTPUT_DIR + "/logs/{SRA_ID}_fastp_paired.log"
    shell: 
        """
        fastp -i {input.r1_f} -I {input.r2_r} \
        -o {output.filtered_r1} -O {output.filtered_r2} \
        --thread {params.threads} \
        --detect_adapter_for_pe \
        -q {params.quality_threshold} \
        --length_required {params.min_length} \
        --json {output.report_json} > {log} 2>&1
        """
