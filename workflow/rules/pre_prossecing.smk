###### Configuration #######
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
        sra_id = "{sra_id}"
    log:
        OUTPUT_DIR + "/logs/{sra_id}_prefetch.log"
    shell:
        """
        prefetch {params.sra_id} --output-file {output.sra_file} > {log} 2>&1
        """
 
rule download_data:
    input:
        sra_file = "results/sra/{sra_id}.sra"
    output: 
        r1 = temp(OUTPUT_DIR + "/down_{sra_id}_1.fastq"),
        r2 = temp(OUTPUT_DIR + "/down_{sra_id}_2.fastq") if config["sra"].get("paired", False) else []
    params:
        sra_id = "{sra_id}",
        threads = config["sra"]["thread"],
        paired = config["sra"].get("paired", False),
    log:
        OUTPUT_DIR + "/logs/{sra_id}_download.log"
    shell:
        """
        set -euo pipefail
        fasterq-dump {params.sra_id} \
            --outdir {OUTPUT_DIR} \
            {{ "--split-files" if params.paired else "" }} \
            --threads {params.threads} > {log} 2>&1
        """
 
rule process_paired_data:
    input:
        r1_f = OUTPUT_DIR + "/down_{sra_id}_1.fastq",
        r2_r = OUTPUT_DIR + "/down_{sra_id}_2.fastq" if config["sra"].get("paired", False) else []
    output: 
        filtered_r1 = OUTPUT_DIR + "/{sra_id}_filtered_1.fastq",
        filtered_r2 = OUTPUT_DIR + "/{sra_id}_filtered_2.fastq" if config["sra"].get("paired", False) else [], 
        report_json = OUTPUT_DIR + "/{sra_id}_fastp_report.json"
    params:
        threads = config["fastp"]["thread"],
        quality_threshold = config["fastp"]["qualified_quality_phred"],
        min_length = config["fastp"]["min_length"],
        detect_adapters = config["fastp"]["detect_adapters"]
    log:
        OUTPUT_DIR + "/logs/{sra_id}_fastp_paired.log"
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
