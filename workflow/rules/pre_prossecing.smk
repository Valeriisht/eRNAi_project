###### Configuration #######
configfile: "config/config.yaml"
 
# for output files
OUTPUT_DIR = config["output_dir"]
SRA_ID = config["sra"]["sra_id"]

 
# create dir
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)
 
# pre-download
rule prefetch_data:
    output:
        sra_file = "metagenome/sra/{SRA_ID}.sra"
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
        sra_file = "metagenome/sra/{SRA_ID}.sra"
    output: 
        f1 = OUTPUT_DIR + "/raw_{SRA_ID}_1.fastq",
        r1 = OUTPUT_DIR + "/raw_{SRA_ID}_2.fastq" 
    params:
        sra_id = "{SRA_ID}",
        threads = config["sra"]["thread"],
        paired = config["sra"]["paired"]
    log:
        OUTPUT_DIR + "/logs/{SRA_ID}_download.log"
    shell: 
        """
        set -euo pipefail
        strace -e trace=file fasterq-dump {input.sra_file} \
        --outdir {OUTPUT_DIR} \
        --split-files  \
        --threads {params.threads} > {log} 2>&1

        # Переименование файлов
        mv {OUTPUT_DIR}/{params.sra_id}_1.fastq {output.f1}
        mv {OUTPUT_DIR}/{params.sra_id}_2.fastq {output.r1}
        """
 
rule process_paired_data:
    input:
        f1 = OUTPUT_DIR + "/raw_{SRA_ID}_1.fastq",
        r1 = OUTPUT_DIR + "/raw_{SRA_ID}_2.fastq" 
    output: 
        filtered_f1 = OUTPUT_DIR + "/{SRA_ID}_filtered_1.fastq",
        filtered_r2 = OUTPUT_DIR + "/{SRA_ID}_filtered_2.fastq", 
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
        fastp -i {input.f1} -I {input.r1} \
        -o {output.filtered_f1} -O {output.filtered_r2} \
        --thread {params.threads} \
        --detect_adapter_for_pe \
        -q {params.quality_threshold} \
        --length_required {params.min_length} \
        --json {output.report_json} > {log} 2>&1

        rm {input.f1} {input.r1}
        """