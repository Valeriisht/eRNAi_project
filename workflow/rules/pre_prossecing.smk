
####### Configuration #######
configfile: "enviromental.yaml"

SRA_ID = $1 
OUTPUT_DIR = "processed_data" 

# Создаем директорию для выходных файлов, если она не существует
shell:
    (
    """

    mkdir -p {OUTPUT_DIR}

    """
    )

rule download_data:
    output: 
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq"
        r2 = OUTPUT_DIR + "/{sra_id}_2.fastq"
    params:
        sra_id = SRA_ID
    log:
        OUTPUT_DIR + "/logs/download.log"
    run: 
        fasterq-dump {params.sra_id} --outdir {OUTPUT_DIR}  --split-files > {log} 2>&1 # --threads 8


rule process_paired_data:
    input:
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq"
        r2 = OUTPUT_DIR + "/{sra_id}_2.fastq"
    output: 
        filtered_r1 = OUTPUT_DIR + "/{sra_id}_filtered_1.fastq"
        filtered_r2 = OUTPUT_DIR + "/{sra_id}_filtered_2.fastq"
        report_json = OUTPUT_DIR + "/{sra_id}_fastp_report.json"
    params:
        sra_id = SRA_ID
    log:
        OUTPUT_DIR + "/logs/fastp_paired.log"
    shell: 
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.filtered_r1} -O {output.filtered_r2} \ 
              --threads 8
              --detect_adapter_for_pe
              -q 15
              --length_required 50
              --correction
              --json {output.report_json} > {log} 2>&1

        """

rule process_single_data:
    input:
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq"
    output: 
        filtered_r1 = OUTPUT_DIR + "/{sra_id}_filtered_1.fastq"
        report_json = OUTPUT_DIR + "/{sra_id}_fastp_report.json"
    params:
        sra_id = SRA_ID
    log:
        OUTPUT_DIR + "/logs/fastp_single.log"
    shell: 
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.filtered_r1} -O {output.filtered_r2} \ 
              --threads 8
              --detect_adapter_for_pe
              -q 15
              --length_required 50
              --correction
              --json {output.report_json} > {log} 2>&1
        """

rule clean_temp_files:
    input: 
        r1 = OUTPUT_DIR + "/{sra_id}_1.fastq"
        r2 = OUTPUT_DIR + "/{sra_id}_2.fastq"
    params:
        sra_id = SRA_ID
    shell: 
        """

        rm -f {input.r1} {input.r2} 

        """

