####### Configuration #######
configfile: "config/config.yaml"

OUTPUT_DIR = config["output_dir"]
SRA_ID = config["sra"]["sra_id"]

rule prefetch_data:
    output:
        sra_file = "results/sra/{sra_id}.sra"
    params:
        sra_id = "{sra_id}"
    log:
        OUTPUT_DIR + "/logs/{sra_id}_prefetch.log"
    shell:
        "prefetch {params.sra_id} --output-file {output.sra_file} > {log} 2>&1"

rule download_data:
    input:
        rules.prefetch_data.output.sra_file
    output: 
        r1 = temp(OUTPUT_DIR + "/{sra_id}_1.fastq"),
        r2 = temp(OUTPUT_DIR + "/{sra_id}_2.fastq") if lambda w: config["sra"].get("paired", False) else []
    params:
        sra_id = "{sra_id}",
        threads = config["sra"]["thread"]
    log:
        OUTPUT_DIR + "/logs/{sra_id}_download.log"
    shell:
        """
        fasterq-dump {input} \
            --outdir {OUTPUT_DIR} \
            --split-files \
            --threads {params.threads} > {log} 2>&1
        """

rule process_paired_data:
    input:
        r1 = rules.download_data.output.r1,
        r2 = rules.download_data.output.r2
    output: 
        filtered_r1 = OUTPUT_DIR + "/{sra_id}_filtered_1.fastq",
        filtered_r2 = OUTPUT_DIR + "/{sra_id}_filtered_2.fastq",
        report_json = OUTPUT_DIR + "/{sra_id}_fastp_report.json"
    params:
        threads = config["fastp"]["thread"],
        quality_threshold = config["fastp"]["qualified_quality_phred"],
        min_length = config["fastp"]["min_length"],
        detect_adapters = config["fastp"].get("detect_adapters", False)
    log:
        OUTPUT_DIR + "/logs/{sra_id}_fastp_paired.log"
    shell: 
        """
        fastp -i {input.r1} -I {input.r2} \
        -o {output.filtered_r1} -O {output.filtered_r2} \
        --thread {params.threads} \
        --detect_adapter_for_pe \
        -q {params.quality_threshold} \
        --length_required {params.min_length} \
        --json {output.report_json} > {log} 2>&1
        """
