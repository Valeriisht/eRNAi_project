####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/pre_prossecing.smk"

rule all:
     input:
        filtered_r1 = OUTPUT_DIR + "/{SRA_ID}_filtered_1.fastq",
        filtered_r2 = OUTPUT_DIR + "/{SRA_ID}_filtered_2.fastq",
        report_json = OUTPUT_DIR + "/{SRA_ID}_fastp_report.json"
