####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/pre_prossecing.smk"
 
rule all:
    input:
        expand(
            os.path.join(OUTPUT_DIR, "{SRA_ID}_filtered_1.fastq"),
            sra=SRA_ID 
        ),
        expand(
            os.path.join(OUTPUT_DIR, "{SRA_ID}_filtered_2.fastq"),
            sra=SRA_ID
        )
        
