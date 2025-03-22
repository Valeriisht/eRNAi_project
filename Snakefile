####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/pre_prossecing.smk"
 
rule all:
    input:
        expand(
            os.path.join(OUTPUT_DIR, "{SRA_ID}_filtered_1.fastq"),
            SRA_ID = "SRR8265535" 
        ),
        expand(
            os.path.join(OUTPUT_DIR, "{SRA_ID}_filtered_2.fastq"),
            SRA_ID = "SRR8265535"
        )
        
