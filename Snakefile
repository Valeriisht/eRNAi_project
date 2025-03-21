####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/pre_prossecing.smk"
 
rule all:
    input:
        expand(
            os.path.join(OUTPUT_DIR, "{sra_id}_filtered_1.fastq"),
            sra_id = "SRR8265535" 
        ),
        expand(
            os.path.join(OUTPUT_DIR, "{sra_id}_filtered_2.fastq"),
            sra_id = "SRR8265535"
        )
        
