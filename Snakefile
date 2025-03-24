####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/transcriptome.smk"
 
rule all:
    input:
        expand(OUTPUT_DIR + "/{taxid}/{sra_id}_quant_results", taxid=config["taxid"], sra_id=SRA_ID)
        
