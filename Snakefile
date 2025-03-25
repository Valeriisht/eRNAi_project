####### Configuration #######
configfile: "config/config.yaml"


include: "workflow/rules/pre_prossecing.smk"

SRA_ID = config["sra"]["sra_id"]
rule all:
    input:
        # Генерируем ВСЕ выходные файлы для всех SRA_ID
        expand(
            OUTPUT_DIR + "/{sra_id}_filtered_1.fastq", 
            sra_id=SRA_IDS
        ),
        expand(
            OUTPUT_DIR + "/{sra_id}_filtered_2.fastq", 
            sra_id=SRA_IDS
        )



#rule all:
#    input:
#        expand(OUTPUT_DIR + "/{taxid}/{sra_id}_quant_results", taxid=config["taxid"], sra_id=SRA_ID)
        
