####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/download_genome.smk"

rule all:
    input:
        directory(OUTPUT_DIR),  # Директория transcriptome_kallisto
        expand(OUTPUT_DIR + "/{taxid}_transcriptome.idx", taxid=config["taxid"]),  # Индекс транскриптома
        expand(OUTPUT_DIR + "/{taxid}/{sra_id}_quant_results", taxid=config["taxid"], sra_id=SRA_ID)  # Результаты квази-выравнивания
        
