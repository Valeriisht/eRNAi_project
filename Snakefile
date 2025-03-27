####### Configuration #######

configfile: "config/config.yaml"

# Переменные из конфига можно загрузить глобально (опционально)
OUTPUT_DIR = config["output_dir"]
SRA_IDS = config["sra"]["sra_id"]


include: "workflow/rules/meta_genome.smk"

rule all:
    input:
        # Для Kraken2 + Bracken требуем финальный отчет и конвертированный файл
        expand(
            "{OUTPUT_DIR}/bracken_output.txt",
            out_dir=config["output_dir"]
        ),
        expand(
            "{OUTPUT_DIR}/{sample}_report.tsv",
            out_dir=config["output_dir"],
            sample=config["sample_name"]
        )


#rule all:
#    input:
#        GENOME_FILE

#rule all:
#    input:
#        expand(OUTPUT_DIR + "/{sra_id}_filtered_1.fastq", sra_id=SRA_IDS),
#        expand(OUTPUT_DIR + "/{sra_id}_filtered_2.fastq", sra_id=SRA_IDS)

#rule all:
#    input:
#        expand(OUTPUT_DIR + "/{taxid}/{sra_id}_quant_results", taxid=config["taxid"], sra_id=SRA_ID)
        
