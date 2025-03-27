####### Configuration #######

configfile: "config/config.yaml"

# Переменные из конфига можно загрузить глобально (опционально)
OUTPUT_DIR = config["output_dir"]
SRA_IDS = config["sra"]["sra_id"]


include: "workflow/rules/meta_genome.smk"
rule all:
    input:
        expand("{out_dir}/bracken_output.txt", out_dir=OUT_DIR),
        expand("{out_dir}/{sample}_report.tsv", out_dir=OUT_DIR, sample=SAMPLE)



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
        
