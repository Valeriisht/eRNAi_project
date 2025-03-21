####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/download_genome.smk"
 
rule all:
    input:
        expand(
            config["GENOME_FILE"],  # Путь к геному из конфига
            sample=config["sra"]["sra_id"])
        )
