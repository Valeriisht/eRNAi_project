####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/download_genome.smk"
 
rule all:
    input:
        # Используем expand, если GENOME_FILE содержит wildcards (например, {sample})
        expand(
            config["GENOME_FILE"],  # Путь к геному из конфига
            # Дополнительные параметры, если требуется (например, sample=config["samples"])
        )
