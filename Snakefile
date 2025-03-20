####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/download_genome.smk"
 
rule all:
    input:
        GENOME_FILE
