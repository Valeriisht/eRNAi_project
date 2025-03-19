####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/meta_genome.smk"


rule all:
     input:
         f"{OUT_DIR}/{SAMPLE}_report.tsv" 
