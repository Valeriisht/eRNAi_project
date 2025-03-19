####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/meta_genome.smk"


rule all:
     input:
         f"{config["output_dir"]}/{config["sample_name"]}_report.tsv" 
