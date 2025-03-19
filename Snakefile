####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/pre_prossecing.smk"


rule all:
     input:
         f"{config['output_dir']}/{config['sample_name']}_report.tsv" 
