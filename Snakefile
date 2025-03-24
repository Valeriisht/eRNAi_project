####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/pre_prossecing.smk"
 
 rule all:
     input:
         expand(
             os.path.join(OUTPUT_DIR, "{sra}_filtered_1.fastq"),
             sra=config["sra"]["sra_id"] 
         ),
         expand(
             os.path.join(OUTPUT_DIR, "{sra}_filtered_2.fastq"),
             sra=config["sra"]["sra_id"] 
         )
        
