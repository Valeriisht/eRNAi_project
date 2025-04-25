####### Configuration #######
configfile: "config/config.yaml"

include: "workflow/rules/transcriptome.smk"


OUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"]
SRA_ID = config["sra"]["sra_id"] 


# include: "workflow/rules/transcriptome.smk"
include: "workflow/rules/host_community.smk"

rule all:
    input:
        expand(
            "{out_dir}/{sra_id}_{type}_reads.fastq.gz",
            out_dir=OUT_DIR,
            sra_id=SRA_ID,
            type=["host", "metagenome"]
        )
# rule all:
#     input:
#         expand(
#             os.path.join(OUTPUT_DIR, "{SRA_ID}_filtered_1.fastq"),
#             SRA_ID=config["sra"]["sra_id"]  # Список SRA_ID из config.yaml
#         ),
#         expand(
#             os.path.join(OUTPUT_DIR, "{SRA_ID}_filtered_2.fastq"),
#             SRA_ID=config["sra"]["sra_id"]
#         )


#rule all:
#    input:
#        expand(OUT_DIR + "/kraken2_{sra_id}.report",sra_id=SRA_ID),
#        expand(OUT_DIR + "/kraken2_output_{sra_id}.report",sra_id=SRA_ID),
#        expand(OUT_DIR + "/bracken_{sra_id}_output_{level}.report",level = LEVEL,sra_id=SRA_ID)
        



        #expand("{out_dir}/{sample}_report.tsv", 
        #      out_dir=OUT_DIR,
        #      sample=SAMPLE)
        



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

# rule all:
#    input:
#        host_reads = expand("{output_dir}/host_reads.fastq.gz", output_dir=config["output_dir"]),
#        metagenome_reads = expand("{output_dir}/metagenome_reads.fastq.gz", output_dir=config["output_dir"])

# # INPUT_R1 = expand(
#      "{out_dir}/{sra_id}_filtered_1.fastq",
#      out_dir=OUT_DIR,
#      sra_id=SRA_ID
#      )

# INPUT_R2 = expand(
#      "{out_dir}/{sra_id}_filtered_2.fastq",
#      out_dir=OUT_DIR,
#      sra_id=SRA_ID
#      )
        
