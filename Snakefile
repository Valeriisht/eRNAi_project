configfile: "config/config.yaml"

# Global params
INP_DIR = config["input_dir"]
OUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"]
SRA_IDS = config["sra"]["sra_id"]
DB = config["database"]
ALGO = config["algorithm"]
READ_LEN = config["read_length"]
LEVELS = config["taxonomic_level"]


# =============================================  
rule for meta_genome
include: "workflow/rules/meta_genome.smk"
rule all:
    input:
        expand(f"{OUT_DIR}/kraken2_{{sra_id}}.report", sra_id=SRA_IDS),
        expand(f"{OUT_DIR}/bracken_{{sra_id}}_output_{{level}}.report", 
               sra_id=SRA_IDS, level=LEVELS),
        f"{OUT_DIR}/{SAMPLE}_report.tsv"
# =============================================        
rule for host_community 
include: "workflow/rules/host_community.smk"
rule all:
        input:
                expand(
                        "{out_dir}/{sra_id}_{type}_reads.fastq.gz",
                        out_dir=OUT_DIR,
                        sra_id=SRA_IDS,
                        type=["host", "metagenome"])
# =============================================  
rule for pre_prossecing
include: "workflow/rules/pre_prossecing.smk"
rule all:
    input:
        expand(OUTPUT_DIR + "/{sra_id}_filtered_1.fastq", sra_id=SRA_IDS),
        expand(OUTPUT_DIR + "/{sra_id}_filtered_2.fastq", sra_id=SRA_IDS)
# =============================================
 rule for transcriptome
include: "workflow/rules/transcriptome.smk"
rule all:
        input:
                expand(OUTPUT_DIR + "/{taxid}/{sra_id}_quant_results", taxid=config["taxid"], sra_id=SRA_ID)
# =============================================
rule for meta_genome
wildcard_constraints:
  level="|".join(LEVELS)  

include: "workflow/rules/meta_genome.smk"

rule all:
    input:
        expand(f"{OUT_DIR}/kraken2_{sra_id}.report", sra_id=SRA_IDS),
        expand(f"{OUT_DIR}/bracken_{sra_id}_output_S.report", sra_id=SRA_IDS),
        f"{OUT_DIR}/{SAMPLE}_report.tsv"


        expand("{out_dir}/{sample}_report.tsv", 
              out_dir=OUT_DIR,
              sample=SAMPLE)
        
# =============================================
rule for download_genome
include: "workflow/rules/download_genome.smk"
rule all:
        input:
                GENOME_FILE
# =============================================



        
