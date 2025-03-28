####### Configuration #######

configfile: "config/config.yaml"

# Переменные из конфига можно загрузить глобально (опционально)
OUT_DIR = config["output_dir"]
# OUTPUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"]
SRA_ID = config["sra"]["sra_id"]
LEVEL = config['taxonomic_level']

INPUT_R1 = expand(
    "{out_dir}/{sra_id}_filtered_1.fastq",
    out_dir=OUT_DIR,
    sra_id=SRA_ID
)

INPUT_R2 = expand(
    "{out_dir}/{sra_id}_filtered_2.fastq",
    out_dir=OUT_DIR,
    sra_id=SRA_ID
)
include: "workflow/rules/meta_genome.smk"

rule all:
    input:
        expand({out_dir}/kraken2_{sra_id}.report"),sra_id=SRA_ID,out_dir=OUT_DIR),
        expand({out_dir}/kraken2_output_{sra_id}.report"),sra_id=SRA_ID,out_dir=OUT_DIR)
        #expand("{out_dir}/bracken_{sra_id}_output_{level}.report", out_dir=OUT_DIR,level = LEVEL,sra_id=SRA_ID)


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
        
