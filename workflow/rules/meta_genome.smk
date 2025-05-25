# Snakefile
configfile: "config/config.yaml"

# config 
SRA_IDS = config["sra"]["sra_id"]
ALGO = config["algorithm"]       # "kraken2" or "metaphlan"
DB = config.get("database", "")  # path to Kraken2 


INP_DIR = config["input_dir"]
OUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"]   # for output file

wildcard_constraints:
    level="S"  # levels: S, G, P

def INPUT_R1(wildcards):
    return f"{INP_DIR}/{wildcards.sra_id}_metagenome_reads.fastq"


### 1: Kraken2 + Bracken ###
if ALGO == "kraken2":
    if not DB:
        raise ValueError("Kraken2 database required")

    rule kraken2_classify:
        input:
            r1 = INPUT_R1,
        output:
            report = temp(f"{OUT_DIR}/kraken2_{{sra_id}}.report"),
            raw = temp(f"{OUT_DIR}/kraken2_output_{{sra_id}}.report")
        threads: 8
        log:
            f"{OUT_DIR}/logs/kraken2_{{sra_id}}.log"
        shell:
            """
            kraken2 --db {DB} --threads {threads} {input.r1} \
                   --output {output.raw} --report {output.report} > {log} 2>&1
            """

    rule bracken_abundance:
        input:
            rules.kraken2_classify.output.report
        output:
            f"{OUT_DIR}/bracken_{{sra_id}}_output_{{level}}.report"
        params:
            db = DB,
            readlen = READ_LEN,
            threshold = 10
        log:
            f"{OUT_DIR}/logs/bracken_{{sra_id}}_{{level}}.log"
        shell:
            """
            bracken -d {params.db} -i {input} -o {output} \
                    -r {params.readlen} -l {wildcards.level} \
                    -t {params.threshold} > {log} 2>&1
            """

    rule convert_to_mpa:
        input:
            expand(f"{OUT_DIR}/bracken_{{sra_id}}_output_S.report", sra_id=SRA_IDS)
        output:
            f"{OUT_DIR}/{SAMPLE}_report.tsv"
        shell:
            "kreport2mpa.py -r {input} -o {output} --display-header"

### 2: Metaphlan ###
elif ALGO == "metaphlan":
    rule metaphlan_profile:
        input:
            r1 = INPUT_R1,
            r2 = INPUT_R2
        output:
            f"{OUT_DIR}/{SAMPLE}_report.tsv"
        log:
            f"{OUT_DIR}/logs/metaphlan.log"
        shell:
            """
            # taxonomic profile
            metaphlan {input.r1},{input.r2} --input_type fastq \
                      --nproc 8 -o {output} > {log} 2>&1
            """

else:
    raise ValueError("Unsupported algorithm. Valid values: kraken2, metaphlan")
