# Snakefile
configfile: "config/config.yaml"

TAXID = config["taxid"]
OUTPUT_DIR = config["output_dir"]
GENOME_FILE = f"{OUTPUT_DIR}/{TAXID}.fna"

# download archive
rule download_genome:
    output:
        zip = f"{OUTPUT_DIR}/{TAXID}_genome.zip"
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/download.log"
    shell:
        "datasets download genome taxon {params.taxid} --reference --filename {output.zip} > {log} 2>&1"

# unpacking
rule extract_genome:
    input:
        rules.download_genome.output.zip
    output:
        directory(f"{OUTPUT_DIR}/{TAXID}_genome")
    log:
        f"{OUTPUT_DIR}/logs/extract.log"
    shell:
        "unzip {input} -d {output} > {log} 2>&1"

# rename and find
rule find_rename_genome:
    input:
        rules.extract_genome.output
    output:
        GENOME_FILE
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/rename.log"
    shell:
        """
        # Searching for .fna file in the directory
        echo "Searching for .fna file in directory {input}" > {log}
        GENOME_FILE=$(find "{input}" -name '*.fna' -print -quit)
        if [[ -z "$GENOME_FILE" ]]; then
            echo "Genome for taxid {params.taxid} not found" >> {log}
            exit 1
        fi
        echo "Found file: $GENOME_FILE" >> {log}
        mv "$GENOME_FILE" {output} >> {log} 2>&1
        """

rule clean_temp:
    input:
        GENOME_FILE
    output:
        touch(f"{OUTPUT_DIR}/clean.done")
    shell:
        "rm -rf {OUTPUT_DIR}/{TAXID}_genome {OUTPUT_DIR}/{TAXID}_genome.zip"
