configfile: "config/config.yaml"

TAXID = config["taxid"]
OUTPUT_DIR = config["output_dir"]
GENOME_FILE = f"{OUTPUT_DIR}/{TAXID}.fna"

rule create_dirs:
    output:
        directory(f"{OUTPUT_DIR}/logs")
    shell:
        "mkdir -p {output}"

rule download_genome:
    output:
        zip = temp(f"{OUTPUT_DIR}/{TAXID}_genome.zip")
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/download.log"
    shell:
        "datasets download genome taxon {params.taxid} --reference --filename {output.zip} > {log} 2>&1"

rule extract_genome:
    input:
        download_genome.output.zip  # Используем output из правила download_genome
    output:
        directory(f"{OUTPUT_DIR}/{TAXID}_genome")
    log:
        f"{OUTPUT_DIR}/logs/extract.log"
    shell:
        "unzip {input} -d {output} > {log} 2>&1"

rule find_rename_genome:
    input:
        genome_dir = extract_genome.output,
        logs = f"{OUTPUT_DIR}/logs"
    output:
        GENOME_FILE
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/rename.log"
    shell:
        """
        GENOME_FILE=$(find {input.genome_dir} -name '*.fna' -print -quit)
        if [[ -z "$GENOME_FILE" ]]; then
            echo "Геном для taxid {params.taxid} не найден" > {log}
            exit 1
        fi
        mv "$GENOME_FILE" {output} >> {log} 2>&1
        """

rule clean_temp:
    input:
        GENOME_FILE
    output:
        touch(f"{OUTPUT_DIR}/clean.done")
    shell:
        "rm -rf {OUTPUT_DIR}/{TAXID}_genome {OUTPUT_DIR}/{TAXID}_genome.zip"
