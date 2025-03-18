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
        extract_genome.output,  # Используем output из правила extract_genome
        f"{OUTPUT_DIR}/logs"
    output:
        GENOME_FILE
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/rename.log"
    run:
        import subprocess
        import os
        genome = subprocess.getoutput(f"find {input[0]} -name '*.fna' -print -quit")
        if not genome:
            raise ValueError(f"Геном для taxid {params.taxid} не найден")
        subprocess.run(f"mv '{genome}' '{output}' > {log} 2>&1", shell=True)

rule clean_temp:
    input:
        GENOME_FILE
    output:
        touch(f"{OUTPUT_DIR}/clean.done")
    shell:
        "rm -rf {OUTPUT_DIR}/{TAXID}_genome {OUTPUT_DIR}/{TAXID}_genome.zip"
