# Snakefile
configfile: "config/config.yaml"

TAXID = config["taxid"]
OUTPUT_DIR = config["output_dir"]
GENOME_FILE = f"{OUTPUT_DIR}/{TAXID}.fna"

rule create_dirs:
    output:
        directory(f"{OUTPUT_DIR}/logs")
    shell:
        "mkdir -p {output}"

#Загрузка архива
rule download_genome:
    output:
        zip = temp(f"{OUTPUT_DIR}/{TAXID}_genome.zip")
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/download.log"
    shell:
        "datasets download genome taxon {params.taxid} --reference --filename {output.zip} > {log} 2>&1"
#Распаковка
rule extract_genome:
    input:
        rules.download_genome.output.zip
    output:
        directory(f"{OUTPUT_DIR}/{TAXID}_genome")
    log:
        f"{OUTPUT_DIR}/logs/extract.log"
    shell:
        "unzip {input} -d {output} > {log} 2>&1"

# Поиск и переименование 
# Snakefile
configfile: "config/config.yaml"

TAXID = config["taxid"]
OUTPUT_DIR = config["output_dir"]
GENOME_FILE = f"{OUTPUT_DIR}/{TAXID}.fna"

rule create_dirs:
    output:
        directory(f"{OUTPUT_DIR}/logs")
    shell:
        "mkdir -p {output}"

#Загрузка архива
rule download_genome:
    output:
        zip = temp(f"{OUTPUT_DIR}/{TAXID}_genome.zip")
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/download.log"
    shell:
        "datasets download genome taxon {params.taxid} --reference --filename {output.zip} > {log} 2>&1"
#Распаковка
rule extract_genome:
    input:
        rules.download_genome.output.zip
    output:
        directory(f"{OUTPUT_DIR}/{TAXID}_genome")
    log:
        f"{OUTPUT_DIR}/logs/extract.log"
    shell:
        "unzip {input} -d {output} > {log} 2>&1"

# Поиск и переименование 
rule find_rename_genome:
    input:
        rules.extract_genome.output,
        f"{OUTPUT_DIR}/logs"  
    output:
        GENOME_FILE
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/rename.log"
    run:
        import subprocess
        genome = subprocess.getoutput(f"find {input} -name '*.fna' -print -quit")
        if not genome:
            raise ValueError(f"Геном для taxid {params.taxid} не найден")
        shell(f"mv '{genome}' '{output}' > {log} 2>&1")

rule clean_temp:
    input:
        GENOME_FILE
    output:
        touch(f"{OUTPUT_DIR}/clean.done")
    shell:
        "rm -rf {OUTPUT_DIR}/{TAXID}_genome {OUTPUT_DIR}/{TAXID}_genome.zip"
rule clean_temp:
    input:
        GENOME_FILE
    output:
        touch(f"{OUTPUT_DIR}/clean.done")
    shell:
        "rm -rf {OUTPUT_DIR}/{TAXID}_genome {OUTPUT_DIR}/{TAXID}_genome.zip"
