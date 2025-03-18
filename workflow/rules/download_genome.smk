# Snakefile
configfile: "config/config.yaml"

TAXID = config["taxid"]
OUTPUT_DIR = config["output_dir"]
GENOME_FILE = f"{OUTPUT_DIR}/{TAXID}.fna"

#Загрузка архива
rule download_genome:
    output:
        zip = f"{OUTPUT_DIR}/{TAXID}_genome.zip"
    params:
        taxid = TAXID
    log:
        f"{OUTPUT_DIR}/logs/download.log"
    shell:
        "datasets download genome taxon {params.taxid} --reference --filename {output.zip} > {log} 2>&1"
#Распаковка
rule extract_genome:
    input:
        download_genome.output.zip
    output:
        directory(f"{OUTPUT_DIR}/{TAXID}_genome")
    log:
        f"{OUTPUT_DIR}/logs/extract.log"
    shell:
        "unzip {input} -d {output} > {log} 2>&1"
#Поиск и переименование 
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
        # Поиск файла .fna в директории
        GENOME_FILE=$(find {input} -name '*.fna' -print -quit)
        if [[ -z "$GENOME_FILE" ]]; then
            echo "Геном для taxid {params.taxid} не найден" >> {log}
            exit 1
        fi
        # Переименовываем файл
        mv "$GENOME_FILE" {output} >> {log} 2>&1
        """

rule clean_temp:
    input:
        GENOME_FILE
    output:
        touch(f"{OUTPUT_DIR}/clean.done")
    shell:
        "rm -rf {OUTPUT_DIR}/{TAXID}_genome {OUTPUT_DIR}/{TAXID}_genome.zip"
