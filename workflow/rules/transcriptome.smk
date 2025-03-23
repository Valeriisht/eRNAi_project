TRANSCRIPTOME_FASTA = config["output_dir"] + "/{taxid}.fna"  # Файл с транскриптомом
SRA_ID = config["sra"]["sra_id"]
INPUT_FASTQ_R1 = config["output_dir"] + "/{SRA_ID}_filtered_1.fastq"  # Файл с ридами (R1)
INPUT_FASTQ_R2 = config["output_dir"] + "/{SRA_ID}_filtered_2.fastq"  # Файл с ридами (R2) 
OUTPUT_DIR = config["output_dir"] + "/transcriptome_kallisto" # Выходная директория

# Правило для создания директорий
rule create_directories:
    output:
        directory(OUTPUT_DIR),
        directory("logs_kallisto")
    shell:
        """
        mkdir -p {output}
        """

# Индексирование транскриптома
rule kallisto_index:
    input:
        transcriptome = TRANSCRIPTOME_FASTA
    output:
        transcriptome_index = OUTPUT_DIR + "/{taxid}_transcriptome.idx"  # Используем {taxid}
    params:
        kmer_size = 31  # Опционально
    log:
        "logs_kallisto/{taxid}_kallisto_index.log"  # Лог-файл использует {taxid}
    shell:
        """
        kallisto index -i {output.transcriptome_index} -k {params.kmer_size} {input.transcriptome} > {log} 2>&1
        """

# Квази-выравнивание
rule kallisto_quant:
    input:
        index = OUTPUT_DIR + "/{taxid}_transcriptome.idx",  # Используем {taxid}
        r1 = INPUT_FASTQ_R1,
        r2 = INPUT_FASTQ_R2
    output:
        directory(OUTPUT_DIR + "/{taxid}/{SRA_ID}_quant_results")  # Используем {taxid} и {SRA_ID}
    params:
        bootstrap = 100
    log:
        "logs_kallisto/{taxid}/{SRA_ID}_kallisto_quant.log"  # Лог-файл использует {taxid} и {SRA_ID}
    run:
        if os.path.exists(input.r2):
            shell(
                """
                kallisto quant -i {input.index} -o {output} -b {params.bootstrap} {input.r1} {input.r2} > {log} 2>&1
                """
            )
        else:
            shell(
                """
                kallisto quant -i {input.index} -o {output} -b {params.bootstrap} {input.r1} > {log} 2>&1
                """
            )
