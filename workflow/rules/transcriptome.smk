TRANSCRIPTOME_FASTA = config["output_dir"] + "/genome_ref/{taxid}.fna"  # Файл с транскриптомом
INPUT_FASTQ_R1 = config["output_dir"] + "/transcriptome/{sra_id}_filtered_1.fastq"  # Файл с ридами (R1)
INPUT_FASTQ_R2 = config["output_dir"] + "/transcriptome/{sra_id}_filtered_2.fastq"  # Файл с ридами (R2) 
OUTPUT_DIR = "transcriptome_kallisto" # Выходная директория
SRA_ID = config["sra"]["sra_id"]

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
        index = "transcriptome_kallisto/9031_transcriptome.idx",  # Используем {taxid}
        r1 = INPUT_FASTQ_R1,
        r2 = INPUT_FASTQ_R2
    output:
        directory(OUTPUT_DIR + "/{taxid}/{sra_id}_quant_results")  # Используем {taxid} и {sra_id}
    params:
        bootstrap = 100
    log:
        "logs_kallisto/{taxid}/{sra_id}_kallisto_quant.log"  # Лог-файл использует {taxid} и {sra_id}
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
