# Snakefile
configfile: "config/config.yaml"

# Конфигурация
SRA_ID = config["sra"]["sra_id"]
ALGO = config["algorithm"]       # "kraken2" или "metaphlan"
DB = config.get("database", "")  # Путь к базе Kraken2 (обязательно для алгоритма kraken2)
 # Парные fastq.gz файлы
OUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"]   # Имя образца (для выходных файлов)

input = expand(
    "{out_dir}/{sra_id}_filtered_1.fastq",
    "{out_dir}/{sra_id}_filtered_2.fastq",
    out_dir=OUT_DIR, 
    sra_id=SRA_IDS
)

INPUT_R1, INPUT_R2 = config["input"] 

### Вариант 1: Kraken2 + Bracken ###
if ALGO == "kraken2":
    if not DB:
        raise ValueError("Для алгоритма kraken2 требуется база данных (database)")

    # Шаг 1: Таксономическая классификация Kraken2
    rule kraken2_classify:
        input:
            r1 = INPUT_R1,
            r2 = INPUT_R2
        output:
            report = temp(f"{OUT_DIR}/kraken2_report.txt"),
            raw = temp(f"{OUT_DIR}/kraken2_output.txt")
        log:
            f"{OUT_DIR}/logs/kraken2.log"
        shell:
            """
            # Классификация ридов
            kraken2 --db {DB} --threads 8 --paired {input.r1} {input.r2} \
                    --output {output.raw} --report {output.report} > {log} 2>&1
            """

    # Шаг 2: Оценка обилия Bracken
    rule bracken_abundance:
        input:
            rules.kraken2_classify.output.report
        output:
            f"{OUT_DIR}/bracken_output.txt"
        log:
            f"{OUT_DIR}/logs/bracken.log"
        shell:
            """
            # Оценка обилия на уровне видов
            bracken -d {DB} -i {input} -o {output} -l S -t 10 > {log} 2>&1
            # Проверка и исправление файла
            awk 'NR == 1 || $4 ~ /^-?[0-9]+(\\.[0-9]+)?$/ {{print}}' {output} > {output}.tmp
            mv {output}.tmp {output}
            """


    # Шаг 3: Конвертация в формат Metaphlan
    rule convert_to_mpa:
        input:
            rules.bracken_abundance.output
        output:
            f"{OUT_DIR}/{SAMPLE}_report.tsv"
        shell:
            """
            # Генерация Metaphlan-подобного отчета
            kreport2mpa.py -r {input} -o {output} --display-header
            """


### Вариант 2: Metaphlan ###
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
            # Генерация таксономического профиля
            metaphlan {input.r1},{input.r2} --input_type fastq \
                      --nproc 8 -o {output} > {log} 2>&1
            """

else:
    raise ValueError("Неподдерживаемый алгоритм. Допустимые значения: kraken2, metaphlan")
