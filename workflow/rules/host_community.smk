configfile: "config/config.yaml"

# Глобальные переменные
SRA_ID = config["sra"]["sra_id"]  # Список SRA ID
OUT_DIR = config["output_dir"]
INP_DIR = config["input_dir"]
TAXID = config["taxid"]
REFERENCE = config["ref_dir"]

# Создаем директории
os.makedirs(f"{OUT_DIR}/logs", exist_ok=True)

# Правило для индексирования (один раз)
rule bwa_index:
    input:
        host_genome = f"{REFERENCE}/{TAXID}.fa"
    output: 
        expand("{reference}/{taxid}.fa.{ext}",
               reference=REFERENCE,
               taxid=TAXID,
               ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        threads = config["threads"]["bwa_index"]
    log:
        f"{OUT_DIR}/logs/bwa_index.log"
    shell:
        "bwa index  {input.host_genome} 2> {log}"

# Правило для выравнивания (для каждого SRA ID)
rule bwa_align:
    input:
        r1 = f"{INP_DIR}/{{sra_id}}_filtered_1.fastq",
        r2 = f"{INP_DIR}/{{sra_id}}_filtered_2.fastq",
        reference_host_genome = f"{REFERENCE}/{TAXID}.fa",
        index_files = rules.bwa_index.output
    output:
        alignment_file = f"{OUT_DIR}/{{sra_id}}_aligned.sam"
    params:
        threads = config["threads"]["bwa_align"],
        extra_params = "-v 1"
    log:
        f"{OUT_DIR}/logs/{{sra_id}}_bwa_align.log"
    shell:
        "bwa mem -t {params.threads} {params.extra_params} {input.reference_host_genome} {input.r1} {input.r2} > {output.alignment_file} 2> {log}"

rule sam_to_bam:
    input:
        sam = f"{OUT_DIR}/{{sra_id}}_aligned.sam"
    output:
        bam = f"{OUT_DIR}/{{sra_id}}_aligned_sorted.bam",
        bai = f"{OUT_DIR}/{{sra_id}}_aligned_sorted.bam.bai"
    params:
        threads = 4,  # Используем 4 потока, как в рабочем тесте
        memory = "2G"  # Консервативный объем памяти
    log:
        f"{OUT_DIR}/logs/{{sra_id}}_sam_to_bam.log"
    shell:
        """
        # Шаг 1: Конвертация SAM → BAM (как в рабочем тесте)
        samtools view -@ {params.threads} -b -o {OUT_DIR}/temp_{wildcards.sra_id}.bam {input.sam} 2>> {log} || exit 1

        # Шаг 2: Сортировка BAM
        samtools sort -@ {params.threads} -m {params.memory} \
            -o {output.bam} {OUT_DIR}/temp_{wildcards.sra_id}.bam 2>> {log} || exit 1

        # Шаг 3: Индексация
        samtools index -@ {params.threads} {output.bam} 2>> {log} || exit 1

        # Очистка временного файла
        rm -f {OUT_DIR}/temp_{wildcards.sra_id}.bam
        """

rule split_reads:
    input:
        bam = f"{OUT_DIR}/{{sra_id}}_aligned_sorted.bam",
        bai = f"{OUT_DIR}/{{sra_id}}_aligned_sorted.bam.bai"
    output:
        host_reads = f"{OUT_DIR}/{{sra_id}}_host_reads.fastq.gz",
        metagenome_reads = f"{OUT_DIR}/{{sra_id}}_metagenome_reads.fastq.gz"
    params:
        threads = config["threads"]["split_reads"]
    log:
        f"{OUT_DIR}/logs/{{sra_id}}_split_reads.log"
    shell:
        """
        samtools fastq -@ {params.threads} -f 4 -1 {output.metagenome_reads} -2 /dev/null {input.bam} 2>> {log}
        samtools fastq -@ {params.threads} -F 4 -1 {output.host_reads} -2 /dev/null {input.bam} 2>> {log}
        """

