# Конфигурация
configfile: "config.yaml"

SRA_ID = config["sra"]["sra_id"]
OUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"] 
INPUT_R1 = f"{OUT_DIR}/{SRA_ID}_filtered_1.fastq"
INPUT_R2 = f"{OUT_DIR}/{SRA_ID}_filtered_2.fastq"
TAXID = config["taxid"]
GENOME_FILE = f"{OUT_DIR}/{TAXID}.fna"

# индексирование генома хозяина 
rule bwa_index:
    input:
        host_genome = GENOME_FILE
    output: 
        file_index = f("{host_genome}.{ext}", host_genome=GENOME_FILE, ext=["amb", "ann", "bwt", "pac", "sa"]) # файлы на выходе
    params:
        threads = config["threads"]["bwa_index"]
    log:
        f"{OUT_DIR}/logs/bwa_index.log"
    shell:
        "bwa index -p {input.host_genome} --threads {params.threads} > {log} 2>&1"

rule bwa_align:
    input:
        r1 = INPUT_R1,
        r2 = INPUT_R2, 
        reference_host_genome = GENOME_FILE,  
        index_files = expand("{host_genome}.{ext}", host_genome=GENOME_FILE, ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        threads = config["threads"]["bwa_align"]
    output:
        aligment_file = f"{OUT_DIR}/aligned.sam"
    log:
        f"{OUT_DIR}/logs/bwa_align.log"
    shell:
        "bwa mem -t {params.threads} {input.reference_host_genome} {input.r1} {input.r2} > {output.aligment_file} 2> {log}"
        

# постпроцессинг 

# конвертация и сортинг
rule sam_to_bam:
    input:
        sam = rules.bwa_align.output
    output:
        bam = f("aligned_sorted.bam"),
        bai = f("aligned_sorted.bam.bai")
    params:
        threads = config["threads"]["sam_to_bam"]
        temp_bam = "aligned.tmp.bam"  # временный файл
    log:
        f"{OUT_DIR}/logs/sam_to_bam.log"
    shell:
        """
        samtools view -@ {params.threads} -b -o {params.temp_bam} {input.sam} 2>> {log}
        samtools sort -@ {params.threads} -o {output.bam} {params.temp_bam} 2>> {log}
        samtools index -@ {params.threads} {output.bam} 2>> {log}
        rm -f aligned.bam
        """

# разделение прочтений 
rule split_reads:
    input:
        bam = rules.sam_to_bam.output.bam
    output:
        host_reads=f("{OUT_DIR}/host_reads.fastq.gz"),
        metagenome_reads= f("{OUT_DIR}/metagenome_reads.fastq.gz")
    shell:
        """
        samtools view -b -F 4 {input.bam} > host_reads.bam
        samtools fastq host_reads.bam | gzip > {output.host_reads}
        samtools view -b -f 4 {input.bam} > metagenome_reads.bam
        samtools fastq metagenome_reads.bam | gzip > {output.metagenome_reads}
        rm -f host_reads.bam metagenome_reads.bam
        """






