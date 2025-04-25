configfile: "config/config.yaml"

# Глобальные переменные
SRA_ID = config["sra"]["sra_id"]
OUT_DIR = config["output_dir"]
SAMPLE = config["sample_name"] 
TAXID = config["taxid"]
REFERENCE = config["ref_dir"]

# Правило для индексирования генома
rule bwa_index:
    input:
        host_genome = "{REFERENCE}/{TAXID}.fa".format(REFERENCE=REFERENCE, TAXID=TAXID)
    output: 
        expand("{reference}/{taxid}.fa.{ext}",
               reference=REFERENCE,
               taxid=TAXID,
               ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        threads = config["threads"]["bwa_index"]
    log:
        "{OUT_DIR}/logs/bwa_index.log".format(OUT_DIR=OUT_DIR)
    shell:
        "bwa index -p {input.host_genome} 2> {log}"

# Правило для выравнивания
rule bwa_align:
    input:
        r1 = "{OUT_DIR}/{SRA_ID}_filtered_1.fastq".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID),
        r2 = "{OUT_DIR}/{SRA_ID}_filtered_2.fastq".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID),
        reference_host_genome = "{REFERENCE}/{TAXID}.fa".format(REFERENCE=REFERENCE, TAXID=TAXID),
        index_files = rules.bwa_index.output
    output:
        alignment_file = "{OUT_DIR}/{SRA_ID}_aligned.sam".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID)
    params:
        threads = config["threads"]["bwa_align"],
        extra_params = "-v 1"
    log:
        "{OUT_DIR}/logs/{SRA_ID}_bwa_align.log".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID)
    shell:
        "bwa mem -t {params.threads} {params.extra_params} {input.reference_host_genome} {input.r1} {input.r2} > {output.alignment_file} 2> {log}"

# Конвертация и сортировка
rule sam_to_bam:
    input:
        sam = rules.bwa_align.output.alignment_file
    output:
        bam = "{OUT_DIR}/{SRA_ID}_aligned_sorted.bam".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID),
        bai = "{OUT_DIR}/{SRA_ID}_aligned_sorted.bam.bai".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID)
    params:
        threads = config["threads"]["sam_to_bam"],
        memory = "2G"
    log:
        "{OUT_DIR}/logs/{SRA_ID}_sam_to_bam.log".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID)
    shell:
        """
        samtools view -@ {params.threads} -b -o temp.bam {input.sam} 2>> {log}
        samtools sort -@ {params.threads} -m {params.memory} -o {output.bam} temp.bam 2>> {log}
        samtools index -@ {params.threads} {output.bam} 2>> {log}
        rm -f temp.bam
        """

# Разделение чтений
rule split_reads:
    input:
        bam = rules.sam_to_bam.output.bam,
        bai = rules.sam_to_bam.output.bai
    output:
        host_reads = "{OUT_DIR}/{SRA_ID}_host_reads.fastq.gz".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID),
        metagenome_reads = "{OUT_DIR}/{SRA_ID}_metagenome_reads.fastq.gz".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID)
    params:
        threads = config["threads"]["split_reads"]
    log:
        "{OUT_DIR}/logs/{SRA_ID}_split_reads.log".format(OUT_DIR=OUT_DIR, SRA_ID=SRA_ID)
    shell:
        """
        samtools fastq -@ {params.threads} -f 4 -1 {output.metagenome_reads} -2 /dev/null {input.bam} 2>> {log}
        samtools fastq -@ {params.threads} -F 4 -1 {output.host_reads} -2 /dev/null {input.bam} 2>> {log}
        """

