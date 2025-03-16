# Конфигурация
configfile: "config.yaml"

rule all:
    input:
        host_reads = expand("{output_dir}/host_reads.fastq.gz", output_dir=config["output_dir"]),
        metagenome_reads = expand("{output_dir}/metagenome_reads.fastq.gz", output_dir=config["output_dir"])


# индексирование генома хозяина 
rule bwa_index:
    input:
        host_genome = config["host_genome"] 
    output: 
        file_index = expand("{host_genome}.{ext}", host_genome=config["host_genome"], ext=["amb", "ann", "bwt", "pac", "sa"]) # файлы на выходе
    shell:
        """
        bwa index {input.host_genome} 
        """

rule bwa_align:
    input: 
        metagenome = config["metagenome"]
        refference_host_genome = config["host_genome"]  
        index_files = expand("{host_genome}.{ext}", host_genome=config["host_genome"], ext=["amb", "ann", "bwt", "pac", "sa"])
    threads: config["num_threads"]
    output:
        aligment_file = "aligned.sam"
    shell:
        """
        bwa mem -t {threads} {input.host_genome} {input.metagenome} > {output.aligment_file}
        """ 

# постпроцессинг 

# конвертация и сортинг
rule sam_to_bam:
    input:
        "aligned.sam"
    output:
        "aligned_sorted.bam"
    shell:
        """
        samtools view -b -o aligned.bam {input}
        samtools sort -o {output} aligned.bam
        samtools index {output}
        rm -f aligned.bam
        """

# разделение прочтений 
rule split_reads:
    input:
        bam="aligned_sorted.bam"
    output:
        host_reads=expand("{output_dir}/host_reads.fastq.gz", output_dir=config["output_dir"]),
        metagenome_reads=expand("{output_dir}/metagenome_reads.fastq.gz", output_dir=config["output_dir"])
    shell:
        """
        samtools view -b -F 4 {input.bam} > host_reads.bam
        samtools fastq host_reads.bam | gzip > {output.host_reads}
        samtools view -b -f 4 {input.bam} > metagenome_reads.bam
        samtools fastq metagenome_reads.bam | gzip > {output.metagenome_reads}
        rm -f host_reads.bam metagenome_reads.bam
        """








