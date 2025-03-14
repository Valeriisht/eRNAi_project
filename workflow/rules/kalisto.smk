
TRANSCRIPTOME_FASTA = "transcriptome.fa"  # Файл с транскриптомом
INPUT_FASTQ_R1 = "reads_1.fastq"          # Файл с ридами 
INPUT_FASTQ_R2 = "reads_2.fastq"          # Файл с ридами - опционально 
OUTPUT_DIR = "transcriptome"              # output directory


shell:
    """
    mkdir -p {OUTPUT_DIR}
    mkdir -p logs
    """
    

# индексирование генома 

rule kallisto_index: 
    input: 
        transcriptome = TRANSCRIPTOME_FASTA
    output: 
        transcriptome_index = OUTPUT_DIR + "/transcriptome.idx" 

    params:
        kmer_size = 31 # опционально 
    
    log:
        "logs/kallisto_index.log"
        
    shell:
        """
        
        kallisto index -i {output.transcriptome_index} -k {params.kmer_size} {input.transcriptome} > {log} 2>&1

        """

# квази-выравнивание 

rule kallisto_quant:

    input: 
        index = OUTPUT_DIR + "/transcriptome.idx"
        r1 = INPUT_FASTQ_R1 
        r2 = INPUT_FASTQ_R2 
    output: 
        directory(OUTPUT_DIR + "/quant_results")  
    params: 
        bootstrap = 100
    log:
        "logs/kallisto_quant_.log"
    run:
        if os.path.exists(input.r2):
            shell(
                """

                kallisto quant -i {input.index} -o {output} -b {params.bootstrap} {input.r1} {input.r2} > {log} 2>&1
                
                """
            )
        else:
            shell (
                """

                kallisto quant -i {input.index} -o {output} -b {params.bootstrap} {input.r1} > {log} 2>&1

                """
            )
