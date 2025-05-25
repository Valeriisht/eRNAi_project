TRANSCRIPTOME_FASTA = config['ref_dir'] + "/{taxid}.fa"  
SRA_ID = config["sra"]["sra_id"]
INPUT_FASTQ_R1 = config["input_dir"] + "/{SRA_ID}_filtered_1.fastq"  #  (R1)
INPUT_FASTQ_R2 = config["input_dir"] + "/{SRA_ID}_filtered_2.fastq"  #  (R2) 
OUTPUT_DIR =  config["output_dir"]  

# index transcriptome
rule kallisto_index:
    input:
        transcriptome = TRANSCRIPTOME_FASTA
    output:
        transcriptome_index = OUTPUT_DIR + "/{taxid}_transcriptome.idx"  #  {taxid}
    params:
        kmer_size = 31  # extra
    log:
        "logs_kallisto/{taxid}_kallisto_index.log"  # logs {taxid}
    shell:
        """
        kallisto index -i {output.transcriptome_index} -k {params.kmer_size} {input.transcriptome} > {log} 2>&1
        """

# pse-align
rule kallisto_quant:
    input:
        index = OUTPUT_DIR + "/{taxid}_transcriptome.idx",  #  {taxid}
        r1 = INPUT_FASTQ_R1,
        r2 = INPUT_FASTQ_R2
    output:
        directory(OUTPUT_DIR + "/{taxid}/{SRA_ID}_quant_results")  #  {taxid} and {SRA_ID}
    params:
        bootstrap = config["kallisto"]["bootstrap"],
    log:
        "logs_kallisto/{taxid}/{SRA_ID}_kallisto_quant.log"  # logs{taxid} and {SRA_ID}
    run:
        if os.path.exists(input.r2):
            shell(
                """
                kallisto quant -i {input.index} -o {output} -b {params.bootstrap} \
                {input.r1} {input.r2} > {log} 2>&1
                """
            )
        else:
            shell(
                """
                kallisto quant -i {input.index} -o {output} -b {params.bootstrap} {input.r1} > {log} 2>&1
                """
            )