TRANSCRIPTOME_FASTA = config['ref_dir_trans'] + "/{taxid}.fa"  
TRANS_SRA_ID = config["host_transcriptome_sra_id"]
INPUT_FASTQ_R1 = config["input_dir"] + "/{TRANS_SRA_ID}_filtered_1.fastq"  #  (R1)
INPUT_FASTQ_R2 = config["input_dir"] + "/{TRANS_SRA_ID}_filtered_2.fastq"  #  (R2) 
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
        index = OUTPUT_DIR + "/{taxid}_transcriptome.idx",
        r1 = INPUT_FASTQ_R1,
        r2 = INPUT_FASTQ_R2
    output:
        directory(OUTPUT_DIR + "/{taxid}/{TRANS_SRA_ID}_quant_results")
    params:
        bootstrap = config["kallisto"]["bootstrap"],
    log:
        "logs_kallisto/{taxid}/{TRANS_SRA_ID}_kallisto_quant.log"
    shell:
        """
        kallisto quant -i {input.index} -o {output} -b {params.bootstrap} \
            {input.r1} {input.r2} > {log} 2>&1
        """