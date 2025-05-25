configfile: "config.yaml"


rule fastp_preprocessing:
    input: 
        fastq = config["input_fastq"]
    output: 
        processed_fastq = "output/processed.fastq.gz",
        html_report = "output/fastp_report.html",
        json_report = "output/fastp_report.json"
    threads: 6
    shell:
        """
        fastp -i {input.fastq} -o {output.processed_fastq} \
        -h {output.html_report} \
        -j {output.json_report} \
        --thread {threads} \
        --qualified_quality_phred 30 \
        --length_required 50
        """ 
        
# analysis with MetaPhlAn
rule analysis_MetaPhlAn:
    input: 
        processed_fastq = "output/processed.fastq.gz"
    output: 
        metaphlan = "output/metaphlan.tsv"
    threads: 8
    shell:
        """
        fastp -i {input.fastq} -o {output.processed_fastq} \
        -h {output.html_report} \
        -j {output.json_report} \
        --thread {threads} \
        --qualified_quality_phred 30 \
        --length_required 50
        """

# analysis withHUMAnN
rule humann_analysis:
    input: 
        processed_fastq = "output/processed.fastq.gz"
    output: 
        pathabundance = "output/humann_output/humann_pathabundance.tsv"
    threads: 8
    shell: 
        shell: 
        """
        humann {input.processed_fastq} \
        --input_type fastq \
        --output output/humann_output \
        --threads {threads}
        """
        
# rename
rule rename_human:
    input: 
        pathabundance = "output/humann_output/humann_pathabundance.tsv"
    output: 
        table_renamed = "output/humann_output/pathabundance_renamed.tsv"
    shell: 
        """
        humann_rename_table {input.pathabundance} \
        --output {output.table_renamed} \
        --names uniref90
        """

# group by
rule regroup_table_humman:
    input: 
        table_renamed = "output/humann_output/pathabundance_renamed.tsv"
    output: 
        regroup_table_humman =  "output/humann_output/pathabundance_regroup.tsv"
    shell:
        """
        humann_regroup_table {input.table_renamed} \
        --output {output.regroup_table_humman} \
        --groups metacyc
        """ 

# final rename 
rule final_humann_table:
    input: 
        regroup_table_humman =  "output/humann_output/pathabundance_regroup.tsv"
    output: 
        final_table = "output/humann_output/pathabundance_final.tsv"
    shell:
        """
        humann_rename_table {input.regroup_table_humman} \
        --output {output.final_table} \
        --groups metacyc
        """

# tmp file cleaning
rule clean_temporary_files:
    input: 
        processed_fastq = "output/processed.fastq.gz",
    output: 
        touch("output/cleanup_done.txt")
    shell: 
        """
        rm -f {input.processed_fastq}
        """


rule all:
    input:
        metaphlan_output = "output/metaphlan.tsv",
        humann_output = "output/humann_output/pathabundance_final.tsv"