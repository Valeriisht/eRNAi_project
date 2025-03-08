
####### Configuration #######
configfile: "enviromental.yaml"

SRA_ID = config.get("SRA_ID", "SRR123456") 
OUTPUT_DIR = "processed_data" 

# Создаем директорию для выходных файлов, если она не существует
shell(f"mkdir -p {OUTPUT_DIR}")


result_dir:  OUTPUT_DIR/        # should be kept
working_dir: temp/     

rule fastp_trim:
    input:
    output: 
    shell: 
    

