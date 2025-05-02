## подготовим файлы для анализа

# for file in *species.report; do
# echo "Обрабатываю $file..."
# /home/v.ishtuganova/eRNAi_project_last/eRNAi_project/kraken2/KrakenTools/kreport2mpa.py -r "$file" -o "${file%.*}.mpa" || echo "Ошибка в $file"
# done
# 

################################## metagenome

# файлы 
file_path <- "results/kraken_res/bracken/"
file_pattern <- "*_species.report"  


# Получаем список файлов

file_list <- list.files(path = file_path, 
                        pattern = file_pattern, 
                        full.names = TRUE)

process_metaphlan_file <- function(file_path) {
  
  sample_name <- gsub(".*/(.*).nt_bracken_species\\.report", "\\1", file_path)
  
  df <- read.delim(file_path, 
                   sep = "\t",
                   header = FALSE,
                   quote = "",  # игнорировать кавычки - аномальные назвния 
                   comment.char = "", 
                   col.names = c("Percentage", "Value_num", "Null", "Rank", "Taxonomy_ID", "Species"),
                   stringsAsFactors = FALSE)
  df_processed <- df %>% mutate(
    Sample = sample_name
  )
  
  return(df_processed)
}

all_data <- bind_rows(lapply(file_list, process_metaphlan_file))
all_data_filltered <- all_data[all_data$Rank == "S", ]

metagenome_filltered <- all_data_filltered[, c("Value_num", "Species", "Sample")]



####################################################################################################################################  


################################# transcriptome ################################################################################################### 


# Указываем путь к корневой папке
root_path <- "results/output/transcriptome_kallisto/mrna/" # во всех подпапках надо

# Находим все CSV-файлы в папках и подпапках
pattern <- "*.tsv" 



# Получаем список файлов (рекурсивно)
trans_list <- list.files(path = root_path, 
                         pattern = pattern, 
                         full.names = TRUE,
                         recursive = TRUE) 

process_kallisto_file <- function(file_path) {
  # Извлекаем имя образца из названия файла
  sample_name <- gsub(".*/([^/]+)_quant_results/abundance\\.tsv", "\\1", file_path)
  
  df <- read.table(file_path, 
                   sep = "\t", 
                   header = TRUE,
                   stringsAsFactors = FALSE)
  
  df_processed <- df %>%
    select(target_id, tpm) %>%  # Используем TPM значения
    rename(Gene = target_id, Abundance = tpm) %>%
    mutate(Sample = sample_name)
  
  return(df_processed)
}


all_transcriptome <- bind_rows(lapply(trans_list, process_kallisto_file))

# Преобразуем в матрицу образцов × генов
transcriptome_matrix <- all_transcriptome %>%
  pivot_wider(
    names_from = Sample,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  as.data.frame()

# Устанавливаем гены как названия строк
rownames(transcriptome_matrix) <- transcriptome_matrix$Gene
transcriptome_matrix <- transcriptome_matrix[, -1]

# write.table(transcriptome_matrix, "transcriptome.tsv", sep="\t", row.names = TRUE)

#################################################################################################################################################################### 


# save.image(file = "Preprocess_environment.RData")




