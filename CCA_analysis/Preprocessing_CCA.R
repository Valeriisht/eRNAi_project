## Preprocessing for CCA


##################################metagenome################################## 

# files loading 
file_path <- "kraken/"
file_pattern <- "*_species.report"  


# list of files

file_list <- list.files(path = file_path, 
                        pattern = file_pattern, 
                        full.names = TRUE)

process_metaphlan_file <- function(file_path) {
  
  sample_name <- gsub(".*_(.*)_bracken_species\\.report", "\\1", file_path)
  
  df <- read.delim(file_path, 
                   sep = "\t",
                   header = FALSE,
                   quote = "",  # abnormal name 
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


# root folder
root_path <- "kallisto_res/mrna/" # replace


#TSV-file 
pattern <- "*.tsv" 



# list - recursion 
trans_list <- list.files(path = root_path, 
                         pattern = pattern, 
                         full.names = TRUE,
                         recursive = TRUE) 

process_kallisto_file <- function(file_path) {
  # sample_name extraction
  sample_name <- gsub(".*/([^/]+)_quant_results/abundance\\.tsv", "\\1", file_path)
  
  df <- read.table(file_path, 
                   sep = "\t", 
                   header = TRUE,
                   stringsAsFactors = FALSE)
  
  df_processed <- df %>%
    select(target_id, tpm) %>%  #  TPM value
    rename(Gene = target_id, Abundance = tpm) %>%
    mutate(Sample = sample_name)
  
  return(df_processed)
}


all_transcriptome <- bind_rows(lapply(trans_list, process_kallisto_file))

# mtarix sample × genes
transcriptome_matrix <- all_transcriptome %>%
  pivot_wider(
    names_from = Sample,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  as.data.frame()

# set row names
rownames(transcriptome_matrix) <- transcriptome_matrix$Gene
transcriptome_matrix <- transcriptome_matrix[, -1]

#################################################################################################################################################################### 

# write.table(transcriptome_matrix, "transcriptome.tsv", sep="\t", row.names = TRUE)
save.image(file = "Preprocess_environment.RData")




