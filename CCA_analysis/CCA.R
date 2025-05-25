# Загрузка необходимых библиотек - не устанавоиваем 
library(dplyr)
library(tidyr)
library(SpiecEasi)
library(limma)
library(tibble) 

# Многострочный комментарий
# # Command + Shift + C
# Ctrl + Shift + C

### =========================================================================================================================
# Канонический корреспондентный анализ (CCA)
# Основан на корреспондентном анализе, то есть вместо исходных сырых, центрированных данных (𝒀) используется матрица 𝑸
# (Матрица вкладов в отклонение от нулевой модели, основанной на независимости между признаками и объектами).

 # - CCA - канонический кореспондентный анализ в vegan
 #   - Зависимые переменные (отклики) - обилие видов
 #   - Независимые переменные (предикторы) - транскриптом ???
 #  


# Транспонирование таблиц для удобства вычисления корреляций

metagenome_raw <- metagenome_filltered %>% 
  
  group_by(Species) %>%
  filter(sum(Value_num == 0) / n() < 0.9) %>%
  ungroup() %>%
  
  pivot_wider(
    names_from = Species,
    values_from = Value_num,
    values_fill = 0  #  пропуски нулями сразу при pivoting
  ) %>% 
  
  column_to_rownames(var = names(.)[1]) %>% 
  
  .[, colSums(. > 0) >= max(3, 0.25 * nrow(.))] %>% # минимум 3  образцов видов 25% видов
  
  
  .[, colSums(.) >= 0.0005 * sum(.)] %>%
  
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  apply(1, function(x) x / sum(x)) 

metagenome <- metagenome_raw %>% 
  decostand("hellinger") %>% 
  .[!rownames(.) %in% c("Mus Musculus", "Gallus gallus", "Meleagris gallopavo"), ] #долой эукариот 




### =========================================================================================================================

# Удаляем строки (организмы или транскрипты), где более 90% значений равны нулю
# Для метагенома - выше

# Для транскриптома (если данные не нормализованы)
transcriptome <- transcriptome_matrix[rowSums(transcriptome_matrix == 0) / ncol(transcriptome_matrix) < 0.9, ]

# плавим таблицу траснкриптома 


### =========================================================================================================================
# нормализация с помощью программы lima  


# transcriptome <- data.frame(lapply(transcriptome[, -1], function(x) as.numeric(as.character(x))))
# transcriptome <- sweep(transcriptome, 2, colSums(transcriptome), "/")


# норм
metagenome <- sweep(metagenome, 2, colSums(metagenome, na.rm = TRUE), "/")
transcriptome <- sweep(transcriptome, 2, colSums(transcriptome), "/")


metagenome <- t(metagenome)
transcriptome <- t(transcriptome)


####### переопределение строк - соответсвие 

desired_order <- c("SRR8265549", "SRR8265548", "SRR8265622", "SRR8265550",
                   "SRR8265657", "SRR8265655", "SRR8265631", "SRR8265632",
                   "SRR8265651", "SRR8265658", "SRR8265652", "SRR8265594",
                   "SRR8265621", "SRR8265629")

transcriptome <- transcriptome[desired_order, ]
transcriptome <- transcriptome[!rownames(transcriptome) %in% "SRR8265631", ]



rownames(metagenome) <- c("SRR8265549", "SRR8265548", "SRR8265622", "SRR8265550",
                          "SRR8265657", "SRR8265655", "SRR8265632",
                          "SRR8265651", "SRR8265658", "SRR8265652", "SRR8265594",
                          "SRR8265621", "SRR8265629")

rownames(metagenome) == rownames(transcriptome)


# transcriptome <- data.frame(lapply(transcriptome, function(x) as.numeric(as.character(x))))

# T_reads.cca = cca(transcriptome~., metagenome, dist="bray", direction = 'forward', permutations = 999)


### =========================================================================================================================


# Как работает SPARCC?
# 
# - Логарифмирование соотношений — вместо абсолютных значений используются логарифмы отношений между компонентами.
# 
# - Итеративная оценка корреляций — на каждом шаге исключаются слабые связи.
# 
# - Бутстреп-проверка — оценка значимости корреляций через случайные подвыборки.
# 

### =========================================================================================================================

# Функция для вычисления корреляций


calculate_correlations <- function(metagenome, transcriptome, method = "sparcc") {
  if (method == "sparcc") {
    # Используем SparCC для вычисления корреляций
    sparcc_results <- SpiecEasi::sparcc(metagenome, transcriptome)
    correlations <- sparcc_results$cor
    p_values <- sparcc_results$pval
  } else {
    # Используем стандартную функцию cor
    correlations <- cor(metagenome, transcriptome, method = method)
    p_values <- matrix(NA, nrow = ncol(metagenome), ncol = ncol(transcriptome))  # p-значения не вычисляются в cor
  }
  return(list(correlations = correlations, p_values = p_values))
}

#corr_res <- calculate_correlations(metagenome_reduced, transcriptome_reduced)

# сопоставление - у нас 154 таксона - вывберем гены по экспрессии

library(vegan) 

# Удаление низкоэкспрессируемых генов
keep_genes <- rowSums(transcriptome > 0) >= 0.5*ncol(transcriptome) 

transcriptome_filt <- transcriptome[, keep_genes]

# 
# 
# Выбор топ-вариабельных генов (например, 500 самых изменчивых)
top_genes <- names(sort(apply(transcriptome_filt, 2, sd), decreasing = TRUE)[1:100])
# 
# 
# 

transcriptome_reduced <- transcriptome_filt[, top_genes] %>% 
  as.data.frame() %>%
  {rownames(.) <- rownames(transcriptome_filt); .} %>%  # Сохраняем rownames
  apply(2, as.numeric) %>% 
  as.data.frame() %>%
  {rownames(.) <- rownames(transcriptome_filt); .}  


# correlations <- calculate_correlations(metagenome, transcriptome_reduced)

##################### CCA ###########################

reduction <- function(df, mincols = 3) {
  df %>%
    subset(apply(., 1, function(x) sum(x > 0, na.rm = TRUE) >= mincols)) # убираем вырожденные организмы
}

metagenome_reduced <- as.data.frame(reduction(metagenome))


metagenome_filtered <- metagenome_reduced %>%
  .[, colSums(. > 0) >= 0.10 * nrow(.)] %>% #  по распространённости (10% образцов)
  .[, colSums(.) >= 0.0001 * sum(.)]   #  по abundance (0.01% от общего числа reads)

# число нулей 

summary(colSums(metagenome_filtered == 0))

############################################################
# первый вариант 

library(tidyverse)
library(ggrepel)
library(ggpubr)

# 1
cca_result <- cca(metagenome_filtered ~ ., data=transcriptome_reduced)

#data manipulation
sites <- scores(cca_result, display = "sites")
species <- scores(cca_result, display = "species")
env <- scores(cca_result, display = "bp")

df_sites <- as.data.frame(sites)
df_species <- as.data.frame(species) 

#add info
df_species <- df_species %>% 
  mutate(value = apply(metagenome_raw, 2, mean)[rownames(df_species)]) %>% 
  arrange(-value) %>% 
  mutate(label = NA)
df_species$label[1:10] <- rownames(df_species)[1:10]


df_env <- as.data.frame(env * 0.8) 

df_env <- df_env %>% 
  mutate(value = apply(transcriptome, 2, mean)[rownames(df_env)]) %>% 
  arrange(-value) %>% 
  rownames_to_column("label")

# plots

g1 <- ggplot() +
  geom_point(data=df_sites, aes(x=CCA1, y=CCA2), 
             color = "darkblue", size = 2, alpha = 0.7) +
  
  #geom_text_repel(data=df_sites, aes(x=CCA1, y=CCA2), label = rownames(df_sites), color = "steelblue", size = 3) +
  
  geom_point(data = df_species, aes(x = CCA1, y = CCA2, color = -log10(value)),
             shape = 17, size = 3) + 
  scale_color_gradientn("Metagenome,\nmean value, -lg", colors = c("lightgreen","pink", "purple"), n.breaks = 3)+
  ggnewscale::new_scale_color() +
  
  geom_segment(data = df_env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = -log10(value)),
               arrow = arrow(
                 length = unit(0.3, "cm"), 
                 angle = 20,  
                 type = "closed"
               )) +
  
  scale_color_gradientn("Transcriptome,\nmean value, -lg", colors = c("#F09590", "yellow","cadetblue4"), n.breaks = 3)+
  
  theme_minimal() +
  
  labs(x = paste0("CA1: ", round(cca_result$CCA$eig[1]/sum(cca_result$CCA$eig)*100, 1), "%"),
       y = paste0("CA2:", round(cca_result$CCA$eig[2]/sum(cca_result$CCA$eig)*100, 1), "%"),
       title = "CCA",
       subtitle = "Linkage of transcriptome and metagenome data") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom")

g11 <- g1 + geom_text_repel(data = df_species, aes(x = CCA1, y = CCA2, label = label), 
                            color = c("#51107A"), fontface = "italic")
g12 <- g1 + geom_text_repel(data = df_env, aes(x = CCA1, y = CCA2, label = label),
                            color = "#115557", size = 4, fontface = "bold")

gg1 <- ggarrange(g11+xlab(""), g12 +labs(title = NULL, subtitle = NULL), ncol = 1, common.legend = T, legend = "bottom", align = "h")

# 2
cca_result <- cca(transcriptome_reduced ~ ., data=metagenome_filtered)

#data manipulation
sites <- scores(cca_result, display = "sites")
species <- scores(cca_result, display = "bp")
env <- scores(cca_result, display = "species")

df_sites <- as.data.frame(sites)
df_species <- as.data.frame(species) 

#add info
df_species <- df_species %>% 
  mutate(value = apply(metagenome_raw, 2, mean)[rownames(df_species) %>% str_remove_all("`")]) %>% 
  arrange(-value) %>% 
  rownames_to_column("label") %>% 
  mutate(label = label %>% str_remove_all("`"))

df_env <- as.data.frame(env * 0.8) 
df_env <- df_env %>% 
  mutate(value = apply(transcriptome, 2, mean)[rownames(df_env)]) %>% 
  arrange(-value) %>% 
  mutate(label = NA)

df_env$label[1:10] <- rownames(df_env)[1:10]

# plots

g2 <- ggplot() +
  geom_point(data=df_sites, aes(x=CCA1, y=CCA2), 
             color = "darkblue", size = 2, alpha = 0.7) +
  
  #geom_text_repel(data=df_sites, aes(x=CCA1, y=CCA2), label = rownames(df_sites), color = "steelblue", size = 3) +
  
  geom_point(data = df_env, aes(x = CCA1, y = CCA2, color = -log10(value)),
             shape = 17, size = 3) + 
  scale_color_gradientn("Transcriptome,\nmean value, -lg", colors = c("#F09590","yellow", "cadetblue4"), n.breaks = 3)+
  ggnewscale::new_scale_color() +
  
  geom_segment(data = df_species, aes(x = 0, y = 0, xend = CCA1, yend = CCA2, color = -log10(value)),
               arrow = arrow(
                 length = unit(0.3, "cm"), 
                 angle = 20,  
                 type = "closed"
               )) +
  
  
  scale_color_gradientn("Metagenome,\nmean value, -lg", colors = c("lightgreen","pink", "purple"), n.breaks = 3)+
  
  theme_minimal() +
  
  labs(x = paste0("CA1: ", round(cca_result$CCA$eig[1]/sum(cca_result$CCA$eig)*100, 1), "%"),
       y = paste0("CA2:", round(cca_result$CCA$eig[2]/sum(cca_result$CCA$eig)*100, 1), "%"),
       title = "CCA",
       subtitle = "Linkage of metagemome and transcriptome data") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom")

g21 <- g2 + geom_text_repel(data = df_env, aes(x = CCA1, y = CCA2, label = label), 
                            color = c("#51107A"), fontface = "bold")
g22 <- g2 +  geom_text_repel(data = df_species, aes(x = CCA1, y = CCA2, label = label),
                             color = "#115557", size = 4, fontface = "italic") 

gg2 <- ggarrange(g21+xlab(""), g22 +labs(title = NULL, subtitle = NULL), ncol = 1, common.legend = T, legend = "bottom", align = "h")

ggpubr::ggarrange(gg1, gg2, align = "hv")


############################################################
# извлечение нагрузок 

species_loadings <- scores(cca_result, display = "species", choices=1:2)
env_loadings <- scores(cca_result, display = "bp", choices=1:2)

write.csv(species_loadings, "species_loadings.csv")
write.csv(env_loadings, "env_loadings.csv")

############################################################
# корреляция 


library(pheatmap)

loadings_matrix <- rbind(cca_result$CCA$v[, 1:2], 
                         cca_result$CCA$biplot[, 1:2])

rownames(loadings_matrix) <- c(colnames(metagenome_filtered), colnames(transcriptome_reduced))


pheatmap(loadings_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("blue", "white", "brown"))(50),
         main = "Main HeatMap"
         )


correlations <- calculate_correlations(metagenome_filtered[, 1:80], transcriptome_reduced[, 1:80])


library(corrplot)
corrplot(correlations[1:10, 1:10])


vector_genera <- colnames(metagenome_reduced)



# write.table(metagenome_reduced, "metagenome_cca.tsv", sep="\t", row.names = TRUE)
# write.table(transcriptome_reduced, "transcriptome_cca.tsv", sep="\t", row.names = TRUE)

save.image(file = "CCA_environment.RData")

############################
# старый код для построения графика 

# (plot <- ggplot() +
#    geom_point(data=df_sites, aes(x=CCA1, y=CCA2), 
#               color = "steelblue", size = 3, alpha = 0.7) +
#    
#    geom_text_repel(data=df_sites, aes(x=CCA1, y=CCA2), label = rownames(df_sites), 
#                    color = "steelblue", size = 3) +
#    
#    geom_point(data = df_species, aes(x = CCA1, y = CCA2),
#               color = "firebrick", shape = 17, size = 3) + 
#    
#    geom_segment(data = df_env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
#                 arrow = arrow(
#                   length = unit(0.3, "cm"), 
#                   angle = 20,              
#                 ), color = "purple") +
#    
#    geom_text_repel(data = df_env, aes(x = CCA1, y = CCA2, label = rownames(df_env)),
#                    color = "purple", size = 4, fontface = "bold") +
#    
#    theme_minimal() +
#    
#    labs(x = paste("CCA1 (", round(cca_result$CCA$eig[1]/sum(cca_result$CCA$eig)*100, 1), "%)"),
#         y = paste("CCA2 (", round(cca_result$CCA$eig[2]/sum(cca_result$CCA$eig)*100, 1), "%)"),
#         title = "CCA",
#         subtitle = "Linkage of transcriptome and metagemome data") + coord_fixed(ratio = 1))
# 

# cca_result <- cca(transcriptome_reduced ~ ., data=metagenome_filtered)
# summary(cca_result)
# 
# plot(cca_result)
# 
# # оценка значимости 
# 
# perm_test <- anova.cca(cca_result, permutations = 9999)
# print(perm_test)
# 
# 
# sites <- scores(cca_result, display = "sites")
# species <- scores(cca_result, display = "species")
# env <- scores(cca_result, display = "bp")
# 
# 
# # в датафрейм
# df_sites <- as.data.frame(sites)
# df_species <- as.data.frame(species)
# df_env <- as.data.frame(env * 0.8) 
# 
# library(ggplot2)
# library(ggrepel) 


