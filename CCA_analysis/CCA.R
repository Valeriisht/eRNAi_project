# necessary library 

library(dplyr)
library(tidyr)
library(SpiecEasi)
library(limma)
library(tibble) 

# Multiline comment
# # Command + Shift + C
# Ctrl + Shift + C

### =========================================================================================================================
# Canonical Correspondence Analysis (CCA)
# Based on correspondent analysis, that is, a matrix 𝑸 is used instead of the original raw, centered data (𝒀)
# (A matrix of contributions to deviation from the null model based on independence between traits and objects).

# - CCA - canonical co-reference analysis in vegan
# - Dependent - samples 
# - Independent variables (responses) - abundance of species or transcriptome 
# - Environmental factors - species abundance or transcriptome 
### =========================================================================================================================

# Transpose tables for easy calculation of correlations

metagenome_raw <- metagenome_filltered %>% 
  
  group_by(Species) %>%
  filter(sum(Value_num == 0) / n() < 0.9) %>%
  ungroup() %>%
  
  pivot_wider(
    names_from = Species,
    values_from = Value_num,
    values_fill = 0  #  0 gaps pivoting
  ) %>% 
  
  column_to_rownames(var = names(.)[1]) %>% 
  
  .[, colSums(. > 0) >= max(3, 0.25 * nrow(.))] %>% # min 3 on 25% species 
  
  
  .[, colSums(.) >= 0.0005 * sum(.)] %>%
  
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  apply(1, function(x) x / sum(x)) 

metagenome <- metagenome_raw %>% 
  #decostand("hellinger") %>% 
  .[!rownames(.) %in% c("Mus Musculus", "Gallus gallus", "Meleagris gallopavo"), ]




### =========================================================================================================================

# Delete rows (organisms or transcripts) where more than 90% of the values are zero
# For metagenome - above

# For transcriptome (if data are not normalized)
transcriptome <- transcriptome_matrix[rowSums(transcriptome_matrix == 0) / ncol(transcriptome_matrix) < 0.9, ]


### =========================================================================================================================

# transcriptome <- data.frame(lapply(transcriptome[, -1], function(x) as.numeric(as.character(x))))
# transcriptome <- sweep(transcriptome, 2, colSums(transcriptome), "/")


# normaization
metagenome <- sweep(metagenome, 2, colSums(metagenome, na.rm = TRUE), "/")
transcriptome <- sweep(transcriptome, 2, colSums(transcriptome), "/")


metagenome <- t(metagenome)
transcriptome <- t(transcriptome)


####### string override - matching 

desired_order <- read.table("mapping_table.tsv", sep="\t", header=TRUE)
desired_order <- desired_order[desired_order$meta %in% rownames(metagenome), ]

metagenome <- metagenome[desired_order$meta, ]
transcriptome <- transcriptome[desired_order$trans, ]

rownames(transcriptome) <- rownames(metagenome)

rownames(metagenome) == rownames(transcriptome)


### =========================================================================================================================

# How does SPARCC work?
# 
# - Logarithm of ratios - instead of absolute values, logarithms of ratios between components are used.
# 
# - Iterative estimation of correlations - weak correlations are eliminated at each step.
# 
# - Bootstrap validation - assessing the significance of correlations through random subsamples.
# 

### =========================================================================================================================

# Function for calculating correlations


calculate_correlations <- function(metagenome, transcriptome, method = "sparcc") {
  if (method == "sparcc") {
    # Use SparCC to calculate correlations
    sparcc_results <- SpiecEasi::sparcc(metagenome, transcriptome)
    correlations <- sparcc_results$cor
    p_values <- sparcc_results$pval
  } else {
    # use the standard cor function
    correlations <- cor(metagenome, transcriptome, method = method)
    p_values <- matrix(NA, nrow = ncol(metagenome), ncol = ncol(transcriptome))  # p-значения не вычисляются в cor
  }
  return(list(correlations = correlations, p_values = p_values))
}

#corr_res <- calculate_correlations(metagenome_reduced, transcriptome_reduced)

### =========================================================================================================================

library(vegan) 

# Deleting low-expressed genes

keep_genes <- colSums(transcriptome > 0) >= 0.5 * nrow(transcriptome)
transcriptome_filt <- transcriptome[, keep_genes]



# Selection of top-variable genes (e.g., 500 most variable genes)
top_genes <- names(sort(apply(transcriptome_filt, 2, sd), decreasing = TRUE)[1:200])


transcriptome_reduced <- transcriptome_filt[, top_genes] %>% 
  as.data.frame() %>%
  {rownames(.) <- rownames(transcriptome_filt); .} %>%  # save rownames
  apply(2, as.numeric) %>% 
  as.data.frame() %>%
  {rownames(.) <- rownames(transcriptome_filt); .}  


##################### CCA ###########################

reduction <- function(df, mincols = 3) {
  df %>%
    subset(apply(., 1, function(x) sum(x > 0, na.rm = TRUE) >= mincols)) # removing degenerate organisms
}

metagenome_reduced <- as.data.frame(reduction(metagenome))


metagenome_filtered <- metagenome_reduced %>%
  .[, colSums(. > 0) >= 0.10 * nrow(.)] %>% #  by prevalence (10% of samples)
  .[, colSums(.) >= 0.0001 * sum(.)]   #  by abundance (0.01% of total reads)

transcriptome_reduced <-  transcriptome_reduced[!rownames(transcriptome_reduced) %in% "SRR8265551", ]

############################################################
# first option 

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
  mutate(
    value = apply(metagenome_raw, 1, mean)[rownames(df_species)],  # sample average
    abs_value = abs(value)  #  significance for importance
  ) %>% 
  arrange(-abs_value) %>%  #  descending |value|
  mutate(label = NA)
df_species$label <- rownames(df_species)



df_env <- as.data.frame(env * 0.8) 


df_env <- df_env %>% 
  mutate(value = apply(transcriptome, 2, mean)[rownames(df_env)]) %>% 
  arrange(-value) %>% 
  rownames_to_column("label")

df_species$label[-(1:10)] <- NA  # only the top 10 tags
df_env$label[-(1:5)] <- NA

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

g1 <- g1 + theme(aspect.ratio = 0.5)

g11 <- g1 + geom_text_repel(data = df_species, aes(x = CCA1, y = CCA2, label = label), 
                            color = c("#51107A"), fontface = "italic", segment.color = NA)
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
  mutate(value = apply(metagenome_raw, 1, mean)[rownames(df_species) %>% str_remove_all("`")]) %>% 
  arrange(-value) %>% 
  rownames_to_column("label") %>% 
  mutate(label = label %>% str_remove_all("`"))


df_env <- as.data.frame(env * 0.8) 
df_env <- df_env %>% 
  mutate(value = apply(transcriptome, 2, mean)[rownames(df_env)]) %>% 
  arrange(-value) %>% 
  mutate(label = NA)


df_env$label <- rownames(df_env)

df_species$label[-(1:10)] <- NA  #  only the top 10 tags
df_env$label[-(1:5)] <- NA


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
       subtitle = "Linkage of metagenome and transcriptome data") + 
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom")

g2 <- g2 + theme(aspect.ratio = 0.5)

g21 <- g2 + geom_text_repel(data = df_env, aes(x = CCA1, y = CCA2, label = label), 
                            color = c("#51107A"), fontface = "bold", segment.color = NA)
g22 <- g2 +  geom_text_repel(data = df_species, aes(x = CCA1, y = CCA2, label = label),
                             color = "#115557", size = 4, fontface = "italic") 
gg2 <- ggarrange(g21+xlab(""), g22 +labs(title = NULL, subtitle = NULL), ncol = 1, common.legend = T, legend = "bottom", align = "h")

g_res <- ggpubr::ggarrange(g11, g21, align = "hv")

ggsave("CCA_final.png", plot = g_res, device = "png", 
       width = 10, height = 8)

############################################################
# load extraction 

species_loadings <- scores(cca_result, display = "species", choices=1:2)
env_loadings <- scores(cca_result, display = "bp", choices=1:2)

write.csv(species_loadings, "species_loadings.csv")
write.csv(env_loadings, "env_loadings.csv")

############################################################

save.image(file = "CCA_environment.RData")

############################


