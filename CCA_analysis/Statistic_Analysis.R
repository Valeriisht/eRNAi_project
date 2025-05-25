# loadings CCA

cca_result <- cca(transcriptome_reduced ~ ., data=metagenome_filtered)

species_loadings <- as.data.frame(scores(cca_result, display = "species", choices=1:2))
env_loadings <- as.data.frame(scores(cca_result, display = "bp", choices=1:2))



# dot product - 45 degrees 
is_within_45deg <- function(point, ref_vector) {
  cos_theta <- sum(point * ref_vector) / (sqrt(sum(point^2)) * sqrt(sum(ref_vector^2)))
  cos_theta >= cos(pi/4)
}

################################## 

collinear_results <- list()
against_collinear_results <- list()
non_colinear_results <- list()

# rename 
rownames(env_loadings) <- gsub("`", "", rownames(env_loadings))

rownames(env_loadings) <- trimws(rownames(env_loadings))

rownames(species_loadings) <- gsub("[^[:alnum:]]", "", rownames(species_loadings))


# cycle
for (i in 1:nrow(env_loadings)) {
  org_name <- rownames(env_loadings)[i]
  target_microbe <- as.data.frame(env_loadings[i, , drop=FALSE]) 
  v_norm <- target_microbe / sqrt(sum(target_microbe^2))
  
  # select
  selected <- apply(species_loadings, 1, is_within_45deg, ref_vector=v_norm)
  not_selected <- apply(species_loadings, 1, is_within_45deg, ref_vector=-v_norm)
  collinear_transcripts <- rownames(species_loadings)[selected]
  against_collinear_transcripts <- rownames(species_loadings)[not_selected]
  other_transcripts <- rownames(species_loadings)[!not_selected & !selected]
  
  # write 
  collinear_results[[org_name]] <- collinear_transcripts
  against_collinear_results[[org_name]] <- against_collinear_transcripts
  non_colinear_results[[org_name]] <- other_transcripts
  
  angle <- atan2(v_norm[1,2], v_norm[1,1])
  slope1 <- tan(angle + pi/4)
  slope2 <- tan(angle - pi/4)
  
  # vis
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + 
    
    geom_point(data=as.data.frame(species_loadings), aes(x = CCA1, y = CCA2), col="gray",
               alpha = 0.6) +
    geom_point(data=as.data.frame(species_loadings[selected, ]), aes(x = CCA1, y = CCA2), 
               col="#A47DAB", alpha = 0.8, size=5) +
    geom_point(data=as.data.frame(species_loadings[not_selected, ]), aes(x = CCA1, y = CCA2), 
               col="#FFE5B4", alpha = 0.8, size=5) + 
    geom_segment(aes(x = 0, y = 0, xend = v_norm[1,1], yend = v_norm[1,2]),
                 arrow = arrow(length = unit(0.2, "cm")), color = "aquamarine3") +
    geom_text(aes(x=target_microbe[1, 1], y=target_microbe[1, 2]), label=org_name,
              color = "red", hjust = 1.1, vjust = 1.1, fontface = "italic") +
    geom_abline(intercept = 0, slope = slope1, linetype = "dashed", color = "aquamarine3")+
    geom_abline(intercept = 0, slope = slope2, linetype = "dashed", color = "aquamarine3")+
    labs(title = substitute(paste("Collinear transcripts from ", italic(x)), list(x = org_name)),
         x = "CCA1", y = "CCA2") +
    theme_minimal()
  
  ggsave(filename = paste0("/nfs/home/emarenina/eRNAi_project/cca_image/plot_", org_name, ".png"), 
         plot = p, 
         width = 8, 
         height = 6)
  
}

#############################################################################################

# Statistic Analysis 

bslast_res <- list()

file_path <- "/nfs/home/emarenina/eRNAi_project/database_gallus/blast_res_1e-3/"

file_pattern <- "*.txt"  


# get a list of files
blast_list <- list.files(path = file_path, 
                         pattern = file_pattern, 
                         full.names = TRUE)


# func
process__file <- function(file_path) {
  
  sample_name <- sub("^(.*?)_reference.*$", "\\1", basename(file_path))
  
  df <- read.table(file_path, 
                   sep = "\t",
                   header = FALSE,
                   quote = "",  # ignore quotation marks - abnormal names 
                   comment.char = "", 
                   col.names = c("qseqid",  "sseqid", "pident",  "length",  "mismatch",  "gapopen",  "qstart",  "qend",  "sstart",  "send",    "evalue", "bitscore"),
                   stringsAsFactors = FALSE)
  df_processed <- df %>% mutate(
    Sample = sample_name
  )
  
  return(df_processed)
}

blast_data <- bind_rows(lapply(blast_list, process__file))


##### find the intersection ##### 

blast_id <- blast_data[, c("sseqid", "Sample")]
blast_id <- blast_id[blast_id$sseqid %in% rownames(species_loadings), ]


org <- unique(blast_id[, "Sample"])

blast_list <- list()

for (i in 1:length(org)) {
  blast_list[[org[i]]] <- unique(blast_id[blast_id$Sample==org[i],]$sseqid)
}

blast_names <- names(blast_list)



collinear_names <-names(collinear_results)
against_collinear_names <- names(against_collinear_results)
non_colinear_names <- names(non_colinear_results)


common_list_names <- intersect(blast_names, collinear_names)
common_list_names <- intersect(common_list_names, against_collinear_names)
common_list_names <- intersect(common_list_names, non_colinear_names)

blast_list_common <- blast_list[common_list_names]
collinear_results_common <- collinear_results[common_list_names]
against_collinear_results_common <- against_collinear_results[common_list_names]
non_colinear_results_common <- non_colinear_results[common_list_names]

# data 
org_data_list <- lapply(common_list_names, function(org) {
  list(
    similar = blast_list_common[[org]],
    collinear = collinear_results_common[[org]],
    against_collinear = against_collinear_results_common[[org]],
    non_colinear = non_colinear_results_common[[org]]
  )
}
)

names(org_data_list) <- common_list_names

# overlaps
calculate_overlaps <- function(org_data, func = "binomial") {
  overlap_collinear_similar <- intersect(org_data$similar, org_data$collinear)
  overlap_against_collinear_similar <- intersect(org_data$similar, org_data$against_collinear)
  overlap_non_collinear_similar <- intersect(org_data$similar, org_data$non_colinear)
  
  # group size 
  n_collinear <- length(org_data$collinear)
  n_against_collinear <- length(org_data$against_collinear)
  n_non_collinear <- length(org_data$non_colinear)
  n_similar <- length(org_data$similar)
  n_total <- length(unique(c(org_data$collinear, org_data$against_collinear, org_data$non_colinear)))
  
  p_collinear <- length(overlap_collinear_similar) / n_collinear
  p_against_collinear <- length(overlap_against_collinear_similar) / n_against_collinear
  p_non_collinear <- length(overlap_non_collinear_similar) / n_non_collinear
  
  if (func == "binomial") {
    p0 <- n_similar / n_total
    
    binom_collinear <- binom.test(
      x = length(overlap_collinear_similar),
      n = n_collinear,
      p = p0,
      alternative = "greater"
    )
    
    binom_against_collinear <- binom.test(
      x = length(overlap_against_collinear_similar),
      n = n_against_collinear,
      p = p0,
      alternative = "greater"
    )
    
    binom_non_collinear <- binom.test(
      x = length(overlap_non_collinear_similar),
      n = n_non_collinear,
      p = p0,
      alternative = "greater"
    )
    
    return(list(
      overlaps = list(
        collinear = length(overlap_collinear_similar),
        against = length(overlap_against_collinear_similar),
        non = length(overlap_non_collinear_similar)
      ),
      prop = c(
        collinear = p_collinear,
        against = p_against_collinear,
        non = p_non_collinear
      ),
      binomial_tests = list(
        collinear = binom_collinear,
        against = binom_against_collinear, 
        non = binom_non_collinear
      )
    ))
  }
  
  if (func == "fisher") {
    contingency_table <- matrix(
      c(
        length(overlap_collinear_similar), n_collinear - length(overlap_collinear_similar),
        length(overlap_against_collinear_similar), n_against_collinear - length(overlap_against_collinear_similar),
        length(overlap_non_collinear_similar), n_non_collinear - length(overlap_non_collinear_similar)
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("collinear", "against_collinear", "non_collinear"),
        c("similar", "not_similar")
      )
    )
    print(contingency_table)
    fisher_test <- fisher.test(contingency_table)
    
    return(list(
      overlap_counts = c(
        collinear_similar = length(overlap_collinear_similar),
        against_collinear_similar = length(overlap_against_collinear_similar),
        non_collinear_similar = length(overlap_non_collinear_similar)
      ),
      probabilities = c(
        p_collinear = p_collinear,
        p_against_collinear = p_against_collinear, 
        p_non_collinear = p_non_collinear
      ),
      fisher_test = fisher_test
    ))
  }
  
  stop("Unsupported function type. Use 'binomial' or 'fisher'.")
} 

overlap_stats <- lapply(org_data_list, calculate_overlaps)
  


###### vis ###### for each organism - binomial  ###### 
library(ggplot2)
library(gridExtra)

visualize_overlaps_one <- function(results, title = "Overlap Statistics") {
  
  overlap_data <- data.frame(
    Group = c("Collinear", "Against Collinear", "Non Collinear"),
    Count = c(results$overlaps$collinear,
              results$overlaps$against,
              results$overlaps$non),
    
    Proportion = c(results$prop["collinear"],
                   results$prop["against"],
                   results$prop["non"]),
    p_value = c(results$binomial_tests$collinear$p.value,
                results$binomial_tests$against$p.value,
                results$binomial_tests$non$p.value))
  
  overlap_data$p_label <- ifelse(overlap_data$p_value < 0.001, "***",
                                 ifelse(overlap_data$p_value < 0.01, "**",
                                        ifelse(overlap_data$p_value < 0.05, "*", "ns")))
  
  
  p1 <- ggplot(overlap_data,aes(x=Group, y=Count, fill = Group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Count, "\n", p_label)), vjust = -0.5, size = 4) +
    labs(title = paste(title, "- Counts"),
         x = "Group",
         y = "Overlapping Transcripts") +
    theme_minimal() +
    theme(legend.position = "none")
  
  p2 <- ggplot(overlap_data,aes(x=Group, y=Proportion, fill = Group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(Proportion, 3), "\n", p_label)), vjust = -0.5, size = 4) +
    labs(title = paste(title, "- Proportion"),
         x = "Group",
         y = "Proportion of Overlapping Transcripts") +
    ylim(0, 1) + 
    theme_minimal() +
    theme(legend.position = "none")
  
  grid.arrange(p1, p2, ncol = 2, top=title)
  
}
  
visualize_overlaps_one(overlap_stats$`Alistipes megaguti`)
  
#============
# general

visualize_overlaps <- function(results, title = "Overlap Statistics") {
  all_data <- data.frame()
  
  for (org in names(results)) {
    org_data <- results[[org]]
    
    overlap_data <- data.frame(
      Organism = org,
      Group = c("Co-directional", "Oppositely directed ", "Non Collinear"),
      Count = c(org_data$overlaps$collinear,
                org_data$overlaps$against,
                org_data$overlaps$non),
      
      Proportion = c(org_data$prop["collinear"],
                     org_data$prop["against"],
                     org_data$prop["non"]),
      p_value = c(org_data$binomial_tests$collinear$p.value,
                  org_data$binomial_tests$against$p.value,
                  org_data$binomial_tests$non$p.value))
    
    overlap_data$p_label <- ifelse(overlap_data$p_value < 0.001, "***",
                                   ifelse(overlap_data$p_value < 0.01, "**",
                                          ifelse(overlap_data$p_value < 0.05, "*", "")))
    
    all_data <- rbind(all_data, overlap_data) 
    
  }
  
  
  ggplot(all_data, aes(x = Organism %>% fct_reorder(Count), y = Count, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = paste0(p_label, "\n", Count)), 
              position = position_dodge(width = 0.9),
              vjust = -0.1, hjust = 0.5, size = 5, lineheight = 0.9, color = "gray30") +
    labs(title = "Overlapping Transcripts Across Microorganism by CCA, Binomial Test",
         x = "Microorganism",
         y = "Number of Overlapping Microorganism",
         fill = "Group Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 8),
          legend.position = "bottom",
          axis.title.y = element_text(margin = margin(r = 10)),
          ) +
    scale_fill_manual(values = c("#A47DAB", "#FFE5B4", "aquamarine3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(clip = "off")
}
  
(gg <- visualize_overlaps(overlap_stats))

ggsave("binomial_cca.png", plot = gg, device = "png", 
       width = 10, height = 8)


### fisher vis ### 
overlap_stats_fisher <- lapply(org_data_list, function(x) calculate_overlaps(x, "fisher"))

visualize_fisher <- function(results) {
  
  all_data <- data.frame()
  
  for (org in names(results)) {
    org_data <- results[[org]]
    
    overlap_data <- data.frame(
      Organism = rep(org, 3),
      Group = c("Co-directional", "Oppositely directed ", "Non Collinear"),
      Count = c(org_data$overlap_counts["collinear_similar"],
                org_data$overlap_counts["against_collinear_similar"],
                org_data$overlap_counts["non_collinear_similar"]),
      
      Probabilities = c(org_data$probabilities["p_collinear"],
                     org_data$probabilities["p_against_collinear"],
                     org_data$probabilities["p_non_collinear"]),
      
      
      p_value = rep(org_data$fisher_test$p.value %>% 
                      p.adjust(method = "fdr"), 3) )
     
    
    overlap_data$p_label <- ifelse(overlap_data$p_value < 0.001, "***",
                                   ifelse(overlap_data$p_value < 0.01, "**",
                                          ifelse(overlap_data$p_value < 0.05, "*", "")))
    
    all_data <- rbind(all_data, overlap_data) 
    
  }

  
  ggplot(all_data, aes(x = Organism %>% fct_reorder(Count), y = Count, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75)  +
    geom_text(aes(label = paste0(p_label, "\n", Count)), 
              position = position_dodge(width = 0.9),
              vjust = -0.1, hjust = 0.5, size = 5, lineheight = 0.9, color = "gray30") +
    labs(title = "Overlapping Transcripts Across Microorganism by CCA, Fisher Test",
         x = "Microorganism",
         y = "Number of Overlapping Microorganism",
         fill = "Group Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 8),
          legend.position = "bottom",
          axis.title.y = element_text(margin = margin(r = 10)),
    )+
    scale_fill_manual(values = c("#A47DAB", "#FFE5B4", "aquamarine3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(clip = "off")
  
} 

(fisher_cca <- visualize_fisher(overlap_stats_fisher))

ggsave("fisher_cca.png", plot = fisher_cca, device = "png", 
       width = 10, height = 8)


##==========================================
# statistic with CCA transcripts and not cca 

transcripts_cca <- rownames(species_loadings) 

transcripts <- read.table("/nfs/home/emarenina/eRNAi_project/data/data_ggallus/transcriptome/mrna.fa.ann")
  
transcripts_all <- transcripts$V2[
  grepl("[A-Za-z]", transcripts$V2) & 
    grepl("[0-9]", transcripts$V2)
]

blast_transcripts <- blast_data[, c("sseqid", "Sample")]

org_trans <- unique(blast_id[, "Sample"])

blast_list_transcripts <- list()

for (i in 1:length(org_trans)) {
  blast_list_transcripts[[org_trans[i]]] <- unique(blast_transcripts[blast_transcripts$Sample==org_trans[i],]$sseqid)
}


overlap_cca_transcriptome <-  function(table_blast, transcripts_cca, transcripts_all) {
  
  results <- list()
  
  for (org in names(table_blast)) {
    blast_hits <- table_blast[[org]]
      
    transcriptome_cca_target <- intersect(blast_hits, transcripts_cca)
    transcriptome_cca_other <- setdiff(transcripts_cca, blast_hits)
    transcripts_non_cca_target <- setdiff(blast_hits, transcripts_cca)
    transcripts_non_cca_other <- setdiff(setdiff(transcripts_all, transcripts_cca), blast_hits)
    
    fisher_table <-matrix(
      c(
        length(transcriptome_cca_target), 
        length(transcriptome_cca_other), 
        length(transcripts_non_cca_target), 
        length(transcripts_non_cca_other)
        ),
        nrow = 2,
      
      dimnames = list(
        "Transcript_type" = c("CCA", "non_CCA"), 
        "BLAST_result" = c("targets", "no_targets")     
      )
      )
    print(fisher_table)
    fisher_dif = (fisher_table[1,1]/fisher_table[1,2])/
      (fisher_table[2,1]/fisher_table[2,2])
    
    fisher_test <- fisher.test(fisher_table)

    
    results[[org]] <- list(
      counts = fisher_table,
      fisher_test = fisher_test,
      fisher_dif = fisher_dif,
      cca_targets = transcriptome_cca_target,
      cca_other = transcriptome_cca_other,
      non_cca_targets = transcripts_non_cca_target,
      non_cca_other = transcripts_non_cca_other
    )
    
  } 
  return(results)
}


cca_vs_non_cca <- overlap_cca_transcriptome(blast_list_transcripts, transcripts_cca, transcripts_all)

#### vis ####

p_values <- sapply(cca_vs_non_cca, function(x) x$fisher_test$p.value) %>% 
  p.adjust(method = "fdr")

fisher_dif <- sapply(cca_vs_non_cca, function(x) x$fisher_dif)

# barplot

df <- data.frame(
  organism = names(p_values),
  p_values = p_values,
  fisher_dif = fisher_dif,
  log_pvalue = -log10(p_values)
)

df <- df %>% 
  mutate(significance = case_when(
    log_pvalue > -log10(0.001) ~ "***",
    log_pvalue > -log10(0.05) ~ "**",
    log_pvalue > -log10(0.1) ~ ".",
    TRUE ~ " ",
  ))

g <- ggplot(df, aes(x=organism %>% fct_reorder(-log_pvalue), y=log_pvalue, fill = log_pvalue)) +
  geom_bar(stat="identity", width = 0.8)  +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  labs(title = "Fisher test adjusted p-values *(-log10 scale, fdr)*",
       y = "-log10 p-value",
       x = "Organism") +
  theme_minimal() +
  scale_fill_gradient("log p-value", low = "#FFE5B4", high = "#A47DAB") + 
  theme(axis.text.x = element_text(angle=80, hjust=1, face = "italic"),
        plot.title = ggtext::element_markdown()) +
  geom_text(aes(label = round(log_pvalue, 1)), 
            vjust = -0.5, size = 3.5, color = "black") 
  
print(g)


ggsave("CCA-out.png", plot = g, device = "png", 
       width = 15, height = 8)

#save.image(file = "Statistic_Analysis.RData")


#### within microorganism #### 

# Differences in probabilities of co-occurrence/anti-occurrence/out-of-area by group between species

org_data_list

prob_list <- lapply(org_data_list, prob_data_calc)
# table
prob_df <- do.call(rbind, prob_list)

prob_df$organism <- rownames(prob_df)

tbl <- tibble()
for (i in 1:length(org_data_list)) {
  tmp <- data.frame(
    name = names(org_data_list[i])
  )
  for (j in 1:4) {
    tmp[names(org_data_list[[i]][j])] <- length(org_data_list[[i]][j] %>% unlist)
  }
  tbl <- rbind(tmp, tbl)
}

tbl %>% 
  ggplot(aes(against_collinear/non_colinear, similar/collinear, color = name)) +
  geom_point() +
  theme_minimal() 

tbl %>% 
  ggplot(aes(against_collinear/similar, collinear/similar, color = name)) +
  geom_point() +
  theme_minimal() 

library(scales)

viridis_palette <- viridis_pal(option = "turbo")(13)

colnames(tbl) <- c("name", "DSK targets", "co-directional", "oppositely directional", "non-collinear")



library(ggradar)
ggradar(tbl,
        group.point.size = 0.5,
        group.line.width = 0.5,
        line.alpha = 0.3,
        axis.label.size = 3,
        legend.text = element_text(size = 8, face = "italic")) +
  scale_color_manual(values=viridis_palette) +
  theme_bw()

ggsave("windrose.png", width = 3000, height = 2000, units = "px")


library(rstatix)


tbl_long = tbl %>% pivot_longer(., 
                     cols=-name,
                     names_to = "direction_type", 
                     values_to = "count" )

res_ttest <- tbl_long %>%
  t_test(
  formula = count ~ name,
  p.adjust.method = "fdr"
  )

### vis ### 


library(ggpubr)
vis <- ggplot(tbl_long, aes(x=direction_type, y=count, color=name)) +
  geom_point(position=position_jitterdodge(jitter.width=0.2), size=3) +
  stat_summary(fun.data=mean_cl_normal, geom="errorbar", 
               width=0.3, position=position_dodge(0.8)) +
  stat_summary(fun=mean, geom="point", size=4, 
               position=position_dodge(0.8)) +
  scale_color_brewer(palette="Set1") +
  theme_minimal()


### common_transcript ### upset ### 

library(UpSetR)
library(tidyverse)
library(RColorBrewer)

colset <- brewer.pal(6, "Set3")

unique_transcripts <- unique(unlist(org_data_list))

# matrix
long_data <- imap_dfr(org_data_list, function(categories, org) {
  imap_dfr(categories, function(transcripts, category) {
    data.frame(
      org = org,
      category = category, 
      transcript = transcripts
    )
  })
})



binary_matrix <- long_data %>%
  mutate(present = 1) %>%
  unite("org_category", org, category, sep = "_") %>%
  pivot_wider(
    names_from = org_category,
    values_from = present,
    values_fill = 0,
    id_cols = transcript
  ) 


binary_matrix <- binary_matrix %>% 
  column_to_rownames("transcript") %>%  
  mutate(across(everything(), as.logical)) %>% t()



upset_data <- as.data.frame(binary_matrix) %>%
  rownames_to_column("transcript") %>%
  mutate(across(-transcript, as.logical)) %>%
  column_to_rownames("transcript")


gg <- ComplexUpset::upset(
  upset_data,
  intersect = colnames(upset_data)[1:20],
  name = "Transcripts",
  width_ratio = 0.2, 
  
  base_annotations = list(
    'Intersection size' = (
      ggplot2::ggplot() +
        ggplot2::geom_bar(
          stat = 'count',
          fill = "#A47DAB",
          width = 0.7 
        ) +
        ggplot2::geom_text(
          stat = 'count',
          aes(label = after_stat(count)),
          vjust = -0.5,
          hjust = -0.2,
          size = 3, 
          fontface = "bold"
        ) +
        ylab("Number of organism") +
        expand_limits(y = max(table(apply(upset_data, 1, paste, collapse = ""))) * 1.2) +
        theme(
          axis.title.y = element_text(size = 12), 
          axis.text.y = element_text(size = 12)
        )
    )
  ), 
  
  set_sizes = ComplexUpset::upset_set_size() +
    geom_text(
      aes(label = after_stat(count)),
      hjust = 1.1,
      stat = 'count',
      size = 5
      ) +
    ylab("Set size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
      
  themes = ComplexUpset::upset_modify_themes(
    list(
      "intersections_matrix" = theme(text = element_text(size = 12)),
      "Intersection size" = theme(axis.text.y = element_text(size = 12))
    )
  ),
  min_size = 1,
  sort_sets = "descending"
)   
  
print(gg)     

ggsave("Upset.png", plot = gg, device = "png", 
       width = 15, height = 8)


