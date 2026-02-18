library(MotrpacHumanPreSuspension)
library(TMSig)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(purrr)
library(tidyr)
library(ggh4x)
library(stringr)
library(RColorBrewer)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(readxl)

motrpac_qc <- load_qc()
motrpac_da <- load_differential_analysis()

# This helps to simplify naming: For transcriptomics and metbabolomics
subset_list <- list(
  RE_4hr = "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise",
  EE_4hr = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise",
  Control_4hr = "group_timepointADUControl.post_3.5_4_hr - group_timepointADUControl.pre_exercise",
  EE_45min = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise",
  RE_45min = "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise",
  EE_45min_vs_Control = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  RE_45min_vs_Control = "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  EE_45min_vs_RE = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise",
  EE_24hr = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise",
  RE_24hr = "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise",
  Control_45min = "group_timepointADUControl.post_15_30_45_min - group_timepointADUControl.pre_exercise",
  RE_4hr_vs_Control = "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  EE_4hr_vs_Control = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  Control_24hr = "group_timepointADUControl.post_24_hr - group_timepointADUControl.pre_exercise",
  RE_24hr_vs_Control = "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  EE_24hr_vs_Control = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  EE_24hr_vs_RE = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_24_hr + group_timepointADUResist.pre_exercise",
  RE_vs_Control = "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise",
  EE_vs_RE = "group_timepointADUEndur.pre_exercise - group_timepointADUResist.pre_exercise",
  EE_vs_Control = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise",
  EE_4hr_vs_RE = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_3.5_4_hr + group_timepointADUResist.pre_exercise"
)
# For proteomics and phosphoproteomics
prot_list <- list(
  RE_4hr = "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise",
  EE_4hr = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise",
  Control_4hr = "group_timepointADUControl.post_3.5_4_hr - group_timepointADUControl.pre_exercise",
  RE_4hr_vs_Control = "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  EE_4hr_vs_Control = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  RE_vs_Control = "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise",
  EE_vs_RE = "group_timepointADUEndur.pre_exercise - group_timepointADUResist.pre_exercise",
  EE_vs_Control = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise",
  EE_4hr_vs_RE = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_3.5_4_hr + group_timepointADUResist.pre_exercise"
)

# 1. Transcriptomics
precawg_trans_da <- as.data.frame(motrpac_da$adipose$`transcript-rna-seq`)
precawg_trans_da <- precawg_trans_da %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")

t_vol <- map_df(
  names(subset_list),
  function(name) {
    contrast_value <- subset_list[[name]]
    
    precawg_trans_da %>%
      filter(contrast == contrast_value) %>%
      mutate(
        facet_label = gsub("_", " ", name),
        is_ee_vs_re = grepl("EE_vs_RE", name),
        color = case_when(
          adj_p_value >= 0.05 ~ "grey",
          is_ee_vs_re & logFC > 0 ~ "#d95f02",
          is_ee_vs_re & logFC < 0 ~ "#1b9e77",
          grepl("EE", name) ~ "#d95f02",
          grepl("RE", name) ~ "#1b9e77",
          TRUE ~ "#7570b3"
        )
      ) %>%
      group_by(facet_label) %>%
      mutate(
        gene_symbol = as.character(gene_symbol),
        p_rank = rank(p_value),
        label = ifelse(adj_p_value < 0.05 & p_rank <= 10, gene_symbol, NA)
      ) %>%
      ungroup()
  }
)
##Transriptomics DA PCA (Figure 2B)
logfc_matrix <- t_vol  %>%
  filter(facet_label %in% c("EE 45min vs Control", "EE 4hr vs Control", "EE 24hr vs Control", "RE 45min vs Control", "RE 4hr vs Control", "RE 24hr vs Control" )) %>%
  group_by(gene_symbol, facet_label) %>%
  dplyr::summarise(logFC = mean(logFC, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = facet_label, values_from = logFC) %>%
  drop_na() 
#PCA 
pca_result <- prcomp(t(logfc_matrix %>% dplyr::select(-gene_symbol)), scale = TRUE)
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100 # to label PC percentage

pca_df <- as.data.frame(pca_result$x)
pca_df$contrast <- colnames(logfc_matrix)[-1]
pca_df <- pca_df %>%
  mutate(contrast = case_when(
    contrast == "EE 45min vs Control" ~ "EE 45min",
    contrast == "EE 4hr vs Control" ~ "EE 4hr",
    contrast == "EE 24hr vs Control" ~ "EE 24hr",
    contrast == "RE 45min vs Control" ~ "RE 45min",
    contrast == "RE 4hr vs Control" ~ "RE 4hr",
    contrast == "RE 24hr vs Control" ~ "RE 24hr",
    TRUE ~ contrast  # Keep other values unchanged
  ))
# Figure 2B (top)
ggplot(pca_df, aes(x = PC1, y = PC2, label = contrast, color = contrast)) +
  geom_point(size = 4) +
  geom_text_repel(
    max.overlaps = 20,
    size = 3,
    box.padding = 0.5,
    segment.color = "gray60",
    force = 1
  ) +
  theme_minimal() +
  labs(
    title = "PCA of Transcriptomic responses",
    x = paste0("PC1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 2), "%)")
  ) +
  scale_color_manual(
    values = c("EE 45min" = "#d95f02", 
               "EE 4hr" = "#d95f02", 
               "EE 24hr" = "#d95f02",
               "RE 45min" = "#1b9e77", 
               "RE 4hr" = "#1b9e77", 
               "RE 24hr" = "#1b9e77"),
    name = "Contrast"
  ) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# 2. Proteomics
precawg_prot_da <- motrpac_da$adipose$`prot-pr` %>%
  as.data.frame() %>%
  left_join(
    HUMAN_FEATURE_TO_GENE %>% select(feature_id, uniprot, gene_symbol),
    by = "feature_id"
  )
pr_vol <- map_df(
  names(prot_list),
  function(name) {
    contrast_value <- prot_list[[name]]
    
    precawg_prot_da %>%
      filter(contrast == contrast_value) %>%
      mutate(
        facet_label = gsub("_", " ", name),
        is_ee_vs_re = grepl("EE_vs_RE", name),
        color = case_when(
          adj_p_value >= 0.05 ~ "grey",
          is_ee_vs_re & logFC > 0 ~ "#d95f02",
          is_ee_vs_re & logFC < 0 ~ "#1b9e77",
          grepl("EE", name) ~ "#d95f02",
          grepl("RE", name) ~ "#1b9e77",
          TRUE ~ "#7570b3"
        )
      ) %>%
      group_by(facet_label) %>%
      mutate(
        gene_symbol = as.character(gene_symbol),
        p_rank = rank(p_value),
        label = ifelse(adj_p_value < 0.05 & p_rank <= 10, gene_symbol, NA)
      ) %>%
      ungroup()
  }
)

# 3. Phosphoproteomics
precawg_phos_da <- as.data.frame(motrpac_da$adipose$`prot-ph`)
precawg_phos_da <- precawg_phos_da %>%
  # Join with HUMAN_FEATURE_TO_GENE to map feature_id to gene_symbol
  left_join(
    HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol),
    by = "feature_id"
  ) %>%
  # Process the feature_id after mapping
  mutate(
    # Extract base_feature_id (protein ID) before the first underscore
    base_feature_id = sub("_.*$", "", feature_id),  
    
    # Extract phosphosite(s) from everything after the first underscore
    phosphosites = sub("^.*_", "", feature_id) %>%
      gsub("s$", "", .) %>%  # Remove trailing 's'
      gsub("([STY]\\d+)[a-zA-Z]", "\\1;", .) %>%  # Replace any extra letters after phosphosites with ';'
      gsub(";$", "", .),  # Remove trailing semicolon
    
    # Create gene_symbol_with_phosphosite
    gene_symbol_with_phosphosite = ifelse(
      !is.na(gene_symbol),
      paste0(gene_symbol, "-", phosphosites),  # Combine gene symbol and cleaned phosphosites
      feature_id  # Fallback to feature_id if gene_symbol is NA
    )
  )
ph_vol <- map_df(
  names(prot_list),
  function(name) {
    contrast_value <- prot_list[[name]]
    
    precawg_phos_da %>%
      filter(contrast == contrast_value) %>%
      mutate(
        facet_label = gsub("_", " ", name),
        is_ee_vs_re = grepl("EE_vs_RE", name),
        color = case_when(
          adj_p_value >= 0.05 ~ "grey",
          is_ee_vs_re & logFC > 0 ~ "#d95f02",
          is_ee_vs_re & logFC < 0 ~ "#1b9e77",
          grepl("EE", name) ~ "#d95f02",
          grepl("RE", name) ~ "#1b9e77",
          TRUE ~ "#7570b3"
        ),
        
        # Top-10 labels using gene_symbol_with_phosphosite
        label = ifelse(
          adj_p_value < 0.05 & rank(p_value) <= 10,
          gene_symbol_with_phosphosite,
          NA
        )
      )
  }
)

# 4. Metabolomics
precawg_metab_da <- motrpac_da$adipose %>%
  keep(~ grepl("metab", .x$assay[1])) %>% # Select elements containing "metab" in assay column
  bind_rows()
# Combine all dataframes by rows, allowing for differing columns
qc_norm_list <- list()

# Loop through all objects starting with "metab"
for (metab_name in names(motrpac_qc$adipose)) {
  if (startsWith(metab_name, "metab")) {
    # Extract qc_norm and sample_metadata for the current metab object
    qc_norm <- motrpac_qc$adipose[[metab_name]]$qc_norm
    sample_metadata <- motrpac_qc$adipose[[metab_name]]$sample_metadata
    
    # Filter sample_metadata to keep only pre-exercise samples
    pre_exercise_samples <- sample_metadata$Timepoint == "pre_exercise" & sample_metadata$visitcode == "ADU_BAS"
    filtered_metadata <- sample_metadata[pre_exercise_samples, ]
    
    # Subset qc_norm to include only pre-exercise vialLabels
    pre_exercise_qc_norm <- qc_norm[, colnames(qc_norm) %in% filtered_metadata$vialLabel, drop = FALSE]
    
    # Rename vialLabels (columns) to their corresponding pid from sample_metadata
    colnames(pre_exercise_qc_norm) <- filtered_metadata$pid[match(colnames(pre_exercise_qc_norm), filtered_metadata$vialLabel)]
    
    # Add the filtered qc_norm to the list
    qc_norm_list[[metab_name]] <- pre_exercise_qc_norm
  }
}

m_vol <- map_df(
  names(subset_list),
  function(name) {
    contrast_value <- subset_list[[name]]
    
    precawg_metab_da %>%
      filter(contrast == contrast_value) %>%
      mutate(
        facet_label = gsub("_", " ", name),
        is_ee_vs_re = grepl("EE_vs_RE", name),
        color = case_when(
          adj_p_value >= 0.05 ~ "grey",
          is_ee_vs_re & logFC > 0 ~ "#d95f02",
          is_ee_vs_re & logFC < 0 ~ "#1b9e77",
          grepl("EE", name) ~ "#d95f02",
          grepl("RE", name) ~ "#1b9e77",
          TRUE ~ "#7570b3"
        ),
        p_rank = rank(p_value),
        
        # use feature_id (s label
        label = ifelse(adj_p_value < 0.05 & p_rank <= 6, 
                       as.character(feature_id), 
                       NA)
      )
  }
)
# Metab PCA
logfc_matrix_metab <- m_vol  %>%
  filter(facet_label %in% c("EE 45min vs Control", "EE 4hr vs Control", "EE 24hr vs Control", "RE 45min vs Control", "RE 4hr vs Control", "RE 24hr vs Control" )) %>%
  group_by(feature_id, facet_label) %>%
  dplyr::summarise(logFC = mean(logFC, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = facet_label, values_from = logFC) %>%
  drop_na() 
#PCA 
pca_result <- prcomp(t(logfc_matrix_metab %>% dplyr::select(-feature_id)), scale = TRUE)
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100 # to label PC percentage

pca_df <- as.data.frame(pca_result$x)
pca_df$contrast <- colnames(logfc_matrix_metab)[-1]
pca_df <- pca_df %>%
  mutate(contrast = case_when(
    contrast == "EE 45min vs Control" ~ "EE 45min",
    contrast == "EE 4hr vs Control" ~ "EE 4hr",
    contrast == "EE 24hr vs Control" ~ "EE 24hr",
    contrast == "RE 45min vs Control" ~ "RE 45min",
    contrast == "RE 4hr vs Control" ~ "RE 4hr",
    contrast == "RE 24hr vs Control" ~ "RE 24hr",
    TRUE ~ contrast  # Keep other values unchanged
  ))

# Figure 2B (bottom)
ggplot(pca_df, aes(x = PC1, y = PC2, label = contrast, color = contrast)) +
  geom_point(size = 4) +
  geom_text_repel(
    max.overlaps = 20,
    size = 3,
    box.padding = 0.5,
    segment.color = "gray60",
    force = 1
  ) +
  theme_minimal() +
  labs(
    title = "PCA of Metab responses",
    x = paste0("PC1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 2), "%)")
  ) +
  scale_color_manual(
    values = c("EE 45min" = "#d95f02", 
               "EE 4hr" = "#d95f02", 
               "EE 24hr" = "#d95f02",
               "RE 45min" = "#1b9e77", 
               "RE 4hr" = "#1b9e77", 
               "RE 24hr" = "#1b9e77"),
    name = "Contrast"
  ) +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Plot volcano plots for each omic layer
selected_facets <- c(
  "EE 45min vs Control", "EE 4hr vs Control", "EE 24hr vs Control", 
  "RE 45min vs Control", "RE 4hr vs Control", "RE 24hr vs Control"
)
all_vol <- bind_rows(
  t_vol %>% mutate(Data_Type = "Transcriptomics", label = as.character(label)),
  pr_vol %>% mutate(Data_Type = "Proteomics", label = as.character(label)),
  ph_vol %>% mutate(Data_Type = "Phosphoproteomics", label = as.character(label)),
  m_vol %>% mutate(Data_Type = "Metabolomics", label = as.character(label))
) %>%
  filter(facet_label %in% selected_facets) %>%
  mutate(
    Data_Type = factor(Data_Type, levels = c("Transcriptomics", "Proteomics", "Phosphoproteomics", "Metabolomics")),
    facet_label = factor(facet_label, levels = selected_facets)
  )

all_vol <- all_vol %>%
  mutate(
    Exercise = case_when(
      str_starts(facet_label, "EE") ~ "EE",
      str_starts(facet_label, "RE") ~ "RE",
      TRUE ~ NA_character_
    ),
    Clean_Time = facet_label %>% 
      str_remove(" vs Control") %>%  
      str_remove("^EE ") %>%          
      str_remove("^RE "),              
    Clean_Time = factor(Clean_Time, levels = c("45min", "4hr", "24hr"))
  )

# Custom function for volcano plots
plot_volcano <- function(omic_type, data) {
  
  df <- data %>% filter(Data_Type == omic_type)
  
  df <- df %>%
    mutate(
      color2 = case_when(
        adj_p_value >= 0.05 ~ "grey80",
        
        # EE vs RE: use direction
        Exercise == "EE vs RE" & logFC > 0 ~ "#d95f02",  # EE
        Exercise == "EE vs RE" & logFC < 0 ~ "#1b9e77",  # RE
        
        # EE vs CON (non-directional)
        Exercise == "EE" ~ "#d95f02",
        
        # RE vs CON (non-directional)
        Exercise == "RE" ~ "#1b9e77",
        
        # Control only
        Exercise == "CON" ~ "#7570b3",
        
        TRUE ~ "#7570b3"
      )
    )
  
  ggplot(df, aes(x = logFC, y = -log10(adj_p_value), color = color2)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_identity() +
    facet_nested(
      cols = vars(Clean_Time),
      rows = vars(Exercise),
      nest_line = element_line(color = "black", linewidth = 1),
      scales = "free_y",
      strip = strip_themed(
        background_x = elem_list_rect(fill = "white", color = "black"),
        background_y = elem_list_rect(fill = "white", color = "black")
      )
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_label_repel(
      aes(label = label),
      size = 2,
      max.overlaps = 22,
      box.padding = 0.1,
      point.padding = 0.1,
      label.padding = 0.1,
      segment.size = 0.15,
      segment.color = "grey50",
      label.size = 0.2,
      fill = "white",
      color = "black"
    ) +
    labs(title = omic_type, x = "Log Fold Change", y = "-log10 Adjusted P-Value") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_line(size = 0.1),
      legend.position = "none",
      strip.text.y = element_text(size = 16, face = "bold"),
      strip.text.x = element_text(size = 16, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# Figure 2A
p_trans <- plot_volcano("Transcriptomics", all_vol)
p_prot  <- plot_volcano("Proteomics", all_vol)
p_phos  <- plot_volcano("Phosphoproteomics", all_vol)
p_metab <- plot_volcano("Metabolomics", all_vol)

## Make volcano for Figure S2: CON-unadjusted DA
selected_facets2 <- c(
  "Control 45min", "Control 4hr", "Control 24hr",
  "EE 45min", "EE 4hr", "EE", 
  "RE 45min", "RE 4hr", "RE 24hr"
)

all_vol2 <- bind_rows(
  t_vol %>% mutate(Data_Type = "Transcriptomics"),
  pr_vol %>% mutate(Data_Type = "Proteomics"),
  ph_vol %>% mutate(Data_Type = "Phosphoproteomics"),
  m_vol %>% mutate(Data_Type = "Metabolomics")
) %>%
  filter(facet_label %in% selected_facets2) %>%
  mutate(
    Data_Type = factor(Data_Type, levels = c("Transcriptomics", "Proteomics", "Phosphoproteomics", "Metabolomics")),
    facet_label = factor(facet_label, levels = selected_facets2)  
  )
all_vol2 <- all_vol2 %>%
  mutate(
    Exercise = case_when(
      str_starts(facet_label, "EE") ~ "EE",
      str_starts(facet_label, "RE") ~ "RE",
      str_starts(facet_label, "Control") ~ "CON",
      TRUE ~ NA_character_
    ),
    Clean_Time = facet_label %>% 
      str_remove(" vs Control") %>%  
      str_remove("^EE ") %>%          
      str_remove("^RE ") %>%          
      str_remove("^Control "),         
    Clean_Time = factor(Clean_Time, levels = c("45min", "4hr", "24hr"))  
  ) %>%
  drop_na(Clean_Time) 

# Figure S2A
p_trans2 <- plot_volcano("Transcriptomics", all_vol2)
p_prot2  <- plot_volcano("Proteomics", all_vol2)
p_phos2  <- plot_volcano("Phosphoproteomics", all_vol2)
p_metab2 <- plot_volcano("Metabolomics", all_vol2)

###Make volcano for Figure3: Comparing EE vs. RE
selected_facets3 <- c(
  "EE 45min vs RE",  "EE 4hr vs RE" , "EE 24hr vs RE"
)

all_vol3 <- bind_rows(
  t_vol %>% mutate(Data_Type = "Transcriptomics"),
  pr_vol %>% mutate(Data_Type = "Proteomics"),
  ph_vol %>% mutate(Data_Type = "Phosphoproteomics"),
  m_vol %>% mutate(Data_Type = "Metabolomics")
) %>%
  filter(facet_label %in% selected_facets3) %>%
  mutate(
    Data_Type = factor(Data_Type, levels = c("Transcriptomics", "Proteomics", "Phosphoproteomics", "Metabolomics")),
    facet_label = factor(facet_label, levels = selected_facets3)  
  )
all_vol3 <- all_vol3 %>%
  mutate(
    Exercise = ifelse(str_detect(facet_label, "EE .* vs RE"), "EE vs RE", NA_character_),  # ✅ Detects "EE (any time) vs RE"
    Clean_Time = facet_label %>% 
      str_remove(" vs RE") %>%   
      str_remove("^EE "),        
    Clean_Time = factor(Clean_Time, levels = c("45min", "4hr", "24hr"))  
  ) %>%
  drop_na(Clean_Time)  

# Figure 3A
p_trans3 <- plot_volcano("Transcriptomics", all_vol3)
p_prot3  <- plot_volcano("Proteomics", all_vol3)
p_phos3  <- plot_volcano("Phosphoproteomics", all_vol3)
p_metab3 <- plot_volcano("Metabolomics", all_vol3)


# Integrated z.std comparison of ee vs re hits
eere_hits <- all_vol3 %>%
  filter(contrast_type == "Endur_vs_Resist", adj_p_value < 0.05) %>%
  select(feature_id, Clean_Time)
eere_exercise_filtered <- all_vol %>%
  inner_join(eere_hits, by = c("feature_id", "Clean_Time"))
eere_wide <- eere_exercise_filtered %>%
  mutate(
    label_name = case_when(
      assay == "transcript-rna-seq" ~ gene_symbol,
      assay == "prot-pr" ~ gene_symbol,
      assay == "prot-ph" ~ gene_symbol_with_phosphosite,
      assay == "metab" ~ feature_id,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(contrast_category %in% c("EE-CON", "RE-CON")) %>%
  select(feature_id, label_name, assay, Clean_Time, contrast_category, z.std) %>%
  pivot_wider(names_from = contrast_category, values_from = z.std) %>%
  drop_na(`EE-CON`, `RE-CON`)  # only keep features with both values

assay_colors <- c(
  "metab" = "#6D4B08",
  "prot-pr" = "#228833",
  "prot-ph" = "#F3A02B",
  "transcript-rna-seq" = "#377EB8"
)
label_df <- bind_rows(
  eere_wide %>%
    filter((`EE-CON` < 0 & `RE-CON` > 0) | (`EE-CON` > 0 & `RE-CON` < 0)) %>%
    filter(assay == "prot-ph") %>%
    slice_max(order_by = abs(`EE-CON` - `RE-CON`), n = 15, with_ties = FALSE),
  
  eere_wide %>%
    filter((`EE-CON` < 0 & `RE-CON` > 0) | (`EE-CON` > 0 & `RE-CON` < 0)) %>%
    filter(assay != "prot-ph") %>%
    group_by(assay) %>%
    slice_max(order_by = abs(`EE-CON` - `RE-CON`), n = 7, with_ties = FALSE) %>%
    ungroup(),
  eere_wide %>% filter(label_name %in% c("NR4A1", "NDRG1-S336", "SM 42:2;O2"))
  
)
# Figure 3B
ggplot(eere_wide, aes(x = `EE-CON`, y = `RE-CON`, color = assay, shape = Clean_Time)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(
    values = assay_colors,
    labels = c(
      "transcript-rna-seq" = "Transcriptomics",
      "prot-pr" = "Proteomics",
      "prot-ph" = "Phosphoproteomics",
      "metab" = "Metabolomics"
    )
  ) +
  geom_label_repel(
    data = label_df,
    aes(label = label_name),
    box.padding = 0.3,
    max.overlaps = 15,
    size = 4,
    label.size = 0.2,
    fill = "white",
    show.legend = FALSE
  ) +
  labs(
    x = "Z.std (EE-CON)",
    y = "Z.std (RE-CON)",
    color = "Assay",
    shape = "Timepoint",
    title = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

# DA results involving cell type specificity
# Load cell marker list
load("celltype_markers.RData")
ba_svf_list <- ba_svf_list[!duplicated(names(ba_svf_list))]

gene_cell_map <- stack(ba_svf_list) %>%
  rename(Gene = values, CellType = ind) %>%
  mutate(CellType = as.character(CellType))

all_vol_celltype <- all_vol %>%
  filter(Data_Type != "Metabolomics") %>%
  left_join(gene_cell_map, by = c("gene_symbol" = "Gene")) %>%
  mutate(CellType = ifelse(is.na(CellType), "Non-marker", as.character(CellType))) %>%
  filter(facet_label %in% selected_facets) %>%
  mutate(
    Data_Type = factor(Data_Type, levels = c("Transcriptomics", "Proteomics", "Phosphoproteomics")),
    facet_label = factor(facet_label, levels = selected_facets)  
  )

## Apply cell type on previous DA volcano dfs.   
all_vol_celltype_fil <- all_vol_celltype %>%
  filter(adj_p_value<0.05)  %>%
  mutate(
    label = case_when(
      adj_p_value < 0.05 & CellType != "Non-marker" & Data_Type %in% c("Transcriptomics", "Proteomics") ~ gene_symbol,
      adj_p_value < 0.05 & CellType != "Non-marker" & Data_Type == "Phosphoproteomics" ~ gene_symbol_with_phosphosite,
      TRUE ~ NA_character_
    )
  )


# Define colors 
cell_types <-c(
  "Adip_1", "Adip_2", "Vascular", "Pre_Ad", "Stem", "LAM", "Resident", "Mast", "NK.T"
)
#cell_types <- sort(cell_types)
# Assign colors using 'Paired' from RColorBrewer
paired_colors <- brewer.pal(n = length(cell_types), name = "Paired")
# Create color mapping, keeping "Non-marker" grey
color_mapping <- setNames(paired_colors, cell_types)
color_mapping["Non-marker"] <- "grey90"


plot_volcano_celltype <- function(omic_type, data) {
  data <- data %>%
    filter(Data_Type == omic_type) 
  ggplot(data, aes(x = logFC, y = -log10(adj_p_value), color = CellType, alpha = CellType != "Non-marker")) +  
    geom_point(size = 1.5) +  
    
    # Use custom colors for CellType
    scale_color_manual(values = color_mapping) +
    
    # Make Non-marker more transparent
    scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.2), guide = "none") +  
    
    # Nested facets
    facet_nested(
      cols = vars(Clean_Time),
      rows = vars(Exercise),
      nest_line = element_line(color = "black", linewidth = 1),
      scales = "free_y",
      strip = strip_themed(
        background_x = elem_list_rect(fill = "white", color = "black"),
        background_y = elem_list_rect(fill = "white", color = "black")
      )
    ) +
    
    # Vertical dashed line at x=0
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.8) +
    
    # Label only top 20 per facet dynamically!
    geom_label_repel(
      aes(label = label, fill = CellType),
      size = 2.5,
      max.overlaps = 25,  
      force = 3,  
      direction = "both",  
      nudge_y = 0.5,
      box.padding = 0.2,
      label.padding = 0.2,
      segment.size = 0.2,
      segment.color = "grey50",
      label.size = 0.2,
      color = "black"
    ) +  
    scale_fill_manual(values = color_mapping) +  
    
    # Titles and themes
    labs(title = omic_type, x = "Log Fold Change", y = "-log10 Adjusted P-Value") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(size = 0.2),
      panel.grid.minor = element_line(size = 0.1),
      legend.position = "right",
      strip.text.y = element_text(size = 16, face = "bold"),
      strip.text.x = element_text(size = 16, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}
# Figure 4B
p_trans_ct <- plot_volcano_celltype("Transcriptomics", all_vol_celltype_fil)
p_prot_ct  <- plot_volcano_celltype("Proteomics", all_vol_celltype_fil)
p_phos_ct  <- plot_volcano_celltype("Phosphoproteomics", all_vol_celltype_fil)

# Cell type labeling on EE vs. RE
all_vol_celltype_eere <- all_vol3 %>%
  filter(Data_Type != "Metabolomics") %>%
  left_join(gene_cell_map, by = c("gene_symbol" = "Gene")) %>%
  mutate(CellType = ifelse(is.na(CellType), "Non-marker", as.character(CellType))) %>%
  mutate(
    Data_Type = factor(Data_Type, levels = c("Transcriptomics", "Proteomics", "Phosphoproteomics")),
    facet_label = factor(facet_label, levels = selected_facets3)  # ✅ Keep full labels for filtering
  )

all_vol_celltype_eere_fil <- all_vol_celltype_eere %>%
  filter(adj_p_value<0.05)  %>%
  mutate(
    label = case_when(
      adj_p_value < 0.05 & CellType != "Non-marker" & Data_Type %in% c("Transcriptomics", "Proteomics") ~ gene_symbol,
      adj_p_value < 0.05 & CellType != "Non-marker" & Data_Type == "Phosphoproteomics" ~ gene_symbol_with_phosphosite,
      TRUE ~ NA_character_
    )
  )

# Figure S4B
p_prot_ct_eere  <- plot_volcano_celltype("Transcriptomics", all_vol_celltype_eere_fil)
p_phos_ct_eere  <- plot_volcano_celltype("Phosphoproteomics", all_vol_celltype_eere_fil)



### Feature heatmap - ceramide (Metab) 
cer_features <- precawg_metab_da %>%
  filter(str_detect(feature_id, "Cer")) %>%
  pull(feature_id) %>%
  unique()
cer_features <- cer_features[!str_detect(cer_features, "HexCer")]

create_feature_heatmap <- function(data, contrasts, order, features = NULL, gene_col) {
  # Step 1: Filter data for the specified contrasts
  heatmap_data <- data %>%
    filter(contrast %in% contrasts) %>%
    mutate(
      group = ifelse(grepl("^group_timepointADUEndur", contrast), "EE", "RE"),
      timeline = case_when(
        grepl("post_15_30_45_min", contrast) ~ "Post 45 min",
        grepl("post_3.5_4_hr", contrast) ~ "Post 4 hr",
        grepl("post_24_hr", contrast) ~ "Post 24 hr"
      ),
      contrast = factor(contrast, levels = order) # Ensure correct contrast order
    ) %>%
    arrange(contrast)
  
  # If features are provided, filter data to include only those features
  if (!is.null(features)) {
    heatmap_data <- heatmap_data %>%
      filter(!!sym(gene_col) %in% features)
  }
  
  # Step 2: Pivot data to wide format
  logFC_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, logFC) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = logFC,
      values_fill = 0
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  # Create a significance mask for asterisks
  significance_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, adj_p_value) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = adj_p_value,
      values_fill = Inf
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  significance_matrix <- ifelse(significance_matrix < 0.05, "*", "")
  
  # Step 3: Column annotations
  col_annotations <- data.frame(
    contrast = order,
    group = ifelse(grepl("^group_timepointADUEndur", order), "EE", "RE"),
    timeline = rep(c("Post 45 min", "Post 4 hr", "Post 24 hr"), times = 2) # Repeat for EE and RE
  )
  
  # Set order for timeline legend
  col_annotations$timeline <- factor(col_annotations$timeline, levels = c("Post 45 min", "Post 4 hr", "Post 24 hr"))
  
  col_annotation <- HeatmapAnnotation(
    Exercise = col_annotations$group,
    Timeline = col_annotations$timeline,
    col = list(
      Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77"),
      Timeline = c("Post 45 min" = "#AE76A3", "Post 4 hr" = "#882E72", "Post 24 hr" = "#61194F")
    ),
    annotation_legend_param = list(
      Exercise = list(title = "Modality"),
      Timeline = list(title = "Timepoint", at = c("Post 45 min", "Post 4 hr", "Post 24 hr"))
    ),
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm"),
    height=unit(0.5, "mm")
  )
  
  # Step 4: Heatmap creation
  heatmap <- Heatmap(
    logFC_matrix,
    name = "logFC",
    col = colorRamp2(c(-max(abs(logFC_matrix)), 0, max(abs(logFC_matrix))), c("blue", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_dend_side = "right", # Move row dendrogram to the right
    row_names_side = "left", # Move row labels (gene names) to the left
    column_split = col_annotations$group,
    column_title = NULL,
    top_annotation = col_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (significance_matrix[i, j] == "*") {
        # Convert fill (RGB) to luminance
        rgb_vals <- col2rgb(fill)
        luminance <- (0.299 * rgb_vals[1, ] + 0.587 * rgb_vals[2, ] + 0.114 * rgb_vals[3, ]) / 255
        
        # Decide text color based on brightness
        text_color <- ifelse(luminance < 0.5, "white", "black")
        
        grid.text(
          significance_matrix[i, j],
          x,
          y,
          gp = gpar(fontsize = 20, col = text_color)
        )
      }
    }
  )
  
  return(heatmap)
}
contrasts_to_include <- c(
  "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
)
contrast_order <- c(
  "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
)
# Figure 3G
create_feature_heatmap(
  data = precawg_metab_da,
  contrasts = contrasts_to_include,
  order = contrast_order,
  features = cer_features,
  gene_col = "feature_id" # or "feature_id"
)

# Custom heatmap function for proteomics or phosphoproteomics
create_4hr_heatmap <- function(data, contrasts, features = NULL, gene_col) {
  
  # Step 1: Filter data for the specified contrasts
  heatmap_data <- data %>%
    filter(contrast %in% contrasts) %>%
    mutate(
      group = ifelse(grepl("^group_timepointADUEndur", contrast), "EE", "RE"),
      timeline = "4hrPost",
      contrast = factor(contrast, levels = contrasts)
    ) %>%
    arrange(contrast)
  
  # If features are provided, filter to those features
  if (!is.null(features)) {
    heatmap_data <- heatmap_data %>%
      filter(!!sym(gene_col) %in% features)
  }
  
  # STEP to solve duplicate rows: keep the one with smallest adj_p_value
  heatmap_data <- heatmap_data %>%
    group_by(!!sym(gene_col), contrast) %>%
    slice_min(adj_p_value, with_ties = FALSE) %>%
    ungroup()
  
  # Step 2: Pivot to wide format (logFC matrix)
  logFC_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, logFC) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = logFC,
      values_fill = 0
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  # Step 3: Create significance matrix for asterisks
  significance_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, adj_p_value) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = adj_p_value,
      values_fill = Inf
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  significance_matrix <- apply(significance_matrix, c(1, 2), function(p) {
    if (p < 0.05) {
      "*"
    } else if (p < 0.1) {
      "\u2020"  # Unicode for dagger (†)
    } else {
      ""
    }
  })
  
  # Step 4: Column annotations
  col_annotations <- data.frame(
    contrast = contrasts,
    group = ifelse(grepl("^group_timepointADUEndur", contrasts), "EE", "RE"),
    timeline = "4hrPost"
  )
  
  col_annotations$timeline <- factor(col_annotations$timeline, levels = c("4hrPost"))
  
  col_annotation <- HeatmapAnnotation(
    Exercise = col_annotations$group,
    Timeline = col_annotations$timeline,
    col = list(
      Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77"),
      Timeline = c("4hrPost" = "#882E72")
    ),
    annotation_legend_param = list(
      Exercise = list(title = "Modality"),
      Timeline = list(title = "Timepoint")
    ),
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm"),
    height=unit(0.5, "mm")
  )
  
  # Step 5: Heatmap
  heatmap <- Heatmap(
    logFC_matrix,
    name = "logFC",
    col = colorRamp2(
      c(-max(abs(logFC_matrix)), 0, max(abs(logFC_matrix))),
      c("blue", "white", "red")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_dend_side = "right",
    row_names_side = "left",
    column_split = col_annotations$group,
    column_title = NULL,
    top_annotation = col_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      symbol <- significance_matrix[i, j]
      if (symbol != "") {
        font_color <- if (symbol == "*") "white" else "black"
        font_size  <- if (symbol == "*") 20 else 18
        grid.text(
          label = symbol,
          x = x,
          y = y,
          gp = gpar(fontsize = font_size, col = font_color)
        )
      }
    }
  )
  
  return(heatmap)
}


contrasts_to_include2 <- c(
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise"
)
# Figure S2H. A dagger sign was manually added to TNS1-S1343;S1369 (EE) given its adj. p=0.1004
tns1 <- c("TNS1-T1622", "TNS1-S1307;S1314", "TNS1-S752", "TNS1-S1343;S1369")
create_4hr_heatmap(
  data = precawg_phos_da,
  contrasts = contrasts_to_include2,
  features = tns1,
  gene_col = "gene_symbol_with_phosphosite" # or "feature_id"
)

# Correlation between ROCK1 and TNS1
# Transpose and prepare expression data
prot_qc <- motrpac_qc$adipose$`prot-pr`$qc_norm
protein_subset <- prot_qc[c("Q9HBL0", "Q13464"), ]
protein_df <- as.data.frame(t(protein_subset))
colnames(protein_df) <- c("TNS1", "ROCK1")
protein_df$vialLabel <- rownames(protein_df)

# Ensure numeric and remove NA
protein_df <- protein_df %>%
  mutate(across(c(TNS1, ROCK1), as.numeric)) %>%
  drop_na()

# Merge with metadata
meta_df <- motrpac_qc$adipose$`prot-pr`$sample_metadata %>%
  mutate(vialLabel = as.character(vialLabel)) %>%
  select(vialLabel, randomGroupCode, Timepoint)

# Ensure vialLabel is character
protein_df$vialLabel <- as.character(protein_df$vialLabel)
meta_df$vialLabel <- as.character(meta_df$vialLabel)

# Join and filter Timepoint to valid levels
protein_df <- protein_df %>%
  left_join(meta_df, by = "vialLabel") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("pre_exercise", "post_3.5_4_hr")),
    shape_group = factor(
      if_else(Timepoint == "pre_exercise", "Pre", "Post 4hr"),
      levels = c("Pre", "Post 4hr")
    ),
    group_label = recode(
      randomGroupCode,
      "ADUControl" = "CON",
      "ADUEndur"  = "EE",
      "ADUResist" = "RE"
    )
  )
# Compute Spearman correlation
cor_test <- cor.test(protein_df$TNS1, protein_df$ROCK1, method = "spearman")
cor_val <- cor_test$estimate
p_val <- cor_test$p.value
annot_text <- paste0("r = ", round(cor_val, 2),
                     ", p ", format.pval(p_val, digits = 3, eps = .001))

# Define named color vector with new group labels
group_colors_named <- c("CON" = "#7570b3", "EE" = "#d95f02", "RE" = "#1b9e77")

# Figure S2I
ggplot(protein_df, aes(x = TNS1, y = ROCK1, color = group_label, shape = shape_group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_smooth(
    aes(x = TNS1, y = ROCK1), 
    method = "lm", 
    se = TRUE, 
    color = "black", 
    linetype = "dashed", 
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = group_colors_named, name = "Group") +
  scale_shape_manual(values = c("Pre" = 1, "Post 4hr" = 16), name = "Timepoint") +
  annotate("text", 
           label = annot_text, 
           x = Inf, y = -Inf, 
           hjust = 1.05, vjust = -0.5, 
           size = 4.5, 
           fontface = "italic") +
  labs(
    title = NULL,
    x = "TNS1 protein",
    y = "ROCK1 protein"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

## Figure S2C: circadian genes from external datasets
#1. Fabre 2018
fabre_un <-  read_excel("Fabre_2018_untrained.xlsx")
cir_genes <- c("DBP", "PER1", "PER2", "NR1D1", "NR1D2", "CLOCK")

# Filter for those genes
fabre_plot_data <- fabre_un %>%
  filter(external_gene_id %in% cir_genes) %>%
  mutate(
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )

# Barplot
p1 <- ggplot(fabre_plot_data, aes(x = external_gene_id, y = logFC, fill = logFC > 0)) +
  geom_col(width = 0.6, color = "black", show.legend = FALSE) +
  geom_text(aes(label = significance), vjust = ifelse(fabre_plot_data$logFC > 0, -0.5, 1.5), size = 5) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(
    title = "4hr post aerobic exercise (80% VO2max) \n(Fabre et al., 2018)",
    x = NULL,
    y = "Log Fold Change (logFC)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # plot border
    axis.text.x = element_text(size = 12)  # ✅ axis text sizing here
  )

# 2. Ahn 2025
ahn_3x <-  read.csv("Ahn_2025_acute_ex.csv", row.names = 1)
ahn_3x_mod <- ahn_3x %>%
  filter(contrast == "MOD Post - Pre")

ahn_plot_data <- ahn_3x_mod %>%
  filter(GeneID %in% cir_genes) %>%
  mutate(
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )

p2 <- ggplot(ahn_plot_data, aes(x = GeneID, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col(width = 0.6, color = "black", show.legend = FALSE) +  # ✅ black bar outline
  geom_text(aes(label = significance), vjust = ifelse(ahn_plot_data$log2FoldChange > 0, -0.5, 1.5), size = 5) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(
    title = "1.5hr post aerobic exercise (65% VO2max) \n(Ahn et al., 2025)",
    x = NULL,
    y = "Log Fold Change (logFC)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # plot border
    axis.text.x = element_text(size = 12)  # ✅ axis text sizing here
  )

# Figure S2C
p1 + p2
