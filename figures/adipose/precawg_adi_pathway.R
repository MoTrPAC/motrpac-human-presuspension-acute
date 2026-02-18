# load packages
library(MotrpacHumanPreSuspension)
library(ComplexHeatmap)
library(patchwork)
library(ggrepel)
library(UpSetR)
library(dplyr)
library(TMSig)
library(tibble)
library(stringr)
library(tidyr)
library(purrr)

data("CAMERA_RESULTS") 
################### Enrichment related to Figure 2 ###################
### Transcriptomics - Figure 2D 
selected_contrasts <- c(
  "Endur.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)",
  "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
  "Endur.post_24_hr - Control.post_24_hr (delta-delta)",
  "Resist.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)",
  "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
  "Resist.post_24_hr - Control.post_24_hr (delta-delta)"
)

#camera_res_t$set_id <- as.character(camera_res_t$set_id)
selected_set_ids <- c(
  "12575", "11461", "13274", "04103", "11452", "14704", "13236", "11380", "07184", "11662",
  "14120", "14642", "14756", "11443", "13006", "12358", "00146", "11613", "07016", "11419",
  "11394", "14391", "00012", "11663", "00382", "12459", "00180", "11583", "14272", "11528",
  "14019", "00686", "07502", "06369", "13563", "12154", "03896", "00287", "12791", "03031",
  "10607", "00096", "12787"
)

# Filter the dataset for the selected set_id values
camera_res_t <- CAMERA_RESULTS %>%
  filter(tissue == 'adipose') %>%
  filter(assay == "transcript-rna-seq") %>%
  filter(contrast_short %in% selected_contrasts) %>%
  filter(set_id %in% selected_set_ids) %>%
  mutate(
    contrast_short = factor(contrast_short, levels = selected_contrasts)
  )


# Define exercise & control colors
Exercise_colors <- c("EE" = "#d95f02", "RE" = "#1b9e77")
time_colors <- c("45minPost" = "#AE76A3", "4hrPost" = "#882E72", "24hrPost" = "#61194F")  # Timepoint colors

# Extract Exercise Type & Time for annotation
camera_res_t <- camera_res_t %>%
  mutate(
    Timepoint = case_when(
      grepl("24_hr", contrast) ~ "24hrPost",
      grepl("4_hr", contrast) ~ "4hrPost",
      grepl("45_min", contrast) ~ "45minPost",
      TRUE ~ NA_character_  # fallback
    ),
    Exercise = case_when(
      grepl("Endur", contrast_short) ~ "EE",
      grepl("Resist", contrast_short) ~ "RE"
    )
  )
col_anno <- columnAnnotation(
  Exercise = rep(c("EE", "RE"), each = 3),
  Time = rep(c("45minPost", "4hrPost", "24hrPost"), times = 2),
  col = list(
    Exercise = Exercise_colors,
    Time = time_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Modality", at = c("EE", "RE")),
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))
  )
)
# Figure S2D 
camera_res_t_plot <- camera_res_t %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = FALSE,
    set_column = "set_short",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      column_split = rep(c("EE", "RE"), each = 3),
      show_column_names = FALSE
    )
  ) 

# Proteomics
selected_set_ids_pr <- c("08383","01614","05448","02037","02232","00457","08753","09356","04103","08528","10260","05317", "10546")
selected_contrasts_pr <- c(
  "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
  "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)"
)  
camera_res_pr <- CAMERA_RESULTS %>%
  filter(tissue == 'adipose') %>%
  filter(assay == "prot-pr") %>%
  filter(contrast_short %in% selected_contrasts_pr) %>%
  filter(set_id %in% selected_set_ids_pr) %>%
  mutate(
    contrast_short = factor(contrast_short, levels = selected_contrasts)
  )
camera_res_pr <- camera_res_pr %>%
  mutate(
    Timepoint = case_when(
      grepl("4_hr", contrast) ~ "4hrPost",
    ),
    Exercise = case_when(
      grepl("Endur", contrast_short) ~ "EE",
      grepl("Resist", contrast_short) ~ "RE"
    )
  )
col_anno <- columnAnnotation(
  Exercise = rep(c("EE", "RE"), each = 1),
  Time = rep("4hrPost", times = 2),
  col = list(
    Exercise = Exercise_colors,
    Time = time_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Modality", at = c("EE", "RE")),
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))
  )
)
# Figure S2E 
camera_res_pr_plot <- camera_res_pr %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = TRUE,
    set_column = "set_short",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      column_split = rep(c("EE", "RE"), each = 1),
      show_column_names = FALSE
    )
  ) 
#### Phospho ORA #####
motrpac_da <- load_differential_analysis()
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
phos_list <- precawg_phos_da %>%
  filter(
    contrast %in% c(
      "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
      "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise"
    ),
    adj_p_value < 0.1
  ) %>%
  group_by(contrast) %>%
  summarise(gene_symbols = list(unique(gene_symbol))) %>%
  deframe()  # Convert tibble to named list
phos_list <- precawg_phos_da %>%
  filter(
    contrast %in% c(
      "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
      "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise"
    ),
    adj_p_value < 0.1
  ) %>%
  mutate(direction = ifelse(logFC >= 0, "UP", "DOWN")) %>%
  group_by(contrast, direction) %>%
  summarise(gene_symbols = list(unique(gene_symbol)), .groups = "drop") %>%
  mutate(contrast_name = paste(contrast, direction, sep = "_")) %>%
  select(contrast_name, gene_symbols) %>%
  deframe()
bckg <- as.character(unique(precawg_phos_da$gene_symbol))

#run ora
ora_phos <- lapply(phos_list, function(input_i) {
  run_ORA(input = as.character(input_i),
          background = bckg,
          overlap_cutoff = 0)
}) %>% 
  bind_rows(.id = "contrast") 

top_terms <- ora_phos %>%
  group_by(contrast) %>%
  slice_min(adj_p_value, n = 6, with_ties = FALSE) %>%  # Select top 5 lowest adj_p_value per contrast
  pull(set_short) %>%  
  unique()
ora_phos <- ora_phos %>%
  mutate(
    contrast = case_when(
      contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise_UP" ~ "EE_4hrPost_up",
      contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise_DOWN" ~ "EE_4hrPost_down",
      contrast == "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise_UP" ~ "RE_4hrPost_up",
      contrast == "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise_DOWN" ~ "RE_4hrPost_down",
      TRUE ~ contrast  # Keep other contrast names unchanged
    )
  )
ora_phos <- ora_phos %>%
  mutate(
    modality = sub("_.*", "", contrast),                      # EE / RE
    time = sub("^[^_]+_([^_]+)_.*$", "\\1", contrast),        # 4hrPost
    direction = sub(".*_(up|down)$", "\\1", contrast),        # up / down
    neglog10_p = -log10(p_value)                          # just -log10 p-value
  )

direction_colors <- c("up" = "red", "down" = "blue")
contrasts_used <- ora_phos %>%
  filter(set_short %in% top_terms) %>%
  pull(contrast) %>%
  unique()

# Match annotation to those contrasts only
modality_vec <- ora_phos$modality[match(contrasts_used, ora_phos$contrast)]
time_vec     <- ora_phos$time[match(contrasts_used, ora_phos$contrast)]
direction_vec<- ora_phos$direction[match(contrasts_used, ora_phos$contrast)]

col_anno <- ComplexHeatmap::columnAnnotation(
  Modality = modality_vec,
  Time = time_vec,
  Direction = direction_vec,
  col = list(
    Modality = Exercise_colors,
    Time = time_colors,
    Direction = direction_colors
  ),
  annotation_legend_param = list(
    Modality = list(title = "Modality"),
    Time = list(title = "Time"),
    Direction = list(title = "Direction")
  )
)
# Figure S2G
ora_phos_ht <- ora_phos %>% 
  filter(set_short %in% top_terms) %>% 
  enrichmap(
    n_top = Inf,
    set_column = "set_short",
    statistic_column = "neglog10_p",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    padj_cutoff = 0.05,
    colors = c("white", "#543483"),
    heatmap_args = list(
      heatmap_legend_param = list(
        title = latex2exp::TeX("$\\bf{Signed(-\\log_{10}(\\text{Adj P}))}$")
      ),
      top_annotation = col_anno,
      column_split = rep(c("EE", "RE"), each = 2),
      column_title = NULL,
      show_column_names = FALSE
    )) 

# Metabolomics
camera_res_m <- CAMERA_RESULTS %>%
  filter(tissue == "adipose") %>%
  filter(assay == "metab") %>%
  filter(contrast_short %in% selected_contrasts) %>%
  mutate(
    contrast_short = factor(contrast_short, levels = selected_contrasts)
  )

# Extract Exercise Type & Time for annotation
camera_res_m <- camera_res_m %>%
  mutate(
    Timepoint = case_when(
      grepl("45_min", contrast) ~ "45minPost",
      grepl("4_hr", contrast) ~ "4hrPost",
      grepl("24_hr", contrast) ~ "24hrPost"
    ),
    Exercise = case_when(
      grepl("Endur", contrast_short) ~ "EE",
      grepl("Resist", contrast_short) ~ "RE"
    )
  )
col_anno <- columnAnnotation(
  Exercise = rep(c("EE", "RE"), each = 3),
  Time = rep(c("45minPost", "4hrPost", "24hrPost"), times = 2),
  col = list(
    Exercise = Exercise_colors,
    Time = time_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Modality", at = c("EE", "RE")),
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))
  )
)
# Figure S2F
camera_res_m_plot <- camera_res_m %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = TRUE,
    set_column = "set",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      column_split = rep(c("EE", "RE"), each = 3),
      show_column_names = FALSE
    )
  ) 


#### compare enriched terms between with CON vs. withot CON
camera_adipose <- CAMERA_RESULTS %>%
  filter(tissue == "adipose") %>%
  mutate(
    timepoint = case_when(
      grepl("post_15_30_45_min", contrast_short) ~ "45min",
      grepl("post_3.5_4_hr", contrast_short) ~ "4hr",
      grepl("post_24_hr", contrast_short) ~ "24hr"
    ),
    modality = case_when(
      grepl("Endur", contrast_short) ~ "EE",
      grepl("Resist", contrast_short) ~ "RE"
    )
  )


camera_raw <- camera_adipose %>%
  filter(contrast_type == "exercise_no_controls", adj_p_value < 0.05) %>%
  pull(set_short) %>%
  unique()

camera_adj <- camera_adipose %>%
  filter(contrast_type == "exercise_with_controls", adj_p_value < 0.05) %>%
  pull(set_short) %>%
  unique()

sig_terms <- union(camera_raw, camera_adj)

camera_adipose2 <- camera_adipose %>%
  filter(set_short %in% sig_terms) %>%
  filter(collection %in% c("C2", "C5", "REFMET"))

camera_with <- camera_adipose2 %>%
  filter(!database %in% c("WP", "KEGG_MEDICUS", "PID")) %>%
  filter(contrast_type == "exercise_with_controls") %>%
  select(set_short, timepoint, modality, assay, z_with = z.std, p_with = adj_p_value)

camera_no <- camera_adipose2 %>%
  filter(!database %in% c("WP", "KEGG_MEDICUS", "PID")) %>%
  filter(contrast_type == "exercise_no_controls") %>%
  select(set_short, timepoint, modality, assay, z_no = z.std, p_no = adj_p_value)

camera_comp <- full_join(camera_with, camera_no,
                         by = c("set_short", "timepoint", "modality", "assay")) %>%
  filter(!is.na(z_with) | !is.na(z_no)) %>%
  mutate(
    significance = case_when(
      p_with < 0.05 & p_no < 0.05 ~ "Both",
      p_with < 0.05 ~ "CON-adjusted",
      p_no < 0.05 ~ "CON-unadjusted",
      TRUE ~ "None"
    )
  ) %>%
  filter(significance != "None")  # Keep only significant comparisons
camera_comp <- camera_comp %>%
  mutate(
    set_short = gsub("^(GOBP_|GOCC_|GOMF_|REACTOME_|BIOCARTA_)", "", set_short)
  )


plot_camera_scatter <- function(
    df,
    assay_filter,
    manual_labels = NULL,
    color_map = c("Both" = "#640D6B", "CON-adjusted" = "#FF3EA5", "CON-unadjusted" = "grey")
) {
  df <- df %>%
    filter(assay == assay_filter) %>%
    mutate(
      set_short = as.character(set_short),
      timepoint = factor(timepoint, levels = c("45min", "4hr", "24hr")),
      opposing = sign(z_no) != sign(z_with)
    )
  
  # Label manually specified + top 2 opposing terms per facet
  opposing_df <- df %>%
    filter(opposing) %>%
    group_by(timepoint, modality) %>%
    slice_max(order_by = abs(z_no - z_with), n = 0, with_ties = FALSE) %>%
    ungroup()
  
  label_df <- bind_rows(
    df %>% filter(set_short %in% manual_labels),
    opposing_df
  ) %>%
    distinct()
  label_df <- label_df %>%
    mutate(set_short = case_when(
      set_short == "BRANCHED_CHAIN_AMINO_ACID_CATABOLISM" ~ "BCAA_CATABOLISM",
      set_short == "OXIDATIVE_PHOSPHORYLATION" ~ "OXPHOS",
      set_short == "EUKARYOTIC_TRANSLATION_ELONGATION" ~ "TRANSLATION_ELONGATION",
      set_short == "REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS" ~ "SLITS_AND_ROBOS",
      set_short == "TOLL_LIKE_RECEPTOR_BINDING" ~ "TLR_BINDING",
      set_short == "REGULATION_OF_MRNA_PROCESSING" ~ "MRNA_PROCESSING",
      TRUE ~ set_short
    )) %>%
    mutate(set_short = str_replace_all(set_short, c(
      "INTERFERON" = "IFN",
      "ENDOPLASMIC_RETICULUM" = "ER",
      "MITOCHONDRIAL" = "MITO"
    )))
  ggplot(df, aes(x = z_no, y = z_with, color = significance)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = color_map) +
    guides(color = "none") +
    facet_grid(modality ~ timepoint) +
    geom_label_repel(
      data = label_df,
      aes(label = set_short),
      size = 5,
      box.padding = 0.4,
      nudge_x = 1.5,  # try 0.5 to 1
      nudge_y = 1.5,
      segment.color = "gray60",
      segment.curvature = -0.2,
      segment.angle = 20,
      max.overlaps = 8,
      show.legend = FALSE,
      label.size = 0.2,
      force = 3,
      fill = "white"
    ) +
    labs(
      title = paste(assay_filter),
      x = "Z.std (CON-unadjusted)",
      y = "Z.std (CON-adjusted)",
      color = "Contrast"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 20),
      axis.text = element_text(size = 20),
      axis.title.x = element_text(size = 20),  # ← add this
      axis.title.y = element_text(size = 20),  # ← add this
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
}
p1 <- plot_camera_scatter(
  camera_comp,
  "transcript-rna-seq",
  manual_labels = c(
    "INSULIN_RECEPTOR_RECYCLING", "IL6_PATHWAY", "VEGF_PATHWAY",
    "BRANCHED_CHAIN_AMINO_ACID_CATABOLISM", "TRANSLATION", "OXIDATIVE_PHOSPHORYLATION",
    "RHOBTB_GTPASE_CYCLE", "PROTEIN_ACETYLATION",
    "NEUTROPHIL_MIGRATION", "INTEGRIN_CELL_SURFACE_INTERACTIONS", "INTERFERON_ALPHA_BETA_SIGNALING", "REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS", "RIBOSOME",
    "TOLL_LIKE_RECEPTOR_BINDING", "REGULATION_OF_MRNA_PROCESSING"
  )
)

p2 <- plot_camera_scatter(
  camera_comp, "prot-pr",
  manual_labels = c(
    "GMP_BIOSYNTHETIC_PROCESS", "IMP_METABOLIC_PROCESS", "ENDOPLASMIC_RETICULUM_ORGANIZATION",
    "VESICLE_FUSION_WITH_GOLGI_APPARATUS", "OXIDATIVE_PHOSPHORYLATION", "CRISTAE_FORMATION",
    "STRUCTURAL_CONSTITUENT_OF_RIBOSOME", "COMPLEX_I_BIOGENESIS"
  )
)

p3 <- plot_camera_scatter(
  camera_comp, "metab",
  manual_labels = c(
    "Dicarboxylic acids", "Unsaturated FA", "NAE", "O-PE", "Saturated FA", "Acyl CoAs", "PE",
    "SM", "DG", "Hydroxy FA", "TCA acids", "PC", "LPC", "O-PC", "Cer", "TG",
    "Phenolic acids", "Pyrimidine rNDP", "Amino acids", "Acyl carnitines"
  )
)
# Figure 2C
final_plot <- p1 + p2 + p3 + plot_layout(ncol = 3, widths = c(3, 1, 3))


## Highlight couple examples of pathways that are affected by CON
# 1. GOBP: OXPHOS
trajectory_three_lines <- camera_adipose %>%
  filter(
    assay == "transcript-rna-seq",
    set_short == "GOBP_OXIDATIVE_PHOSPHORYLATION",
    contrast_short %in% c(
      "Resist.post_15_30_45_min - Resist.pre_exercise",
      "Resist.post_3.5_4_hr - Resist.pre_exercise",
      "Resist.post_24_hr - Resist.pre_exercise",
      "Resist.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)",
      "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
      "Resist.post_24_hr - Control.post_24_hr (delta-delta)",
      "Control.post_15_30_45_min - Control.pre_exercise",
      "Control.post_3.5_4_hr - Control.pre_exercise",
      "Control.post_24_hr - Control.pre_exercise"
    )
  ) %>%
  mutate(
    # Extract timepoint substring from contrast_short
    timepoint_raw = str_extract(contrast_short, "post_15_30_45_min|post_3.5_4_hr|post_24_hr"),
    
    # Convert to ordered factor for plotting
    timepoint = factor(timepoint_raw,
                       levels = c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"),
                       labels = c("45min", "4hr", "24hr")),
    
    # Significance for point fill
    sig_dot = if_else(adj_p_value < 0.05, "Significant", "Not Significant"),
    
    # Label each contrast type
    contrast_label = case_when(
      str_detect(contrast_short, "Resist.post.* - Resist.pre_exercise") ~ "RE CON-unadjusted",
      str_detect(contrast_short, "Resist.post.* - Control.post.*") ~ "RE CON-adjusted",
      str_detect(contrast_short, "Control.post.* - Control.pre_exercise") ~ "Control"
    ),
    
    # For coloring by modality
    modality_group = case_when(
      str_starts(contrast_label, "RE") ~ "RE",
      str_starts(contrast_label, "Control") ~ "Control",
      TRUE ~ NA_character_
    )
  )

# GOBP OXPHOS line plot for Figure 2C
ggplot(trajectory_three_lines, aes(
  x = timepoint,
  y = z.std,
  color = contrast_label,
  group = contrast_label,
  linetype = contrast_label
)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(fill = sig_dot), shape = 21, size = 3, color = "black", stroke = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "GOBP_OXPHOS",
    x = NULL,
    y = "CAMERA z.std",
    color = "Contrast Type",
    linetype = "Contrast Type",
    fill = "Significance"
  ) +
  scale_color_manual(values = c(
    "RE CON-unadjusted" = "#1b9e77",
    "RE CON-adjusted" = "#1b9e77",
    "Control" = "#7570b3"
  )) +
  scale_linetype_manual(values = c(
    "RE CON-unadjusted" = "dashed",
    "RE CON-adjusted" = "solid",
    "Control" = "dotdash"
  )) +
  scale_fill_manual(values = c(
    "Significant" = "black",
    "Not Significant" = "white"
  )) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# GOBP Neutrophil Migration
trajectory_three_lines2 <- camera_adipose %>%
  filter(
    assay == "transcript-rna-seq",
    set_short == "GOBP_NEUTROPHIL_MIGRATION",
    contrast_short %in% c(
      "Endur.post_15_30_45_min - Endur.pre_exercise",
      "Endur.post_3.5_4_hr - Endur.pre_exercise",
      "Endur.post_24_hr - Endur.pre_exercise",
      "Endur.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)",
      "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
      "Endur.post_24_hr - Control.post_24_hr (delta-delta)",
      "Control.post_15_30_45_min - Control.pre_exercise",
      "Control.post_3.5_4_hr - Control.pre_exercise",
      "Control.post_24_hr - Control.pre_exercise"
    )
  ) %>%
  mutate(
    # Extract timepoint substring from contrast_short
    timepoint_raw = str_extract(contrast_short, "post_15_30_45_min|post_3.5_4_hr|post_24_hr"),
    
    # Convert to ordered factor for plotting
    timepoint = factor(timepoint_raw,
                       levels = c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"),
                       labels = c("45min", "4hr", "24hr")),
    
    # Significance for point fill
    sig_dot = if_else(adj_p_value < 0.05, "Significant", "Not Significant"),
    
    # Label each contrast type
    contrast_label = case_when(
      # EE comparisons
      str_detect(contrast_short, "^Endur\\.post.* - Endur\\.pre_exercise") ~ "EE CON-unadjusted",
      str_detect(contrast_short, "^Endur\\.post.* - Control\\.post.*") ~ "EE CON-adjusted",
      
      # Control comparisons
      str_detect(contrast_short, "^Control\\.post.* - Control\\.pre_exercise") ~ "Control",
      
      TRUE ~ NA_character_
    ),
    
    # For coloring by modality
    modality_group = case_when(
      str_starts(contrast_label, "EE") ~ "EE",
      str_starts(contrast_label, "Control") ~ "Control",
      TRUE ~ NA_character_
    )
  )

#GOBP Neutrophil Migration line plot for Figure 2C
ggplot(trajectory_three_lines2, aes(
  x = timepoint,
  y = z.std,
  color = contrast_label,
  group = contrast_label,
  linetype = contrast_label
)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(fill = sig_dot), shape = 21, size = 3, color = "black", stroke = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "GOBP_NEUTROPHIL_MIGRATION",
    x = NULL,
    y = "CAMERA z.std",
    color = "Contrast Type",
    linetype = "Contrast Type",
    fill = "Significance"
  ) +
  scale_color_manual(values = c(
    "EE CON-unadjusted" = "#d95f02",
    "EE CON-adjusted" = "#d95f02",
    "Control" = "#7570b3"
  )) +
  scale_linetype_manual(values = c(
    "EE CON-unadjusted" = "dashed",
    "EE CON-adjusted" = "solid",
    "Control" = "dotdash"
  )) +
  scale_fill_manual(values = c(
    "Significant" = "black",
    "Not Significant" = "white"
  )) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )


###  Figure 3: CAMERA from EE vs. RE
# transcriptomics
camera_eere_t <- CAMERA_RESULTS %>%
  filter(tissue == "adipose",
         contrast_type == "Endur_vs_Resist",
         assay == 'transcript-rna-seq')
# Define the terms of interest
terms_of_interest <- c(
  "BIOCARTA_KREB_PATHWAY", "BIOCARTA_ETC_PATHWAY", "BIOCARTA_MTOR_PATHWAY",
  "PID_INTEGRIN1_PATHWAY", "PID_PI3KCI_AKT_PATHWAY", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM", "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_NEGATIVE_REGULATION_OF_MAPK_PATHWAY", "REACTOME_FATTY_ACID_METABOLISM",
  "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION", "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_MITOCHONDRIAL_BIOGENESIS", "WP_FATTY_ACID_BIOSYNTHESIS",
  "WP_PPAR_SIGNALING_PATHWAY", "WP_COMPLEMENT_AND_COAGULATION_CASCADES",
  "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_TRIGLYCERIDE_BIOSYNTHETIC_PROCESS",
  "GOBP_RIBOSOME_BIOGENESIS", "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
  "WP_INFLAMMATORY_RESPONSE_PATHWAY", "WP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_1_BETA_PRODUCTION",
  "GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION", "GOBP_REGULATION_OF_MACROPHAGE_ACTIVATION",
  "GOBP_B_CELL_MEDIATED_IMMUNITY", "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
  "GOCC_IMMUNOGLOBULIN_COMPLEX_CIRCULATING"
)

# Subset the data
camera_eere_t_filtered <- camera_eere_t %>%
  filter(set_short %in% terms_of_interest)

camera_eere_t_filtered <- camera_eere_t_filtered %>%
  mutate(
    contrast = case_when(
      contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_45min",
      contrast == "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_24_hr + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_24hr",
      contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_3.5_4_hr + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_4hr",
      TRUE ~ contrast
    ),
    contrast = factor(contrast, levels = c("EE_vs_RE_45min", "EE_vs_RE_4hr", "EE_vs_RE_24hr"))
  )

# Extract Exercise Type & Time for annotation
camera_eere_t_filtered <- camera_eere_t_filtered %>%
  mutate(
    Timepoint = case_when(
      grepl("45min", contrast) ~ "45minPost",
      grepl("4hr", contrast) ~ "4hrPost",
      grepl("24hr", contrast) ~ "24hrPost"
    )
  )
col_anno <- columnAnnotation(
  Time = c("45minPost", "4hrPost", "24hrPost"),
  col = list(
    Time = time_colors
  ),
  annotation_legend_param = list(
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))  # Set order explicitly
  )
)
# Figure 3D
camera_eere_t_plot <- camera_eere_t_filtered %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = FALSE,
    set_column = "set_short",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    colors = c("#1b9e77", "#d95f02"),
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      show_column_names = FALSE
    )
  )

# Proteomics 
camera_eere_pr <- CAMERA_RESULTS %>%
  filter(tissue == "adipose",
         contrast_type == "Endur_vs_Resist",
         assay == 'prot-pr') %>%
  filter(!database %in% c("KEGG_MEDICUS", "MITOCARTA"))
# Define the terms of interest
terms_of_interest <- c(
  "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
  "GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS",
  "GOCC_ACTIN_FILAMENT_BUNDLE",
  "GOCC_ACTOMYOSIN",
  "GOCC_CATALYTIC_STEP_1_SPLICEOSOME",
  "GOCC_CYTOPLASMIC_EXOSOME_RNASE_COMPLEX",
  "GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",
  "GOMF_ELECTRON_TRANSFER_ACTIVITY",
  "OXPHOS subunits",
  "Mitochondrial ribosome"
)

# Subset the data
camera_eere_pr_filtered <- camera_eere_pr %>%
  filter(set_short %in% terms_of_interest)

camera_eere_pr <- camera_eere_pr %>%
  mutate(contrast = case_when(
    contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_3.5_4_hr + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_4hr",
    TRUE ~ contrast  # Keep unchanged if no match
  ))

# Extract Exercise Type & Time for annotation
camera_eere_pr <- camera_eere_pr %>%
  mutate(
    Timepoint = case_when(
      grepl("4hr", contrast) ~ "4hrPost"
    )
  )
col_anno <- columnAnnotation(
  Time = c("4hrPost"),
  col = list(
    Time = time_colors
  ),
  annotation_legend_param = list(
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))  # Set order explicitly
  )
)

# Figure 3E
camera_eere_pr_plot <- camera_eere_pr %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = TRUE,
    set_column = "set_short",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    colors = c("#1b9e77", "#d95f02"),
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      show_column_names = FALSE
    )
  )

# Metabolomics
camera_eere_m <- CAMERA_RESULTS %>%
  filter(tissue == "adipose",
         contrast_type == "Endur_vs_Resist",
         assay == 'metab')

camera_eere_m <- camera_eere_m %>%
  mutate(contrast = case_when(
    contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_45min",
    contrast == "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_24_hr + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_24hr",
    contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_3.5_4_hr + group_timepointADUResist.pre_exercise" ~ "EE_vs_RE_4hr",
    TRUE ~ contrast),
    contrast = factor(contrast, levels = c("EE_vs_RE_45min", "EE_vs_RE_4hr", "EE_vs_RE_24hr"))
  )

# Extract Exercise Type & Time for annotation
camera_eere_m <- camera_eere_m %>%
  mutate(
    Timepoint = case_when(
      grepl("45min", contrast) ~ "45minPost",
      grepl("4hr", contrast) ~ "4hrPost",
      grepl("24hr", contrast) ~ "24hrPost"
    )
  )
col_anno <- columnAnnotation(
  Time = c("45minPost", "4hrPost", "24hrPost"),
  col = list(
    Time = time_colors
  ),
  annotation_legend_param = list(
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))  # Set order explicitly
  )
)
# Figure 3F
camera_eere_m_plot <- camera_eere_m %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = TRUE,
    set_column = "set_short",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    colors = c("#1b9e77", "#d95f02"),
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      column_title = NULL,
      show_column_names = FALSE
    )
  )


# Scatter plots to compare logFC between EE-CON vs. RE-CON per omic, per time point. 
# Run through line ~373 in precawg_adi_da.R to obtain all_vol
logfc_wide <- all_vol %>%
  filter(contrast_category %in% c("EE-CON", "RE-CON")) %>%
  mutate(group = case_when(
    contrast_category == "EE-CON" ~ "EE",
    contrast_category == "RE-CON" ~ "RE"
  )) %>%
  select(feature_id, assay, Timepoint, gene_symbol, group, logFC, adj_p_value) %>%
  pivot_wider(
    names_from = group,
    values_from = c(logFC, adj_p_value),
    names_sep = "_"
  ) %>%
  drop_na(logFC_EE, logFC_RE) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"),
                            labels = c("45min", "4hr", "24hr")))

# Step 2: Add significance category
logfc_wide <- logfc_wide %>%
  mutate(sig_group = case_when(
    adj_p_value_EE < 0.05 & adj_p_value_RE < 0.05 ~ "Both",
    adj_p_value_EE < 0.05 ~ "EE-only",
    adj_p_value_RE < 0.05 ~ "RE-only",
    TRUE ~ "NS"
  ))

# Figure S3
unique_assays <- unique(logfc_wide$assay)
for (a in unique_assays) {
  assay_df <- logfc_wide %>%
    filter(assay == a) %>%
    mutate(sig_group = case_when(
      adj_p_value_EE < 0.05 & adj_p_value_RE < 0.05 ~ "Both",
      adj_p_value_EE < 0.05 ~ "EE-only",
      adj_p_value_RE < 0.05 ~ "RE-only",
      TRUE ~ "NS"
    )) %>%
    mutate(sig_group = factor(sig_group, levels = c("NS", "RE-only", "EE-only", "Both")))
  
  plot <- ggplot(assay_df, aes(x = logFC_EE, y = logFC_RE)) +
    geom_point(data = subset(assay_df, sig_group == "NS"), aes(color = sig_group),
               alpha = 0.5, size = 1.5) +
    geom_point(data = subset(assay_df, sig_group != "NS"), aes(color = sig_group),
               alpha = 0.9, size = 1.8) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
    facet_wrap(~Timepoint, nrow = 1) +
    scale_color_manual(
      values = c("Both" = "#640D6B", "EE-only" = "#d95f02", "RE-only" = "#1b9e77", "NS" = "gray80")
    ) +
    coord_fixed() +
    labs(
      x = "logFC (EE-CON)",
      y = "logFC (RE-CON)",
      color = "Significance"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_blank(),  # remove title
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  # Save as EPS using cairo_ps
  ggsave(
    filename = paste0("logFC_EE_vs_RE_", gsub("[^A-Za-z0-9]", "_", a), ".eps"),
    plot = plot,
    width = 10, height = 4, device = cairo_ps
  )
}

###Export PTM-SEA results for Figure 2D
ptmsea_pai <- read.delim("n3_ptm-sea-results-combined.gct", skip = 2) 

contrast_map <- tibble(
  contrast = c(
    "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
    "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)"
  ),
  contrast_key = c("Endur", "Resist")
)
reshape_ptmsea <- function(df, contrast_map) {
  df %>%
    select(
      id,
      z.std_Endur   = `adipose.z.std_Endur.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.`,
      z.std_Resist  = `adipose.z.std_Resist.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.`,
      pval_Endur    = `pvalue.adipose.z.std_Endur.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.`,
      pval_Resist   = `pvalue.adipose.z.std_Resist.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.`,
      fdr_Endur     = `fdr.pvalue.adipose.z.std_Endur.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.`,
      fdr_Resist    = `fdr.pvalue.adipose.z.std_Resist.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.`
    ) %>%
    pivot_longer(
      cols = -id,
      names_to = c(".value", "Group"),
      names_pattern = "(z\\.std|pval|fdr)_(Endur|Resist)"
    ) %>%
    left_join(contrast_map, by = c("Group" = "contrast_key")) %>%
    mutate(
      tissue = "adipose",
      assay = "prot-ph",
      contrast_type = "exercise_with_controls"
    ) %>%
    select(tissue, assay, contrast_type, contrast, id, z.std, p_value = pval, adj_p_value = fdr)
}

ptmsea_pai_df <- reshape_ptmsea(ptmsea_pai, contrast_map)
ptmsea_pai_df_filtered <- ptmsea_pai_df %>%
  filter(!grepl("^(PERT|DISEASE)", id))

# Step 1: Filter for significant PTM-SEA results (adj p < 0.1)
ee_df <- ptmsea_pai_df %>%
  filter(grepl("^Endur", contrast)) %>%
  select(kinase = id, z.std_EE = z.std, p_EE = p_value, adj_p_EE = adj_p_value)

re_df <- ptmsea_pai_df %>%
  filter(grepl("^Resist", contrast)) %>%
  select(kinase = id, z.std_RE = z.std, p_RE = p_value, adj_p_RE = adj_p_value)

# Step 2: Join them by kinase
ptm_wide <- full_join(ee_df, re_df, by = "kinase")

# Step 3: Filtering (adj_p < 0.1)
ptm_wide <- ptm_wide %>%
  filter((!is.na(adj_p_EE) & adj_p_EE < 0.1) |
           (!is.na(adj_p_RE) & adj_p_RE < 0.1)) %>%
  mutate(
    log10p_EE = -log10(adj_p_EE),
    log10p_RE = -log10(adj_p_RE)
  ) %>%
  mutate(kinase_symbol = str_extract(kinase, "(?<=_)[^.]+"))
# ptm_wide was used to generate Figure 2D


## cytoskeleton story
sig_features <- ph_vol %>% # if all_vol was obtained, ph_vol should have been obtained as well. 
  filter(contrast_type == "exercise_with_controls",
         adj_p_value < 0.05) %>%
  pull(gene_symbol_with_phosphosite) %>%
  unique()
ph_vol_sig <- ph_vol %>%
  filter(contrast_type == "exercise_with_controls",
         gene_symbol_with_phosphosite %in% sig_features)

## cytoskeleton remodeling terms
cytoske_sets3 <- list(
  Cytoskeleton_Org = MOLECULAR_SIGNATURES$GOBP$GOBP_CYTOSKELETON_ORGANIZATION,
  Cell_motility =  MOLECULAR_SIGNATURES$GOBP$GOBP_CELL_MOTILITY,
  Cell_Matrix_Adhesion = MOLECULAR_SIGNATURES$GOBP$GOBP_CELL_MATRIX_ADHESION,
  Cell_junction_assembly =  MOLECULAR_SIGNATURES$GOBP$GOBP_CELL_CELL_JUNCTION_ASSEMBLY
)
gene_set_membership3 <- lapply(cytoske_sets3, function(gset) {
  unique(gset)
}) %>% 
  fromList()

upset(gene_set_membership3, order.by = "freq")

gene_to_pathway3 <- enframe(cytoske_sets3) %>%
  unnest_longer(value) %>%
  rename(Pathway = name, gene_symbol = value) %>%
  distinct()

pathway_priority3 <- c( "Cell_junction_assembly", "Cell_Matrix_Adhesion", "Cytoskeleton_Org", "Cell_motility")

gene_to_pathway_unique3 <- gene_to_pathway3 %>%
  mutate(Pathway = factor(Pathway, levels = pathway_priority3)) %>%
  arrange(gene_symbol, Pathway) %>%
  distinct(gene_symbol, .keep_all = TRUE)

# Join to DA table
ph_vol_sig_annot3 <- ph_vol_sig %>%
  left_join(gene_to_pathway_unique3, by = "gene_symbol")
ph_vol_sig_annot3 <- ph_vol_sig_annot3 %>%
  filter(base_feature_id != "O43399-2") # pivot wider goes into error with this feature because there are multiple rows for the same contrast. Since it's not significant, we remove it.
# Extract significance info per phosphosite and contrast
sig_status3 <- ph_vol_sig_annot3 %>%
  mutate(is_sig = adj_p_value < 0.05) %>%
  select(gene_symbol_with_phosphosite, contrast_category, is_sig) %>%
  pivot_wider(
    names_from = contrast_category,
    values_from = is_sig,
    names_prefix = "sig_"
  ) %>%
  mutate(
    significance = case_when(
      `sig_EE-CON` & `sig_RE-CON` ~ "Both",
      `sig_EE-CON` & !`sig_RE-CON` ~ "EE only",
      !`sig_EE-CON` & `sig_RE-CON` ~ "RE only",
      TRUE ~ "Not Sig"
    )
  )
ph_vol_sig_annot3 <- ph_vol_sig_annot3 %>%
  left_join(sig_status3, by = "gene_symbol_with_phosphosite")

cytoske_wide3 <- ph_vol_sig_annot3 %>%
  select(gene_symbol_with_phosphosite, Pathway, contrast_category, logFC, significance) %>%
  distinct() %>%
  pivot_wider(
    names_from = contrast_category,
    values_from = logFC,
    names_prefix = "logFC_"
  )
# Label top 10 per pathway
cytoske_labeled3 <- cytoske_wide3 %>%
  mutate(color_pathway = ifelse(is.na(Pathway), "Other", as.character(Pathway))) %>%
  group_by(Pathway) %>%
  arrange(desc(abs(`logFC_RE-CON`))) %>%
  mutate(label = if_else(
    row_number() <= 3 | 
      gene_symbol_with_phosphosite %in% c("TNS1-T1622", "TNS1-S1307;S1314"),
    gene_symbol_with_phosphosite,
    NA_character_
  )) %>%
  ungroup()

# Figure 2E
ggplot(cytoske_labeled3, aes(x = `logFC_EE-CON`, y = `logFC_RE-CON`)) +
  geom_point(
    aes(
      shape = significance,
      fill = color_pathway,
      alpha = ifelse(color_pathway == "Other", 0.3, 1)
    ),
    color = "black", size = 3, stroke = 0.4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_label_repel(
    data = cytoske_labeled3 %>% filter(color_pathway != "Other"),
    aes(label = label, color = color_pathway),
    fill = "white",
    box.padding = 0.8,
    point.padding = 0.5,
    size = 3.5,
    nudge_y = 0.5,
    segment.color = "gray50",
    label.size = 0.4,
    segment.curvature = -0.2,
    segment.angle = 20,
    segment.ncp = 3,
    force = 3,
    na.rm = TRUE,
    show.legend = FALSE,
    max.overlaps = 20
  ) +
  labs(
    x = "logFC (EE Post 4hr)",
    y = "logFC (RE Post 4hr)",
    title = "Phosphosite DA: Cytoskeleton remodeling",
    fill = "Pathway",
    shape = "Significance"
  ) +
  scale_fill_manual(values = c(
    "Cytoskeleton_Org" = "#e41a1c",
    "Cell_motility" = "#377eb8",
    "Cell_junction_assembly" = "#4daf4a",
    "Cell_Matrix_Adhesion" = "#984ea3",
    "Other" = "grey70"
  )) +
  scale_color_manual(values = c(
    "Cytoskeleton_Org" = "#e41a1c",
    "Cell_motility" = "#377eb8",
    "Cell_junction_assembly" = "#4daf4a",
    "Cell_Matrix_Adhesion" = "#984ea3"
  )) +
  scale_shape_manual(values = c(
    "Both" = 21,
    "EE only" = 24,
    "RE only" = 22
  )) +
  scale_alpha_identity() +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, alpha = 1)),
    shape = guide_legend()
  )
