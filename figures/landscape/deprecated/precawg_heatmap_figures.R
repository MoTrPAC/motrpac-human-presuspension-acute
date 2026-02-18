library(here)
here::i_am("figures/precawg_heatmap_figures.R")
source(here("figures","blood","function-feature_heatmap.R"))
folder = here("figures","precawg_figures_dl")


# FIGURE: MUSCLE CURATED METAB LIST
metabs_a <- c("ATP", "ADP", "AMP", "CMP", "CDP", "CTP", "NAD+", "NADH", "IMP",
              "NADP+", "NADPH", "CMP", "CDP", "CTP", "UMP", "UDP", "UTP")
feature_heatmap(
  feature_ids = metabs_a,
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "metab",
  # filename = file.path(
  #   folder,
  #   "muscle_metab_exercise_with_controls_energy_v2.pdf"
  # ),
  post_min = 15,
  post_hr = 3.5
)

# FIGURE: SET ID 13826 MUSCLE
feature_heatmap(
  set_id = "13826",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "muscle_set_id_13826.pdf"
  ),
  post_min = 15,
  post_hr = 3.5
)

# FIGURE: SET ID 13826 BLOOD
feature_heatmap(
  set_id = "13826",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_set_id_13826.pdf"
  ),
  post_min = 15,
  post_hr = 3.5
)

# FIGURE: SET 13826 ADIPOSE
feature_heatmap(
  set_id = "13826",
  contrast_type = "exercise_with_controls",
  selected_tissue = "adipose",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "adipose_set_id_13826.pdf"
  ),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE: SET 13826 MultiTissue BLOOD and MUSCLE
feature_heatmap(
  set_id = "13826",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood"),
  selected_ome = "metab",
  filename = file.path(
    folder,
    "multitissue_set_id_13826.pdf"
  ),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# note: the goal would to be to combine the above ones into a single figure

# FIGURE: BLOOD PROTEOMICS
genes_3e = c("AZU1", "GZMA", "GZMB", "IL2RA", "TNFRSF4", "CD209", "CD207", "LY75", "CLEC7A", "CD244", "EPO", "CD63", "IL17RA", "NCR1", "ITGAM", "IL1RL1", "IL2RA", "IL17RB", "CD99", "CD274", "MMP8", "MMP9", "S100A12", "CSF2RA",
"SELE", "FCER2", "CD70", "CD28", "IFNG", "IL7R", "CCL3", "CCL4", "TNFRSF8",
"CD22", "CD5", "CD4", "KLRD1", "CEACAM8", "GNLY", "CD74", "TNFRSF9",
"CD38", "CD79B", "CD8A", "CD27", "CD69", "KLRC1", "LY9"
)

features_3e <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_3e) %>%
  filter(assay == "prot-ol") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_3e,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "prot-ol",
  filename = file.path(
    folder,
    "blood_prot-ol_exercise_with_controls_custom_Figure_3E.pdf"
  ),
  post_min = 15,
  post_hr = 3.5
)

# FIGURE: MUSCLE TRANSCRIPTOMICS
values = c("CPT1B", "SLC25A20","CPT1A","ACADVL","CPT2",
           "SLC22A4","ACADM","ACADS","CRAT","CROT",
           "ACADL","CPT1C")
features_mt <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% values) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_mt,
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = TRUE
  )

# FIGURE: BLOOD SET ID 12753
feature_heatmap(
  set_id = "12753",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_12753.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = TRUE
)

# FIGURE: BLOOD SET ID 12755
feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_12755.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = TRUE
)

# FIGURE: MUSCLE SET ID 12755
feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_12755.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = TRUE
)

# FIGURE: MUSCLE and BLOOD SET ID 12755
feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood"),
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "multitissue_transcript_12755.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE: TRANSCRIPT SET ID 14323
feature_heatmap(
  set_id = "14323",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_14323.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID 16070
feature_heatmap(
  set_id = "16070",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16070.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE: RNA MUSCLE
# "CPT1B", "CPT1A", "SLC25A20", "ACADVL", "ACADL",
# "ACADM", "CRAT", "CROT", "CPT2", "SLC22A4", "CPT1C", "LIPE",
# "PNPLA2", "MGLL", "CD36", "FABP4", "FABP5", "FABP3"

values = c("CPT1B", "CPT1A", "SLC25A20", "ACADVL", "ACADL",
  "ACADM", "CRAT", "CROT", "CPT2", "SLC22A4", "CPT1C", "LIPE",
  "PNPLA2", "MGLL", "CD36", "FABP4", "FABP5", "FABP3")

features_mt <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% values) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_mt,
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_cpt1b_cpt1a_etc.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE METAB SET IDS: 13845 BLOOD
feature_heatmap(
  set_id = "13845",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_13845.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE METAB SET IDS: 13979 BLOOD
feature_heatmap(
  set_id = "13979",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_13979.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE METAB SET IDS: 13940 BLOOD
feature_heatmap(
  set_id = "13940",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_13940.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 15997 BLOOD
feature_heatmap(
  set_id = "15997",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15997.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 16354 BLOOD
feature_heatmap(
  set_id = "16354",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16354.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 15996 BLOOD
feature_heatmap(
  set_id = "15996",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15996.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 15837 BLOOD
feature_heatmap(
  set_id = "15837",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15837.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 15995 BLOOD
feature_heatmap(
  set_id = "15995",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15995.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 16504 BLOOD
feature_heatmap(
  set_id = "16504",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16504.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 16503 BLOOD
feature_heatmap(
  set_id = "16503",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16503.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 16654 BLOOD
feature_heatmap(
  set_id = "16654",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16654.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 16653 BLOOD
feature_heatmap(
  set_id = "16653",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16653.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS: 16534 BLOOD
feature_heatmap(
  set_id = "16534",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16534.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SET IDS 11544 BLOOD
feature_heatmap(
  set_id = "11544",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_11544.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)


# FIGURE TRANSCRIPT SUBSET BLOOD macrophage subset
genes_m1 <- c("NOS2","IRF5","STAT1","IFNG","HLA-DMB",
              "IL1B","CD86","CLEC7A","FCGR2A","FCGR3A","TLR4",
              "HLA-DRB1","HLA-DRA","TNF","CD74","HLA-DQA1"
              )

features_m1 <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_m1) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_m1,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_macrophage.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SUBSET BLOOD monocyte markers
genes_mo <- c(
  "CD14",
  "CCR2",
  "FCGR3B",
  "ITGAM",
  "FCGR3A",
  "HLA-DRB1",
  "HLA-DRA",
  "HLA-DQA1",
  "CD68",
  "CX3CR1"
)
features_mo <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_mo) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_mo,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_monocyte.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)


feature_heatmap(
  set_id = "12753",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("adipose","blood","muscle"),
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "all_tissues_12753.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

