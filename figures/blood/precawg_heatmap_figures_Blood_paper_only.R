library(here)
library(ComplexHeatmap)
library(gridtext)
library(dplyr)
library(grDevices)
library(tibble)
library(tidyr)
library(TMSig) # enrichmap
library(MotrpacHumanPreSuspension)

# folder for blood figure heatmaps
folder = here("figures","blood","feature_heatmaps")

# FIGURE: MUSCLE CURATED METAB LIST
metabs_a <- c("ATP", "ADP", "AMP", "CMP", "CDP", "CTP", "NAD+", "NADH", "IMP",
              "NADP+", "NADPH", "CMP", "CDP", "CTP", "UMP", "UDP", "UTP")
plot_feature_heatmap(
  feature_ids = metabs_a,
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "muscle_metab_exercise_with_controls_energy.pdf"
  ),
  post_min = 15,
  post_hr = 3.5
)

# FIGURE: SET ID 13826 MUSCLE
plot_feature_heatmap(
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

# FIGURE: MUSCLE CURATED METAB LIST 2
metabs_2 <- c("Methionine", "Serine", "Threonine",
              "2-Hydroxy-3-methylbutyric acid", "2-Hydroxybutyric acid",
              "3-Hydroxybutyric acid", "2-Oxoglutaric acid",
              "2-Hydroxyglutaric acid", "Valine",
              "Glutamic acid", "Isoleucine", "Leucine", "N-Acetylaspartic acid")
plot_feature_heatmap(
  feature_ids = metabs_2,
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "muscle_metab_exercise_with_controls_list2.pdf"
  ),
  post_min = 15,
  post_hr = 3.5
)

#muscle REFMET_Hydroxy FA
plot_feature_heatmap(
  set_id = "13879",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "muscle_metab_REFMET_Hydroxy FA.pdf"
  ),
  post_min = 15,
  post_hr = 3.5
)

# FIGURE: SET ID 13826 MUSCLE
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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

plot_feature_heatmap(
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

plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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

plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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
plot_feature_heatmap(
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

plot_feature_heatmap(
  feature_ids = features_m1,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_macrophage.pdf"),
  post_min = 15,
  post_hr = 4,
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

plot_feature_heatmap(
  feature_ids = features_mo,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_monocyte.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)


# FIGURE TRANSCRIPT SUBSET BLOOD monocyte/macrophage
genes_monocyte <- c("NOS2","IRF5","STAT1","IFNG","HLA-DMB", "IL1B",
                    "CLEC7A","FCGR2A","FCGR3A","TLR4",
                    "HLA-DRB1","HLA-DRA","TNF","CD74","HLA-DQA1", "CD14",
                    "CCR2", "FCGR3B", "ITGAM", "FCGR3A", "HLA-DRB1",
                    "HLA-DRA", "HLA-DQA1", "CD68", "CX3CR1", "IL10",
                    "ARG1", "CHI3L1", "TGFB1", "IL1RN", "SOCS3",
                    "CD14", "TLR2", "TLR4", "TNF", "FCGR1A"
)

features_monocyte <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_monocyte) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_monocyte,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_monocyte_GS_curated.pdf"),
  post_min = 30,
  post_hr = 4,
  full_modality_names = FALSE
)



# FIGURE: SET 12753 MultiTissue (ALL) RNA
plot_feature_heatmap(
  set_id = "12753",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("adipose", "blood", "muscle"),
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "multitissue_set_id_12753.pdf"
  ),
  full_modality_names = FALSE
)


# FIGURE: BLOOD SET ID 16190
plot_feature_heatmap(
  set_id = "16190",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16190_granulocyte.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID 16190
plot_feature_heatmap(
  set_id = "16190",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16190_granulocyte.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID 16434
plot_feature_heatmap(
  set_id = "16434",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16434_monocyte_periph_blood human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CELLMARKER_CD4+ T Cell Blood Human
plot_feature_heatmap(
  set_id = "15823",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15823_CELLMARKER_CD4+ T Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)



# FIGURE: BLOOD SET ID CELLMARKER_Activated CD4+ T Cell Blood Human
plot_feature_heatmap(
  set_id = "15736",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CELLMARKER_Activated CD4+ T Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)


# FIGURE: BLOOD SET ID CELLMARKER_Naive CD4 T Cell Blood
plot_feature_heatmap(
  set_id = "16499",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16499_CELLMARKER_Naive CD4 T Cell Blood.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CELLMARKER_B Cell Blood Human
plot_feature_heatmap(
  set_id = "15766",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15766_CELLMARKER_B Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CELLMARKER_B Cell Blood Human
plot_feature_heatmap(
  set_id = "15766",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_15766_CELLMARKER_B Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)



# FIGURE: BLOOD SET ID CELLMARKER_Naive CD8+ T Cell Peripheral Blood Human
plot_feature_heatmap(
  set_id = "16511",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_16511_Naive CD8+ T Cell Peripheral Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)



# FIGURE: BLOOD SET ID CELLMARKER_Central Memory CD8+ T Cell Blood Human
plot_feature_heatmap(
  set_id = "15911",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_Central Memory CD8+ T Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)


# FIGURE: BLOOD SET ID CellMarker CD4+ Central Memory Like T (Tcm-like)
#Cell Peripheral Blood Human
plot_feature_heatmap(
  set_id = "15819",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker CD4+ Central Memory Like T (Tcm-like).pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CellMarker CD4+ Recently Activated Effector Memory Or Effector T Cell (CTL) Blood Human
plot_feature_heatmap(
  set_id = "15822",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker CD4+ Recently Activated Effector Memory Or Effector
    T Cell (CTL) Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CellMarker Central Memory CD8+ T Cell Blood Human
plot_feature_heatmap(
  set_id = "15911",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker Central Memory CD8+ T Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CellMarker Effector CD8+ T Cell Blood Human
plot_feature_heatmap(
  set_id = "15999",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker Effector CD8+ T Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)


# FIGURE: BLOOD SET ID CellMarker Memory CD8+ T Cell Blood Human
plot_feature_heatmap(
  set_id = "16357",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker Memory CD8+ T Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)


# FIGURE: BLOOD SET ID CellMarker Exhausted CD8+ T Cell Peripheral Blood Human
plot_feature_heatmap(
  set_id = "16122",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker Exhausted CD8+ T Cell Peripheral Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)



# FIGURE: BLOOD SET ID CellMarker Natural Killer Cell Blood Human
plot_feature_heatmap(
  set_id = "16122",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker Natural Killer Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)


# FIGURE: BLOOD SET ID Natural Killer T (NKT) Cell Peripheral Blood Human
plot_feature_heatmap(
  set_id = "16544",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker Natural Killer T (NKT) Cell Peripheral Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CELLMARKER_Dendritic Cell Blood Human
plot_feature_heatmap(
  set_id = "15966",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker CELLMARKER_Dendritic Cell Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE: BLOOD SET ID CELLMARKER_Plasmacytoid Dendritic cell(pDC) Blood Human
plot_feature_heatmap(
  set_id = "16645",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker CELLMARKER_Plasmacytoid Dendritic cell(pDC) Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#Dendritic Cell Blood Human
plot_feature_heatmap(
  set_id = "15966",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker CELLMARKER_Dendritic cell Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#CellMarker M1 Macrophage Peripheral Blood Human
plot_feature_heatmap(
  set_id = "16287",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker M1 Macrophage Peripheral Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#CellMarker M2 Macrophage Peripheral Blood Human
plot_feature_heatmap(
  set_id = "16298",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker M2 Macrophage Peripheral Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SUBSET BLOOD T cell markers GUILLAUME
genes_Tcell <- c(
  "KLRG1", "PTK7", "CD69","CCR4", "CCR5", "CCR6",
  "PECAM1", #CD31
  "PTPRC", #CD45RA
  "CD27", "CD28", "CD44",
  "B3GAT1", #cd57
  "CCR7", "CTLA4", "IL7R",
  "SLC3A2", "SLC7A5", #CD98
  "PDCD1", #PD1
  "CDKN2A", "CDKN1A",
  "CD3E", "CD3G", "CD3D"
)
features_Tcell <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_Tcell) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_Tcell,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_Tcells_curated.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)



# FIGURE TRANSCRIPT SUBSET BLOOD CD4+ EM +GZMK
genes_CD4EF <- c(
  "GZMK", "CD4", "GNLY", "CX3CR1", "TBX21", "PRF1", "NKG7", "S1PR1",
  "GZMH", "CTSW", "S1PR5", "KLRG1"
)
features_CD4EF <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_CD4EF) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_CD4EF,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CD4+_Effector Memory_curated.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# FIGURE NK GS curation
genes_NK_GS <- c(
  "GZMA", "GZMB", "B3GAT1", "CX3CR1", "KLRD1", "KLRK1", "KLRC1", "KLRG1",
  "GZMH", "NCAM1", "CD56")
features_NK_GS <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_NK_GS) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_NK_GS,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_NK_GS_curated.pdf"),
  post_min = 30,
  post_hr = 4,
  full_modality_names = FALSE
)

#REFMET_PI blood METAB
plot_feature_heatmap(
  set_id = "13919",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_REFMET_PI.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#REFMET_C24 Bile acids blood METAB
plot_feature_heatmap(
  set_id = "13845",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_REFMET_C24 Bile acids.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#REFMET_C24 Bile acids METAB
plot_feature_heatmap(
  set_id = "13845",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_muscle_REFMET_C24 Bile acids.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#REFMET_Acyl Carnitines METAB
plot_feature_heatmap(
  set_id = "13826",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood", "muscle",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_muscle_REFMET_Acyl Carnitiness.pdf"),
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SUBSET BLOOD CD4+ GS curation
genes_CD4_GS <- c(
  "KLRB1", "ZFP36", "SELL", "ZFP36L1", "PTPRC", "TRIM21", "IDH1",
  "TNFSF13", "STAT2", "ERAP1", "ERAP2", "IL32", "OAS1", "ZNRF3",
  "CD69", "CD4",  "IL7R", "CD27", "CCR7", "IRF8"
)
features_CD4_GS <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_CD4_GS) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_CD4_GS,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CD4+_GS_curated.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE TRANSCRIPT SUBSET BLOOD CD8 GS curation
genes_CD8_GS <- c(
  "GZMB",	"CX3CR1",	"GPR183", "PRF1", "S100A4", "CD28", "CD27", "B3GAT1", "CCR7",
  "HAVCR2", "SELL", "GZMK",	"FGFBP2",	"GNLY",	"CD8A",	"CXCR4",	"CD38",	"CD69",
  "GZMH",	"KLRG1",	"IL7R"
)
features_CD8_GS <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_CD8_GS) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_CD8_GS,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CD8+_GS_curated.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

#B CELL GS CURATED
genes_B_GS <- c(
  "AICDA", "CD52",	"MME",	"PAX5",	"CDC20", "PTPRC",	"BANK1",	"BCL6",	"CDK1",
  "CD24", "FCER2",	"MS4A1",	"CD22",	"MKI67",	"CD1C",	"JCHAIN",	"NR4A1",	"TNFSF13B",
  "CD79B",	"CD79A",	"FCRLA",	"TCL1A", "VPREB3",	"CXCR4",	"CD19",	"BLNK",
  "CD38",	"JUNB",	"CD99"
)
features_B_GS <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_B_GS) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

plot_feature_heatmap(
  feature_ids = features_B_GS,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_B_CELL_GS_curated.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)


# FIGURE: Blood REFMET_Xanthines
plot_feature_heatmap(
  set_id = "13983",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_set_id_13983_Xanthines.pdf"
  ),
  post_min = 30,
  post_hr = 3.5
)


# FIGURE: Blood REFMET_Purine ribonucleosides
plot_feature_heatmap(
  set_id = "13940",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_set_id_13940_REFMET_Purine ribonucleosides.pdf"
  ),
  post_min = 30,
  post_hr = 3.5
)


#All figures across tissues Landscape Figure 4
# FIGURE: SET 13826 MultiTissue BLOOD and MUSCLE
#PID_PDGFRB_PATHWAY
plot_feature_heatmap(
  set_id = "11663",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_11663_PID_PDGFRB_PATHWAY.pdf"
  ),
  full_modality_names = FALSE
)

#WP_VEGFA_VEGFR2_SIGNALING
plot_feature_heatmap(
  set_id = "14756",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_14756_WP_VEGFA_VEGFR2_SIGNALING
.pdf"
  ),
  full_modality_names = FALSE
)

#REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND

plot_feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_12755_REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND
.pdf"
  ),
  full_modality_names = FALSE
)

#blood_RNA_exercise_REACTOME_HSF1_ACTIVATION.pdf
plot_feature_heatmap(
  set_id = "12753",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_RNA_exercise_REACTOME_HSF1_ACTIVATION.pdf"),
  post_min = 30,
  post_hr = 4,
  full_modality_names = FALSE
)

#REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND blood

plot_feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("blood", "muscle"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "blood_muscle_RNA_12755_REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND
.pdf"
  ),
  full_modality_names = FALSE
)


#REACTOME_AUTOPHAGY
plot_feature_heatmap(
  set_id = "12259",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_12259_REACTOME_AUTOPHAGY.pdf"
  ),
  full_modality_names = FALSE
)

#PID_IL6_7_PATHWAY

plot_feature_heatmap(
  set_id = "11618",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_11618_PID_IL6_7_PATHWAY.pdf"
  ),
  full_modality_names = FALSE
)

#CELLMARKER_Mesenchymal Stem Cell Undefined Human
plot_feature_heatmap(
  set_id = "16399",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_16399_CELLMARKER_Mesenchymal Stem Cell Undefined Human.pdf"
  ),
  full_modality_names = FALSE
)

#REACTOME_CIRCADIAN_CLOCK
plot_feature_heatmap(
  set_id = "12358",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_12358_REACTOME_CIRCADIAN_CLOCK.pdf"
  ),
  full_modality_names = FALSE
)

#WP_KEGG_MEDICUS_REFERENCE_MICROTUBULE_RHOA_SIGNALING_PATHWAY
plot_feature_heatmap(
  set_id = "11096",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_KEGG_MEDICUS_REFERENCE_MICROTUBULE_RHOA_SIGNALING_PATHWAY.pdf"
  ),
  full_modality_names = FALSE
)


#PID_ANGIOPOIETIN_RECEPTOR_PATHWAY
plot_feature_heatmap(
  set_id = "11528",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_11528_PID_ANGIOPOIETIN_RECEPTOR_PATHWAY.pdf"
  ),
  full_modality_names = FALSE
)

#PID_ERBB1_DOWNSTREAM_PATHWAY
plot_feature_heatmap(
  set_id = "11575",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_11575_PID_ERBB1_DOWNSTREAM_PATHWAY.pdf"
  ),
  full_modality_names = FALSE
)

#WP_TGF_BETA_SIGNALING_PATHWAY
plot_feature_heatmap(
  set_id = "14704",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("adipose", "blood", "muscle"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_WP_TGF_BETA_SIGNALING_PATHWAY.pdf"
  ),
  full_modality_names = FALSE
)

#WP_MFAP5_EFFECT_ON_PERMEABILITY_AND_MOTILITY_OF...(14444)
plot_feature_heatmap(
  set_id = "14444",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("adipose", "blood", "muscle"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_MFAP5_EFFECT_ON_PERMEABILITY_AND_MOTILITY_OF.pdf"
  ),
  full_modality_names = FALSE
)

#GOBP_PURINE_NUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHETIC_PROCESS

plot_feature_heatmap(
  set_id = "05487",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("adipose", "muscle"),
  selected_ome = "prot-pr",
  filename = file.path(
    folder,
    "multitissue_GOBP_PURINE_NUCLEOSIDE_MONOPHOSPHATE_BIOSYNTHETIC_PROCESSS.pdf"
  ),
  full_modality_names = FALSE
)

#MITOCARTA_OXPHOS subunits
plot_feature_heatmap(
  set_id = "11463",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood", "adipose"),
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "multitissue_set_id_11463_MITOCARTA_OXPHOS subunits.pdf"
  ),
  full_modality_names = FALSE
)

#CELLMARKER_Eosinophil Blood Human
plot_feature_heatmap(
  set_id = "16070",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CellMarker CELLMARKER_Eosinophil Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

#CELLMARKER_Platelet Peripheral Blood Human
plot_feature_heatmap(
  set_id = "16654",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript_CELLMARKER_Platelet Peripheral Blood Human.pdf"),
  post_min = 15,
  post_hr = 4,
  full_modality_names = FALSE
)

# Muscle WP_UNFOLDED_PROTEIN_RESPONSE
plot_feature_heatmap(
  set_id = "14749",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_WP_UNFOLDED_PROTEIN_RESPONSE
.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)


## Molecular Signature Bubble heatmaps Cellmarker----

set_ids <- c(
  "15998", "16534", "16280", "15956",
  "16000", "15738", "15996", "16233",
  "15820", "16128", "16168", "15821",
  "16799", "16581", "16659", "16438",
  "15837", "16544", "16654", "16345",
  "16190", "16434", "16467", "16070",
  "15926", "16290", "16301", "16445",
  "16713", "16817", "16518", "16352",
  "15781", "15778", "16496", "16504",
  "15740"
)

#MITOCARTA_OXPHOS
plot_feature_heatmap(
  set_id = "11461",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_WP_UNFOLDED_PROTEIN_RESPONSE
.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# device.off()
# plot_enrich_heatmap(
#   x = CAMERA_RESULTS,
#   set_ids = set_ids,
#   selected_tissues = "blood",
#   selected_ome = "transcript-rna-seq",
#   filename = "blood/bubbleheatmap_cellmarker_RNA_blood_10.24.25.pdf"
# )

#GOMF_NUCLEAR_GLUCOCORTICOID_RECEPTOR_BINDING
plot_feature_heatmap(
  set_id = "10026",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_GOMF_NUCLEAR_GLUCOCORTICOID_RECEPTOR_BINDING.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

#WP_GLUCOCORTICOID_RECEPTOR_PATHWAY
plot_feature_heatmap(
  set_id = "14272",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript_GLUCOCORTICOID_RECEPTOR_PATHWAY.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# #Atrogenes muscle
# genes_atrophy <- c("Trim32", "Trim28")
# features_atrophy <- HUMAN_FEATURE_TO_GENE %>%
#   mutate(across(.cols = everything(),
#                 .fns = as.character)) %>%
#   filter(gene_symbol %in% genes_atrophy) %>%
#   filter(assay == "transcript-rna-seq") %>%
#            pull(feature_id)
#
#  plot_feature_heatmap(
#    feature_ids = features_atrophy,
#    contrast_type = "exercise_with_controls",
#    selected_tissue = "muscle",
#    selected_ome = "transcript-rna-seq",
#    filename = file.path(
#      folder,
#      "muscle_transcript_atrophy_curated.pdf"),
#    post_min = 15,
#    post_hr = 3.5,
#    full_modality_names = FALSE
#  )

 # All tissues 12753 set id
plot_feature_heatmap(
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
