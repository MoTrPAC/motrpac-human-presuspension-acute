library(here)

# Source the function file
source(here("updated_feature_heatmaps", "function-feature_heatmap_ver10.25.R"))

# Set output folder path
folder <- here("updated_feature_heatmaps")

# FIGURE S4A; RefMET All Sign

set_ids_FigS4A <- c(
  "13940", "13826", "13968", "13858", "13983", "13855", 
  "13951", "13907", "13979", "13900", "13902", "13879",  
  "13919", "13895", "13894", "13830", "13852", "13845", 
  "13863", "13828", "13969"
)

plot_enrich_heatmap(
  x = CAMERA_RESULTS,
  set_ids = set_ids_FigS4A,
  selected_tissues = "blood",
  selected_ome = "metab",
  filename = "updated_feature_heatmaps/FigureS4A_RefMet.pdf"
)

# FIGURE S4B METAB SET IDS: 13845 BLOOD
plot_feature_heatmap(
  set_id = "13845",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_13845.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE S4C METAB SET IDS: 13940 BLOOD
plot_feature_heatmap(
  set_id = "13940",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "metab",
  filename = file.path(
    folder,
    "blood_metab_13940.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)


# FIGURE 5A: SET 13826 MultiTissue BLOOD and MUSCLE
feature_heatmap(
  set_id = "13826",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("muscle","blood"),
  selected_ome = "metab",
  filename = file.path(
    folder,
    "Figure_5A_Blood_multitissue_13826_Acylcarnitines.pdf"
  ),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE 5B: MUSCLE TRANSCRIPTOMICS
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
    "Figure_5B_Blood_muscle_transcript_metabolism.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = TRUE
)



# FIGURE 6D: REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND
feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = c("blood", "muscle"), 
  selected_ome = "transcript",
  filename = file.path(
    folder,
    "Figure_6D_Blood_muscle_RNA_12755_REACTOME_HSP90_CHAPERONE.pdf"
  ),
  post_hr = 3.5,
  full_modality_names = FALSE
)

#Figure 6C: Blood RNA heatmap custom bubble heatmap
## Molecular Signature Bubble heatmaps Cellmarker----

set_ids_Fig6C <- c(
  "03075", "06183", "11552", "11535", "11618", "11712", "11620", "11605", 
  "11671", "11561", "11637", "00068", "00269", "13168", "12755", "12753", 
  "13025", "12783", "13491", "12560", "13370", "12614", "13539", "12329", 
  "13402", "13448", "13336", "12574", "14621", "14113", "14339", "14183",
  "14541", "14658", "14717", "14138", "14331", "11544", "11496", "10874"
)

plot_enrich_heatmap(
  x = CAMERA_RESULTS,
  set_ids = set_ids_Fig6C,
  selected_tissues = "blood",
  selected_ome = "transcript-rna-seq",
  filename = "updated_feature_heatmaps/Figure6C_Blood_top_curated_RNA.pdf"
)

#FIGURE S6B blood_RNA_exercise_REACTOME_HSF1_ACTIVATION
feature_heatmap(
  set_id = "12753",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "Figure_S6B_blood_transcript_REACTOME_HSF1_ACTIVATION.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)



#Figure 7A: CELLMARKER custom bubble heatmap
## Molecular Signature Bubble heatmaps Cellmarker----

set_ids_Fig7A <- c(
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

plot_enrich_heatmap(
  x = CAMERA_RESULTS,
  set_ids = set_ids_Fig7A,
  selected_tissues = "blood",
  selected_ome = "transcript-rna-seq",
  filename = "updated_feature_heatmaps/Figure7A_CELLMARKER_RNA.pdf"
)


# FIGURE 7B: BLOOD PROTEOMICS
genes_7b = c("AZU1", "GZMA", "GZMB", "IL2RA", "TNFRSF4", "CD209", "CD207", 
             "LY75", "CLEC7A", "CD244", "EPO", "CD63", "IL17RA", "NCR1", "ITGAM", 
             "IL1RL1", "IL2RA", "IL17RB", "CD99", "CD274", "MMP8", "MMP9", 
             "S100A12", "CSF2RA", "SELE", "FCER2", "CD70", "CD28", "IFNG", 
             "IL7R", "CCL3", "CCL4", "TNFRSF8", "CD22", "CD5", "CD4", "KLRD1", 
             "CEACAM8", "GNLY", "CD74", "TNFRSF9", "CD38", "CD79B", "CD8A", 
             "CD27", "CD69", "KLRC1", "LY9"
)

features_7b <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_7b) %>%
  filter(assay == "prot-ol") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_7b,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "prot-ol",
  filename = file.path(
    folder,
      "Figure7B_blood_prot-ol.pdf"
  ),
  post_min = 30,
  post_hr = 3.5
)

#Figure S7A_CELLMARKER_Eosinophil Blood Human
feature_heatmap(
  set_id = "16070",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "Figure_S7Ablood_transcript_CellMarker CELLMARKER_Eosinophil.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

#Figure_S7B_CELLMARKER_Granulocytes
feature_heatmap(
  set_id = "16190",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "Figure_S7B_CELLMARKER_blood_transcript_16190_granulocyte.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE_7C TRANSCRIPT SUBSET BLOOD monocyte/macrophage (curated)
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

feature_heatmap(
  feature_ids = features_monocyte,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "Figure_7C_blood_transcript_monocyte_macrophage_curated.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE_7D_NK GS curation
genes_NK_GS <- c(
  "GZMA", "GZMB", "B3GAT1", "CX3CR1", "KLRD1", "KLRK1", "KLRC1", "KLRG1",
  "GZMH", "NCAM1", "CD56")
features_NK_GS <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_NK_GS) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_NK_GS,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "FIGURE_7D_blood_transcript_NK_GS_curated.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)

#Figure_7F_CELLMARKER_Platelet Peripheral Blood Human
feature_heatmap(
  set_id = "16654",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "Figure_7Fblood_transcript_CELLMARKER_Platelet_Peripheral_Blood.pdf"),
  post_min = 15,
  post_hr = 3.5,
  full_modality_names = FALSE
)

# FIGURE_7E_TRANSCRIPT SUBSET BLOOD T cell markers GUILLAUME
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

feature_heatmap(
  feature_ids = features_Tcell,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "FIGURE_7E_blood_transcript_Tcells_curated_GS.pdf"),
  post_min = 30,
  post_hr = 3.5,
  full_modality_names = FALSE
)
