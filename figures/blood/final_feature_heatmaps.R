feature_heatmap <- NULL # fix warning
source("figures/blood/function-feature_heatmap.R")
library(MotrpacHumanPreSuspension)
library(TMSig)
folder <- "figures/blood/feature_heatmaps"

## Fig 2D ----

# REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE...(12755)
feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_exercise_with_controls_REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE.pdf"
  ),
  post_min = 30,
  post_hr = 4
)

## Fig 2E ----
feature_heatmap(
  set_id = "12755",
  contrast_type = "exercise_with_controls",
  selected_tissue = "muscle",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "muscle_transcript-rna-seq_exercise_with_controls_REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE.pdf"
  )
)


## Fig S2B ----
genes_s2b <- c(
  "MFSD2A",
  "ABCD4",
  "SLC27A1",
  "SLC27A5",
  "FABP5",
  "SLC2A1",
  "ABCD2",
  "SLC27A4",
  "CD36",
  "SLC27A2",
  "ABCD1",
  "SLC27A3",
  "ABCD3"
)

features_s2b <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_s2b) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_s2b,
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_custom_Figure_S2B.pdf"
  )
)


## Fig 3E ----
genes_3e <- c(
  "AZU1", "GZMA", "GZMB", "IL2RA",
  "SELE", "FCER2", "CD70", "CD28",
  "CD22", "CD5", "CD4", "KLRD1",
  "CD38", "CD79B", "CD8A", "CD27"
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
  )
)


# Figure S3B ----
feature_heatmap(
  feature_ids = features_3e,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_exercise_with_controls_custom_Figure_S3B_B_cell.pdf"
  )
)


feature_heatmap(
  set_id = "16352",
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_memory_bcell_peripheral_blood.pdf"
  )
)


feature_heatmap(
  set_id = "15766",
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_bcell_peripheral_blood.pdf"
  )
)

##
feature_heatmap(
  set_id = "16496",
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_naive_bcell_blood.pdf"
  )
)

## Figure S3C ----
genes_s3c <- c(
  "FCGR2A", "CR2", "CD27", "CD79B",
  "CD79A", "MS4A1", "FCER2", "CD22",
  "PAX5"
)

features_s3c <- HUMAN_FEATURE_TO_GENE %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(gene_symbol %in% genes_s3c) %>%
  filter(assay == "transcript-rna-seq") %>%
  pull(feature_id)

feature_heatmap(
  feature_ids = features_s3c,
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_exercise_with_controls_custom_Figure_S3C.pdf"
  )
)


## Naive B Cell peripheral blood human ----
feature_heatmap(
  set_id = "16496",
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_naive_B_cell_peripheral.pdf"
  )
)

## B Cell Lymph Node Human ----
feature_heatmap(
  set_id = "15778",
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_B_cell_lymph_node.pdf"
  )
)


## Activated memory B Cell blood ----
feature_heatmap(
  set_id = "15740",
  contrast_type = "Endur_vs_Resist",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_EE_vs_RE_activated_memory_B_cell_blood.pdf"
  )
)


## B Cell Lymph Node Human ----
feature_heatmap(
  set_id = "15778",
  contrast_type = "exercise_with_controls",
  selected_tissue = "blood",
  selected_ome = "transcript-rna-seq",
  filename = file.path(
    folder,
    "blood_transcript-rna-seq_exercise_with_controls_B_cell_lymph_node.pdf"
  )
)



## Molecular Signature Bubble heatmaps ----

set_ids <- c(
  "15998", "16534", "16280", "15956",
  "16000", "15738", "15996", "16233",
  "15820", "16128", "16168", "15821",
  "16799", "16581", "16659", "16438",
  "15837", "16544", "16654", "16345",
  "16190", "16434", "16467", "16070",
  "15926", "16290", "16301", "16445",
  "16713", "16817", "16518", "16352",
  "15781", "15778", "16496", "16504"
)


plot_enrich_heatmap(
  x = CAMERA_RESULTS,
  set_ids = set_ids,
  selected_tissues = "blood",
  selected_ome = "transcript-rna-seq",
  filename = "figures/blood/feature_heatmaps/blood_transcript-rna-seq_exercise_with_controls_CELLMARKER_custom_bubble_heatmap.pdf"
)



## Figure S3A ----
set_ids <- c(
  "16496", "16352", "15781", "16504",
  "16654", "15814", "16287", "16236",
  "16475", "16070", "16301", "16467",
  "16190", "15926", "16434", "16581",
  "16128", "16445", "16817", "16713",
  "16233", "16168", "16544", "16653",
  "15956", "15738", "16000", "15996",
  "16534", "15998"
)

plot_enrich_heatmap(
  x = CAMERA_RESULTS,
  set_ids = set_ids,
  selected_tissues = "blood",
  selected_ome = "transcript-rna-seq",
  contrast_type = "Endur_vs_Resist",
  filename = "figures/blood/feature_heatmaps/blood_transcript-rna-seq_CELLMARKER_EE_vs_RE_CELLMARKER_custom_bubble_heatmap.pdf"
)

