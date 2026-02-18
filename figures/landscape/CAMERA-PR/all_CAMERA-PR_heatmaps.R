## Author: Tyler Sagendorf
## Date: 2025-05-07
##
## Purpose: Create bubble heatmaps of the CAMERA-PR results for all combinations
## of tissue, ome, collection, and contrast type. Includes multi-tissue
## heatmaps.

library(MotrpacHumanPreSuspensionAnalysis)
library(dplyr)

## Folder to save plots, relative to precovid-analyses
folder <- file.path("figures", "landscape", "CAMERA-PR")

if (!dir.exists(file.path(folder, "plots_single_tissue"))) {
  dir.create(file.path(folder, "plots_single_tissue"))
}

if (!dir.exists(file.path(folder, "plots_multi_tissue"))) {
  dir.create(file.path(folder, "plots_multi_tissue"))
}

## Single tissue heatmaps ----
comb_df <- CAMERA_RESULTS %>%
  distinct(tissue, assay, collection, contrast_type) %>%
  mutate(across(.cols = everything(),
                .fns = as.character))

for (i in seq_len(nrow(comb_df))) {
  tissue_i <- comb_df[i, "tissue"]
  assay_i <- comb_df[i, "assay"]
  collection_i <- comb_df[i, "collection"]
  contrast_type_i <- comb_df[i, "contrast_type"]

  # 10% FDR used for phospho results, 5% FDR used for everything else
  padj_cutoff <- 0.05 + 0.05 * (assay_i == "prot-ph")

  xi <- CAMERA_RESULTS %>%
    filter(tissue == tissue_i,
           assay == assay_i,
           collection == collection_i,
           contrast_type == contrast_type_i) %>%
    filter(.by = set,
           any(adj_p_value < padj_cutoff))

  # Separate folders for each contrast type
  if (!dir.exists(file.path(folder, "plots_single_tissue", contrast_type_i)))
    dir.create(file.path(folder, "plots_single_tissue", contrast_type_i))

  if (nrow(xi) > 0L) {
    plot_enrich_heatmap(
      x = xi,
      # Top 15 most significant terms from each contrast
      n_top = 15L,
      selected_ome = assay_i,
      selected_tissues = tissue_i,
      contrast_type = contrast_type_i,
      padj_cutoff = padj_cutoff,
      filename = file.path(
        folder, "plots_single_tissue", contrast_type_i,
        sprintf("CAMERA-PR_heatmap_%s_%s_%s_%s.pdf",
                contrast_type_i, tissue_i, assay_i, collection_i)
      )
    )
  }
}


## Multi-tissue plots ----
comb_df_multi <- comb_df %>%
  # Keep omes measured in multiple tissues
  filter(.by = c(assay, collection, contrast_type),
         n() > 1L) %>%
  select(-tissue) %>%
  distinct()

for (i in seq_len(nrow(comb_df_multi))) {
  tissue_i <- comb_df_multi[i, "tissue"]
  assay_i <- comb_df_multi[i, "assay"]
  collection_i <- comb_df_multi[i, "collection"]
  contrast_type_i <- comb_df_multi[i, "contrast_type"]
  padj_cutoff <- 0.05 + 0.05 * (assay_i == "prot-ph")

  xi <- CAMERA_RESULTS %>%
    filter(assay == assay_i,
           collection == collection_i,
           contrast_type == contrast_type_i) %>%
    droplevels.data.frame() %>%
    filter(.by = set,
           any(adj_p_value < padj_cutoff))

  # Separate folders for each contrast type
  if (!dir.exists(file.path(folder, "plots_multi_tissue", contrast_type_i)))
    dir.create(file.path(folder, "plots_multi_tissue", contrast_type_i))

  if (nrow(xi) > 0L) {
    plot_enrich_heatmap(
      x = xi,
      # Top 15 most significant terms from each contrast
      n_top = 15L,
      selected_ome = assay_i,
      contrast_type = contrast_type_i,
      padj_cutoff = padj_cutoff,
      filename = file.path(
        folder, "plots_multi_tissue", contrast_type_i,
        sprintf("CAMERA-PR_heatmap_%s_multi-tissue_%s_%s.pdf",
                contrast_type_i, assay_i, collection_i)
      )
    )
  }
}

