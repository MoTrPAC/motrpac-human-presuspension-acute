## Author: Tyler Sagendorf
## Date: 2025-05-09
##
## Purpose: Create bubble heatmaps that display the results of nonparametric
## CAMERA-PR and ORA applied to the results of fuzzy c-means clustering.

library(MotrpacHumanPreSuspension)
library(dplyr)

## Folder to save plots, relative to precovid-analyses/
folder <- file.path("figures", "landscape", "fuzzy_c-means")

# Clear all plots
current_device <- dev.cur()

while (current_device != 1L) {
  dev.off()
  current_device <- dev.cur()
}

## FCM CAMERA-PR bubble heatmaps -----------------------------------------------

if (!dir.exists(file.path(folder, "plots_FCM_CAMERA-PR"))) {
  dir.create(file.path(folder, "plots_FCM_CAMERA-PR"))
}

# FCM_CAMERA is the output of run_cluster_cameraPR() from the data package
comb_df <- distinct(FCM_CAMERA, tissue, assay, collection) %>%
  mutate(across(.cols = everything(),
                .fns = as.character))

# Plots for each combination of tissue, assay, and collection
for (i in seq_len(nrow(comb_df))) {
  comb_df_i <- comb_df[i, ]
  tissue_i <- comb_df_i$tissue
  ome_i <- comb_df_i$assay
  collection_i <- comb_df_i$collection

  n_top <- ifelse(ome_i %in% c("transcript-rna-seq", "prot-pr"), 8L, Inf)

  filename_i <- file.path(folder, "plots_FCM_CAMERA-PR",
                          sprintf("FCM_CAMERA-PR_heatmap_%s_%s_%s.pdf",
                                  tissue_i, ome_i, collection_i))

  # 5% FDR for all omes except phospho, which uses 10% FDR
  padj_cutoff <- 0.05 + 0.05 * (ome_i == "prot-ph")

  column_title <- sprintf(
    paste0(paste(rep(" ", 30L), collapse = ""),
           "CAMERA-PR %s %s %s"),
    tissue_i, ome_i, collection_i
  )

  FCM_CAMERA %>%
    filter(collection == collection_i) %>%
    plot_cluster_enrichment(n_top = n_top,
                            padj_cutoff = padj_cutoff,
                            selected_tissues = tissue_i,
                            selected_ome = ome_i,
                            filename = filename_i,
                            column_title = column_title)
}


## FCM CAMERA-PR bubble heatmaps (multi-tissue) --------------------------------

if (!dir.exists(file.path(folder, "plots_FCM_CAMERA-PR"))) {
  dir.create(file.path(folder, "plots_FCM_CAMERA-PR"))
}

# FCM_CAMERA is the output of run_cluster_cameraPR() from the data package
comb_df_multi <- FCM_CAMERA %>%
  distinct(tissue, assay, collection) %>%
  mutate(across(.cols = everything(),
                .fns = as.character)) %>%
  filter(.by = c(assay, collection),
         n() > 1L) %>%
  select(-tissue) %>%
  distinct()

for (i in seq_len(nrow(comb_df_multi))) {
  comb_df_i <- comb_df_multi[i, ]
  ome_i <- comb_df_i$assay
  collection_i <- comb_df_i$collection

  n_top <- ifelse(ome_i %in% c("transcript-rna-seq", "prot-pr"), 5L, Inf)

  filename_i <- file.path(
    folder, "plots_FCM_CAMERA-PR",
    sprintf("FCM_CAMERA-PR_heatmap_multi-tissue_%s_%s.pdf",
            ome_i, collection_i)
  )

  # 5% FDR for all omes except phospho, which uses 10% FDR
  padj_cutoff <- 0.05 + 0.05 * (ome_i == "prot-ph")

  column_title <- paste("CAMERA-PR", ome_i, collection_i)

  FCM_CAMERA %>%
    filter(collection == collection_i) %>%
    plot_cluster_enrichment(n_top = n_top,
                            padj_cutoff = padj_cutoff,
                            selected_ome = ome_i,
                            filename = filename_i,
                            column_title = column_title)
}


## FCM ORA bubble heatmaps -----------------------------------------------------

if (!dir.exists(file.path(folder, "plots_FCM_ORA"))) {
  dir.create(file.path(folder, "plots_FCM_ORA"))
}

# FCM_ORA is the output of run_cluster_ORA() from the data package
comb_df_ORA <- distinct(FCM_ORA, tissue, assay, collection) %>%
  mutate(across(.cols = everything(),
                .fns = as.character))

# Plots for each combination of tissue, assay, and collection
for (i in seq_len(nrow(comb_df_ORA))) {
  comb_df_i <- comb_df_ORA[i, ]
  tissue_i <- comb_df_i$tissue
  ome_i <- comb_df_i$assay
  collection_i <- comb_df_i$collection

  n_top <- ifelse(ome_i %in% c("transcript-rna-seq", "prot-pr"), 8L, Inf)

  filename_i <- file.path(folder, "plots_FCM_ORA",
                          sprintf("FCM_ORA_heatmap_%s_%s_%s.pdf",
                                  tissue_i, ome_i, collection_i))

  # 5% FDR for all omes except phospho, which uses 10% FDR
  padj_cutoff <- 0.05 + 0.05 * (ome_i == "prot-ph")

  column_title <- sprintf(
    paste0(paste(rep(" ", 30L), collapse = ""),
           "ORA %s %s %s"),
    tissue_i, ome_i, collection_i
  )

  FCM_ORA %>%
    filter(collection == collection_i) %>%
    plot_cluster_enrichment(n_top = n_top,
                            padj_cutoff = padj_cutoff,
                            selected_tissues = tissue_i,
                            selected_ome = ome_i,
                            filename = filename_i,
                            column_title = column_title)
}
