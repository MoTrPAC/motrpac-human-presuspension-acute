## Author: Tyler Sagendorf
## Date: 2025-05-07
##
## Purpose: Create bubble heatmaps of the CAMERA-PR results for specific terms.
## These heatmaps are part of Figure 4.

library(MotrpacHumanPreSuspensionAnalysis)
library(dplyr)

## Folder path relative to precovid-analyses/
folder <- file.path("figures", "landscape", "CAMERA-PR")

if (!dir.exists(file.path(folder, "plots_curated"))) {
  dir.create(file.path(folder, "plots_curated"))
}


## Proteomics and transcriptomics ----

# File with terms to include in the bubble heatmaps
path <- file.path(folder, "data",
                  "Curated_multi-tissue_10.3.2025.xlsx")

sheets <- readxl::excel_sheets(path)

id_list <- lapply(sheets, function(sheet_i) {
  readxl::read_xlsx(path, sheet = sheet_i) %>%
    select(ome, set_id) %>%
    mutate(set_id = as.numeric(set_id),
           set_id = sprintf("%05d", set_id))
}) %>%
  setNames(sheets) %>%
  bind_rows() %>%
  summarise(.by = ome,
            set_id = list(set_id)) %>%
  tibble::deframe()

for (ome_i in names(id_list)) {
  plot_enrich_heatmap(
    x = CAMERA_RESULTS,
    set_ids = id_list[[ome_i]],
    selected_ome = ome_i,
    filename = file.path(
      folder, "plots_curated",
      paste0("curated_CAMERA-PR_multi-tissue_",
             ome_i, "_2025-10-03.pdf")
    )
  )
}


## Metabolomics
plot_enrich_heatmap(
  x = CAMERA_RESULTS,
  selected_ome = "metab",
  # Top 6 terms from each contrast
  n_top = 6L,
  filename = file.path(
    folder, "plots_curated",
    "curated_CAMERA-PR_heatmap_exercise_with_controls_multi-tissue_metab_top6_100325.pdf"
  )
)


