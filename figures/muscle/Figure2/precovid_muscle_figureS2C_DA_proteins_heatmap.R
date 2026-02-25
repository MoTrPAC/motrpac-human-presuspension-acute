#This code uses a configuration file to ensure local file structure and local machine settings can be used properly.

#' Steps for using the MoTrPAC config file:
#' 1. Open the configuration template file in the root folder of this repo called `example_config.json` available in `precovid-analyses/`
#' 2. Save a copy of this template to some other location in your computer. Many choose to place it in their home directory: Example: `~/motrpac_config.json`
#' 3. Customize the fields for your machine and save changes. Remove any "comments" as JSON does not permit comments
#' 4. This script expects the following elements in the config file:
##' a. "gsutil_user" and "gsutil_path": this is the way your machine calls the gsutil, usually "/usr/bin/gsutil"
##' b. "local_path": this is the absolute file path for the `precovid-analyses` repo. Please do not include the final "/"
#' 5. Add the path to this config file below 
# config <- jsonlite::read_json(path = "//mnt/c/Users/dankatz/Documents/motrpac_config.json")
rm(list=ls())
setwd("C:/Users/gayat/Documents/MoTrPAC/precovid_analyses")
config <- rjson::fromJSON(file = 'C:/Users/gayat/Documents/MoTrPAC/motrpac_precovid_config_GI.json')

session_packages = c(
  "circlize",
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "ggpubr",
  "hrbrthemes",
  "magrittr",
  "MotrpacHumanPreSuspension",
  "MotrpacHumanPreSuspensionAnalysis",
  "purrr",
  "RColorBrewer",
  "stringr",
  "tidyverse"
)
for (lib_name in session_packages){
  tryCatch({library(lib_name,character.only = TRUE)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("union","base")
conflict_prefer("intersect","base")
conflict_prefer("setdiff", "base")
conflict_prefer("HUMAN_FEATURE_TO_GENE","MotrpacHumanPreSuspensionAnalysis")
conflict_prefer("CONTRAST_CONVERTER","MotrpacHumanPreSuspensionAnalysis")

repo_local_dir <- config$precovid_repo_path
gsutil_cmd <- config$gsutil_user

###Function to calculate heatmap size
calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}

cell_height <- unit(6, "mm")
cell_width <- unit(8, "mm")

################################################################################
sel_tissue <- "muscle"
sel_ome <- "prot-pr"
sel_epigen <- FALSE
contrasts <- c("EE post 15 min",
               "EE post 3.5 hr",
               "EE post 24 hr",
               "RE post 15 min",
               "RE post 3.5 hr",
               "RE post 24 hr")

prot_da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = sel_ome,
                                                                      selected_tissues = sel_tissue,
                                                                      epigen = sel_epigen,
                                                                      combine_with_featgene = T,
                                                                      single_matrix = T) %>%
  filter(contrast_type=="exercise_with_controls") %>%
  mutate(contrast_short=gsub(" - .*","",contrast_short))

prot_sig <- prot_da %>%
  filter(adj_p_value < 0.05)

prot_sig2 <- prot_sig %>%
  mutate(contrast_short=case_when(contrast_short=="Endur.post_15_30_45_min"~"EE post 15 min",
                                  contrast_short=="Endur.post_3.5_4_hr"~"EE post 3.5 hr",
                                  contrast_short=="Endur.post_24_hr"~"EE post 24 hr",
                                  contrast_short=="Resist.post_15_30_45_min"~"RE post 15 min",
                                  contrast_short=="Resist.post_3.5_4_hr"~"RE post 3.5 hr",
                                  contrast_short=="Resist.post_24_hr"~"RE post 24 hr")) %>%
  mutate(ome=case_when(assay=="prot-ol"~"Proteomics OL",
                       assay=="prot-pr"~"Proteomics PR",
                       assay=="prot-ph"~"Phosphoproteomics",
                       assay=="transcript-rna-seq"~"Transcriptomics",
                       assay=="epigen-atac-seq"~"ATAC-seq",
                       assay=="epigen-methylcap-seq"~"MethylCap-seq",
                       .default = "Metabolomics")) %>%
  mutate(ex_modality=gsub("_.*","",contrast_short))


features_to_plot <- prot_sig2 %>%
  select(feature_id) %>%
  distinct() %>%
  pull()

sel_prot_da <- prot_da %>%
  filter(feature_id %in% features_to_plot) %>%
  separate(contrast_short, into = c("Modality","Timepoint"),
           remove = FALSE, extra = "merge") %>%
  mutate(
    Timepoint = case_when(
      Timepoint == "post_15_30_45_min" ~ "post 15 min",
      Timepoint == "post_3.5_4_hr" ~ "post 3.5 hr", 
      Timepoint == "post_24_hr" ~ "post 24 hr",
      TRUE ~ Timepoint
    ),
    Timepoint = gsub("_", " ", Timepoint),
    feature_id2 = paste0(gene_symbol, " (", feature_id, ")")
  ) %>%
  group_by(feature_id) %>%
  filter(any(adj_p_value < 0.05)) %>%
  ungroup()

sel_prot_da2 <- sel_prot_da %>%
  select(feature_id2, contrast_short, tissue, z.std) %>%
  pivot_wider(names_from = c(tissue, contrast_short), values_from = z.std) %>%
  arrange(feature_id2) %>%
  column_to_rownames(var="feature_id2") %>%
  replace(is.na(.), 0) %>%
  select(muscle_Endur.post_15_30_45_min,
         muscle_Endur.post_3.5_4_hr,
         muscle_Endur.post_24_hr,
         muscle_Resist.post_15_30_45_min,
         muscle_Resist.post_3.5_4_hr,
         muscle_Resist.post_24_hr)

adj_p_mat <- sel_prot_da %>%
  select(feature_id2, contrast_short, tissue, adj_p_value) %>%
  pivot_wider(names_from = c(tissue, contrast_short), values_from = adj_p_value) %>%
  arrange(feature_id2) %>%
  column_to_rownames(var="feature_id2") %>%
  replace(is.na(.), 1) %>%
  select(muscle_Endur.post_15_30_45_min,
         muscle_Endur.post_3.5_4_hr,
         muscle_Endur.post_24_hr,
         muscle_Resist.post_15_30_45_min,
         muscle_Resist.post_3.5_4_hr,
         muscle_Resist.post_24_hr)

col.anno <- sel_prot_da %>%
  distinct(tissue, contrast_short, .keep_all = TRUE) %>%
  select(tissue, contrast_short, Modality, Timepoint) %>%
  mutate(
    Modality = case_when(Modality=="Endur" ~ "EE", Modality=="Resist" ~ "RE")
  ) %>% 
  mutate(Timepoint = factor(Timepoint, levels = c("post 15 min", "post 3.5 hr", "post 24 hr"))) %>%
  group_by(tissue, Modality, Timepoint) %>%
  arrange(Modality, Timepoint, .by_group = FALSE) %>%
  mutate(tissue_contrast = paste0(tissue, "_", contrast_short)) %>%
  column_to_rownames(var="tissue_contrast")

column_ha <- columnAnnotation(
  Modality = col.anno$Modality,
  Timepoint = col.anno$Timepoint,
  col = list(
    Timepoint = c("post 15 min" = "#AE76A3", 
                  "post 3.5 hr" = "#882E72", 
                  "post 24 hr" = "#61194F"),
    Modality = c("EE" = "#d95f02", 
                 "RE" = "#1b9e77")
  ),
  show_annotation_name = FALSE
)


ht1 <- ComplexHeatmap::Heatmap(
  as.matrix(sel_prot_da2),
  height = cell_height * nrow(sel_prot_da2),
  width = cell_width * ncol(sel_prot_da2),
  # col = colorRamp2(c(floor(min(sel_prot_da2)),0,ceiling(max(sel_prot_da2))), 
  #                  brewer.pal(11, "RdBu")[c(10, 6, 2)]),
  top_annotation = column_ha,
  column_title = NULL,
  #column_title_gp = gpar(fontsize = 16, fontface = "italic", col = "grey20"),
  border = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  column_split = factor(
    rep(c("EE", "RE"), each = 3),
    levels = c("EE", "RE")),
  column_gap = unit(2, "mm"),
  row_names_gp = gpar(col = "grey20", fontsize=12),
  column_names_gp = gpar(fontsize=10),
  rect_gp = gpar(col = "#2F4F4F", lwd = 1),
  show_column_names = FALSE,
  name = "Z-score",
  heatmap_legend_param = list(direction = "horizontal"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (adj_p_mat[i, j] < 0.05) {
      grid.points(
        x, y,
        pch = 16,                
        size = unit(6, "pt"),    
        gp = gpar(col = "black")
      )
    }
  }
)

size <- calc_ht_size(ht1)

# jpeg("precovid_muscle_FigureS2C_DA_proteins_heatmap.jpeg",
#      width = size[1], height = size[2]+0.5, units = "in", res = 1200)
pdf("precovid_muscle_figureS2C_DA_proteins_heatmap.pdf",
     width = size[1], height = size[2]+0.5)
draw(ht1,
     heatmap_legend_side = "bottom")
dev.off()
