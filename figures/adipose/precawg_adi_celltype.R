library(MotrpacHumanPreSuspension)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(TMSig)

celltype_camera <- run_cameraPR(selected_omes = "transcript-rna-seq",
                                selected_tissues = "adipose",
                                path_to_gmt = "path/cell_marker_genes.gmt",
                                overlap_cutoff = 0) 

# Define the contrast order (only Control comparisons)
contrast_order <- c(
  "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
)

# Filter to keep only selected contrasts
celltype_camera_filtered <- celltype_camera %>%
  filter(contrast %in% contrast_order) %>%
  mutate(contrast = factor(contrast, levels = contrast_order))  # Ensure correct order

celltype_camera_filtered <- celltype_camera_filtered %>%
  mutate(
    Exercise_Type = case_when(
      grepl("ADUEndur", contrast) ~ "EE",
      grepl("ADUResist", contrast) ~ "RE"
    ),
    Timepoint = case_when(
      grepl("15_30_45", contrast) ~ "Post 45 min",
      grepl("3.5_4_hr", contrast) ~ "Post 4 hr",
      grepl("24_hr", contrast) ~ "Post 24 hr"
    )
  )
contrast_colors <- c("EE" = "#d95f02", "RE" = "#1b9e77")
time_colors <- c(
  "Post 45 min" = "#AE76A3",
  "Post 4 hr"   = "#882E72",
  "Post 24 hr"  = "#61194F"
)
# Build annotation vectors from contrast_order
Exercise_vec <- ifelse(
  grepl("ADUEndur", contrast_order), "EE", "RE"
)

Time_vec <- case_when(
  grepl("15_30_45", contrast_order) ~ "Post 45 min",
  grepl("3.5_4_hr",  contrast_order) ~ "Post 4 hr",
  grepl("24_hr",     contrast_order) ~ "Post 24 hr"
)

# Correct annotation: length = number of CONTRASTS, not number of rows
col_anno <- HeatmapAnnotation(
  Exercise = Exercise_vec,
  Time     = Time_vec,
  col = list(
    Exercise = contrast_colors,
    Time = time_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Exercise Type"),
    Time = list(title = "Timepoint")
  )
)
group_labels <- ifelse(
  grepl("^group_timepointADUEndur", contrast_order), "EE", "RE"
)

# Figure 4A
celltype_camera_filtered %>% 
  enrichmap(
    n_top = Inf,
    plot_sig_only = FALSE,
    set_column = "set",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    heatmap_args = list(
      cluster_rows = FALSE,
      top_annotation = col_anno,
      column_order = contrast_order,
      column_split = group_labels,
      show_column_names = FALSE
    )
  )

