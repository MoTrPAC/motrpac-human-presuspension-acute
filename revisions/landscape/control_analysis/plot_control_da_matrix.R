# modified from figure 2 - da feature matrix
plot_control_da_matrix = function(control_only_analysis) {

  timepoint_levels = c("pre_exercise",
                       "during_20_min",
                       "during_40_min",
                       "post_10_min",
                       "post_15_30_45_min",
                       "post_3.5_4_hr",
                       "post_24_hr")

  timepoint_short = c(
    "pre_exercise" = "Pre",
    "during_20_min" = "D20",
    "during_40_min" = "D40",
    "post_10_min" = "P10",
    "post_15_30_45_min" = "P1545",
    "post_3.5_4_hr" = "P35",
    "post_24_hr" = "P24"
  )

  con_slice = control_only_analysis %>%
    dplyr::mutate(
      Ome = dplyr::case_when(
        assay %in% c("prot-ol", "prot-pr") ~ "Proteomics",
        assay == "prot-ph" ~ "Phosphoproteomics",
        grepl("metab", assay) ~ "Metabolomics",
        grepl("rna", assay) ~ "Transcriptomics",
        grepl("atac", assay) ~ "Chromatin Accessibility (ATAC)",
        grepl("methyl", assay) ~ "Methylation",
        TRUE ~ assay
      ),
      Ome = factor(Ome, levels = c("Chromatin Accessibility (ATAC)",
                                   "Transcriptomics",
                                   "Proteomics",
                                   "Phosphoproteomics",
                                   "Metabolomics",
                                   "Methylation")),
      Timepoint = dplyr::recode(as.character(Timepoint), !!!timepoint_short),
      Timepoint = factor(Timepoint, levels = unname(timepoint_short))
    ) %>%
    dplyr::group_by(Ome, tissue, Timepoint) %>%
    dplyr::summarize(
      num_feat = dplyr::n(),
      percent_sig_fdr = 100 * sum(adj_p_value < 0.05) / dplyr::n(),
      num_sig_fdr005 = sum(adj_p_value < 0.05),
      .groups = "drop"
    )

  tissues_sorted = sort(unique(con_slice$tissue))

  col_order = expand.grid(
    Timepoint = unname(timepoint_short),
    tissue = tissues_sorted,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(col = paste(Timepoint, tissue, sep = "_")) %>%
    dplyr::pull(col)

  con_perc = con_slice %>%
    dplyr::select(Ome, tissue, Timepoint, percent_sig_fdr) %>%
    tidyr::pivot_wider(values_from = percent_sig_fdr,
                       names_from = c("Timepoint", "tissue")) %>%
    dplyr::arrange(Ome) %>%
    dplyr::select(Ome, dplyr::any_of(col_order)) %>%
    tibble::column_to_rownames(var = "Ome") %>%
    as.matrix()

  con_num = con_slice %>%
    dplyr::select(Ome, tissue, Timepoint, num_sig_fdr005) %>%
    tidyr::pivot_wider(values_from = num_sig_fdr005,
                       names_from = c("Timepoint", "tissue")) %>%
    dplyr::arrange(Ome) %>%
    dplyr::select(Ome, dplyr::any_of(col_order)) %>%
    tibble::column_to_rownames(var = "Ome") %>%
    as.matrix()

  col_tissues = stringr::str_split_fixed(colnames(con_perc), "_", 2)[, 2]
  gaps_col = which(diff(match(col_tissues, tissues_sorted)) != 0)

  tp_colors = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  tp_colors = tp_colors[names(tp_colors) %in% names(timepoint_short)]
  names(tp_colors) = timepoint_short[names(tp_colors)]

  ann_colors = list(
    Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS,
    Ome = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS,
    Timepoint = tp_colors
  )

  annotation_data = rownames(con_perc) %>%
    as.data.frame() %>%
    dplyr::rename(Ome = ".") %>%
    dplyr::mutate(Ome = factor(Ome, levels = levels(con_slice$Ome)))

  annotation_col = colnames(con_perc) %>%
    as.data.frame() %>%
    dplyr::rename(col = ".") %>%
    tidyr::separate(col, into = c("Timepoint", "Tissue"), sep = "_") %>%
    dplyr::mutate(Timepoint = factor(Timepoint, levels = unname(timepoint_short)))
  rownames(annotation_col) = colnames(con_perc)

  ComplexHeatmap::pheatmap(
    con_perc,
    heatmap_legend_param = list(title = "% features DA\n(Control)"),
    border_color = "gray3",
    scale = "none",
    color = c("white", "#9e9ac8"),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_row = annotation_data,
    annotation_colors = ann_colors,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_names_col = FALSE,
    display_numbers = as.matrix(con_num),
    number_color = "black",
    fontsize_number = 15,
    angle_col = "90",
    gaps_col = gaps_col,
    na_col = "white"
  )
}
