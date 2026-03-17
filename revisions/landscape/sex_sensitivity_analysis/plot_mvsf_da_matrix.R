plot_mvsf_da_matrix = function(direct_m_vs_f) {

  timepoint_short = c(
    "pre_exercise" = "Pre",
    "during_20_min" = "D20",
    "during_40_min" = "D40",
    "post_10_min" = "P10",
    "post_15_30_45_min" = "P1545",
    "post_3.5_4_hr" = "P35",
    "post_24_hr" = "P24"
  )

  mvsf_annotated = direct_m_vs_f %>%
    dplyr::mutate(
      contrast_short = gsub("sex_group_timepoint|[()|]", "", contrast),
      contrast_left = stringr::str_split_fixed(contrast_short, " - ", 2)[,1],
      randomGroupCode = stringr::str_split_fixed(contrast_left, "\\.", 3)[,2],
      Timepoint = sub("^[^.]+\\.[^.]+\\.", "", contrast_left),
      Timepoint = dplyr::recode(Timepoint, !!!timepoint_short),
      Timepoint = factor(Timepoint, levels = unname(timepoint_short))
    ) %>%
    dplyr::select(-contrast_left, -contrast_short)

  mvsf_slice = mvsf_annotated %>%
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
      Group = dplyr::recode(randomGroupCode,
                            "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON"),
      Group = factor(Group, levels = c("EE", "RE", "CON"))
    ) %>%
    dplyr::group_by(Ome, tissue, Group, Timepoint) %>%
    dplyr::summarize(
      num_feat = dplyr::n(),
      percent_sig_fdr = 100 * sum(adj_p_value < 0.05) / dplyr::n(),
      num_sig_fdr005 = sum(adj_p_value < 0.05),
      .groups = "drop"
    )

  # build explicit column order: Group (EE then RE) > Timepoint (factor order) > tissue
  col_order = expand.grid(
    tissue = sort(unique(mvsf_slice$tissue)),
    Timepoint = unname(timepoint_short),
    Group = c("EE", "RE"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(col = paste(Group, Timepoint, tissue, sep = "_")) %>%
    dplyr::pull(col)

  mvsf_perc = mvsf_slice %>%
    dplyr::select(Ome, tissue, Timepoint, percent_sig_fdr, Group) %>%
    tidyr::pivot_wider(values_from = percent_sig_fdr,
                       names_from = c("Group", "Timepoint", "tissue")) %>%
    dplyr::arrange(Ome) %>%
    dplyr::select(Ome, dplyr::any_of(col_order)) %>%
    tibble::column_to_rownames(var = "Ome") %>%
    as.matrix()

  mvsf_num = mvsf_slice %>%
    dplyr::select(Ome, tissue, Timepoint, num_sig_fdr005, Group) %>%
    tidyr::pivot_wider(values_from = num_sig_fdr005,
                       names_from = c("Group", "Timepoint", "tissue")) %>%
    dplyr::arrange(Ome) %>%
    dplyr::select(Ome, dplyr::any_of(col_order)) %>%
    tibble::column_to_rownames(var = "Ome") %>%
    as.matrix()

  tp_colors = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  tp_colors = tp_colors[names(tp_colors) %in% names(timepoint_short)]
  names(tp_colors) = timepoint_short[names(tp_colors)]

  extra_Group_colors = c(EE = "#d95f02", RE = "#1b9e77", CON = "#7570b3")
  ann_colors = list(
    Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS,
    Group = extra_Group_colors,
    Ome = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS,
    Timepoint = tp_colors
  )

  annotation_data = rownames(mvsf_perc) %>%
    as.data.frame() %>%
    dplyr::rename(Ome = ".") %>%
    dplyr::mutate(Ome = factor(Ome, levels = levels(mvsf_slice$Ome)))

  # column names are Group_Timepoint_tissue, all parts are underscore-free
  annotation_col = colnames(mvsf_perc) %>%
    as.data.frame() %>%
    dplyr::rename(col = ".") %>%
    tidyr::separate(col, into = c("Group", "Timepoint", "Tissue"), sep = "_") %>%
    dplyr::mutate(Timepoint = factor(Timepoint, levels = unname(timepoint_short)))
  rownames(annotation_col) = colnames(mvsf_perc)

  n_ee_cols = sum(startsWith(colnames(mvsf_perc), "EE"))

  ComplexHeatmap::pheatmap(
    mvsf_perc,
    heatmap_legend_param = list(title = "% features DA\n(Male vs Female)"),
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
    display_numbers = as.matrix(mvsf_num),
    number_color = "black",
    fontsize_number = 15,
    angle_col = "90",
    gaps_col = n_ee_cols,
    na_col = "white"
  )
}
