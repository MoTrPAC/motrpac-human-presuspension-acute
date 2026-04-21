plot_precawg_mvsf_matrix = function(precawg_and_sex_diff) {

  timepoint_short = c(
    "pre_exercise" = "Pre",
    "during_20_min" = "D20",
    "during_40_min" = "D40",
    "post_10_min" = "P10",
    "post_15_30_45_min" = "P1545",
    "post_3.5_4_hr" = "P35",
    "post_24_hr" = "P24"
  )

  assay_levels = c("transcript-rna-seq", "prot-ol", "prot-pr", "prot-ph", "metab")

  features = precawg_and_sex_diff %>%
    dplyr::mutate(
      Group = dplyr::recode(randomGroupCode,
                            "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON"),
      assay = factor(assay, levels = assay_levels),
      tp_short = dplyr::recode(as.character(Timepoint), !!!timepoint_short),
      tp_short = factor(tp_short, levels = unname(timepoint_short))
    )

  summary_df = features %>%
    dplyr::group_by(assay, tissue, Group, tp_short) %>%
    dplyr::summarize(
      num_feat = dplyr::n(),
      num_da = sum(adj_p_value_precawg < 0.05, na.rm = TRUE),
      percent_da = 100 * num_da / num_feat,
      num_also_mvsf_da = sum(adj_p_value_precawg < 0.05 & adj_p_value_sex_diff < 0.05,
                             na.rm = TRUE),
      pct_sex_da = 100 * num_also_mvsf_da / num_da,
      .groups = "drop"
    )

  # Build column order: Group (EE, RE) x Tissue x Timepoint
  col_order = expand.grid(
    tp_short = unname(timepoint_short),
    tissue = sort(unique(summary_df$tissue)),
    Group = c("EE", "RE"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(col = paste(Group, tissue, tp_short, sep = "_")) %>%
    dplyr::pull(col)

  make_mat = function(df, value_col) {
    df %>%
      dplyr::select(assay, tissue, Group, tp_short, value = dplyr::all_of(value_col)) %>%
      tidyr::pivot_wider(values_from = value, names_from = c("Group", "tissue", "tp_short")) %>%
      dplyr::arrange(assay) %>%
      dplyr::select(assay, dplyr::any_of(col_order)) %>%
      tibble::column_to_rownames("assay") %>%
      as.matrix()
  }

  da_feat_perc = make_mat(summary_df, "percent_da")
  da_feat_num = make_mat(summary_df, "num_da")
  da_sex_pct = make_mat(summary_df, "pct_sex_da")

  # Annotation rows (assay bar)
  assay_cols = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS
  assay_cols["metab"] = assay_cols[grep("metab", names(assay_cols), value = TRUE)[1]]
  assay_cols = assay_cols[names(assay_cols) %in% rownames(da_feat_perc)]

  annotation_data = rownames(da_feat_perc) %>%
    as.data.frame() %>%
    dplyr::rename(assay = ".") %>%
    dplyr::mutate(assay = factor(assay, levels = assay_levels))
  rownames(annotation_data) = rownames(da_feat_perc)

  # Annotation cols (Group + Tissue + Timepoint bars)
  tp_colors = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  tp_colors = tp_colors[names(tp_colors) %in% names(timepoint_short)]
  names(tp_colors) = timepoint_short[names(tp_colors)]

  extra_Group_colors = c(EE = "#d95f02", RE = "#1b9e77", CON = "#7570b3")

  annotation_col = colnames(da_feat_perc) %>%
    as.data.frame() %>%
    dplyr::rename(col = ".") %>%
    tidyr::separate(col, into = c("Group", "Tissue", "Timepoint"), sep = "_") %>%
    dplyr::mutate(
      Group = factor(Group, levels = c("EE", "RE", "CON")),
      Timepoint = factor(Timepoint, levels = unname(timepoint_short))
    )
  rownames(annotation_col) = colnames(da_feat_perc)

  ann_colors = list(
    Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS,
    Group = extra_Group_colors,
    Timepoint = tp_colors,
    assay = assay_cols
  )

  n_ee_cols = sum(startsWith(colnames(da_sex_pct), "EE"))

  ht = ComplexHeatmap::pheatmap(
    da_sex_pct,
    heatmap_legend_param = list(title = "% precawg DA\nalso sex DA"),
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
    angle_col = "90",
    gaps_col = n_ee_cols,
    na_col = "white",
    cell_fun = function(j, i, x, y, width, height, fill) {
      n = da_feat_num[i, j]
      pct = da_sex_pct[i, j]
      if (!is.na(n)) {
        grid::grid.text(as.character(n), x, y + height * 0.18,
                        gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold"))
      }
      pct_label = if (is.nan(pct)) "0%" else if (is.na(pct)) "NA%" else sprintf("%.0f%%", pct)
      grid::grid.text(pct_label, x, y - height * 0.18,
                      gp = grid::gpar(fontsize = 10, col = "grey30"))
    }
  )

  out_dir = here::here("revisions/landscape/sex_sensitivity_analysis/")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file = file.path(out_dir, "precawg_mvsf_matrix.png")

  png(out_file, width = 14, height = 5, units = "in", res = 150)
  ComplexHeatmap::draw(ht, merge_legends = TRUE,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
  dev.off()

  ComplexHeatmap::draw(ht, merge_legends = TRUE,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
}
