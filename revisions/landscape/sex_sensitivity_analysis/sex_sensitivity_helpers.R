# Helper functions for sex_sensitivity_analysis.Rmd
# All output figures and files should be written to the outputs/ subdirectory.

# Pick a readable (black/white) text color for a given background hex, by luminance.
.text_on = function(hex) {
  rgb = grDevices::col2rgb(hex) / 255
  lum = 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  if (lum > 0.55) "black" else "white"
}

# Build coloured facet strips for the single-feature panels: exercise-group columns
# use the canonical group palette, tissue rows use the canonical tissue palette (both
# shared with the other figure scripts). When every panel shares a single ome, that ome
# is dropped from the right-side row strips and returned in the title instead, since it
# is redundant when identical across all rows.
.sex_feature_strip_spec = function(plot_data, label_map, sc) {

  group_pal = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS
  names(group_pal) = dplyr::recode(names(group_pal),
                                   "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON")
  tissue_pal = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS

  resolve_grp = function(v) {
    v = dplyr::recode(as.character(v), !!!label_map)
    if (v %in% names(group_pal)) unname(group_pal[v]) else "grey85"
  }
  resolve_tissue = function(v) {
    key = tolower(as.character(v))
    if (key %in% names(tissue_pal)) unname(tissue_pal[key]) else "grey85"
  }

  # The single-feature panels key the ome on `assay`; the logFC panels use `platform`.
  ome_col = if ("assay" %in% names(plot_data)) plot_data$assay else plot_data$platform
  omes_present = unique(ifelse(grepl("metab", ome_col), "metab", as.character(ome_col)))
  single_ome = length(omes_present) == 1
  row_var = if (single_ome) "tissue" else "tissue_assay"

  col_levels = levels(droplevels(as.factor(plot_data$randomGroupCode)))
  row_levels = levels(droplevels(as.factor(plot_data[[row_var]])))
  col_cols = vapply(col_levels, resolve_grp, character(1))
  row_cols = if (single_ome) vapply(row_levels, resolve_tissue, character(1))
             else rep("grey85", length(row_levels))

  strip = ggh4x::strip_themed(
    background_x = lapply(col_cols, function(cc) ggplot2::element_rect(fill = cc, colour = "grey20")),
    background_y = lapply(row_cols, function(cc) ggplot2::element_rect(fill = cc, colour = "grey20")),
    text_x = lapply(col_cols, function(cc) ggplot2::element_text(colour = .text_on(cc), face = "bold", size = 8 * sc)),
    text_y = lapply(row_cols, function(cc) ggplot2::element_text(colour = .text_on(cc), face = "bold", size = 8 * sc))
  )

  list(
    strip = strip,
    facet_formula = stats::as.formula(paste0(row_var, " ~ randomGroupCode")),
    title = function(feature_label) {
      if (single_ome) paste0(feature_label, " — ", omes_present) else feature_label
    }
  )
}

# ── Matrix plots ──────────────────────────────────────────────────────────────
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

  assay_cols = MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS
  assay_cols["metab"] = assay_cols[grep("metab", names(assay_cols), value = TRUE)[1]]
  assay_cols = assay_cols[names(assay_cols) %in% rownames(da_feat_perc)]

  annotation_data = rownames(da_feat_perc) %>%
    as.data.frame() %>%
    dplyr::rename(assay = ".") %>%
    dplyr::mutate(assay = factor(assay, levels = assay_levels))
  rownames(annotation_data) = rownames(da_feat_perc)

  tp_colors = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  tp_colors = tp_colors[names(tp_colors) %in% names(timepoint_short)]
  names(tp_colors) = timepoint_short[names(tp_colors)]

  # Match the canonical exercise-group palette used across the other figure scripts,
  # remapped from the ADU* keys to the EE/RE/CON labels used in this annotation.
  extra_Group_colors = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS
  names(extra_Group_colors) = dplyr::recode(names(extra_Group_colors),
                                            "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON")

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

  out_dir = here::here("revisions/landscape/sex_sensitivity_analysis/outputs/")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file = file.path(out_dir, "precawg_mvsf_matrix.pdf")

  pdf(out_file, width = 14, height = 5)
  ComplexHeatmap::draw(ht, merge_legends = TRUE,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
  dev.off()

  ComplexHeatmap::draw(ht, merge_legends = TRUE,
                       heatmap_legend_side = "right",
                       annotation_legend_side = "right")
}


# ── Single-feature plots ──────────────────────────────────────────────────────

sex_differences_single_feature = function(feature,
                                          selected_tissues = "all",
                                          output_file = NULL,
                                          scale_factor = 1,
                                          include_legend = TRUE,
                                          legend_position = "right",
                                          verbose = TRUE) {

  if (all(selected_tissues == "all"))
    selected_tissues = MotrpacHumanPreSuspensionAnalysis::tissue_available_list()

  feature_info = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(tolower(feature_id) == tolower(feature) |
                    tolower(gene_symbol) == tolower(feature) |
                    tolower(refmet_name) == tolower(feature)) %>%
    dplyr::semi_join(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE, by = "gene_symbol") %>%
    dplyr::mutate(feature_id = dplyr::case_when(
      !is.na(refmet_name) ~ refmet_name,
      TRUE ~ feature_id
    ))

  gene_symbol_options = as.character(unique(feature_info$gene_symbol))
  refmet_options = as.character(unique(feature_info$refmet_name))
  label_options = union(gene_symbol_options, refmet_options)
  label_options = label_options[!is.na(label_options)]
  feature_label = label_options[1]

  if (length(label_options) > 1 & verbose) {
    message("if there are multiple gene symbols or refmet names corresponding
            to the input, the first is chosen for the purpose of labeling, use caution")
  }

  qc_data = MotrpacHumanPreSuspensionData::load_qc(
    selected_tissues = selected_tissues
  )

  pheno_data = MotrpacHumanPreSuspensionData::pheno$data %>%
    dplyr::select(vialLabel, Sex, randomGroupCode, Timepoint)

  individual_data = lapply(names(qc_data), function(tissue_name) {
    lapply(names(qc_data[[tissue_name]]), function(ome_name) {
      mat = qc_data[[tissue_name]][[ome_name]][["qc_norm"]]
      if (is.null(mat)) return(NULL)

      feature_rows = which(rownames(mat) %in% feature_info$feature_id)
      if (length(feature_rows) == 0) return(NULL)

      mat[feature_rows, , drop = FALSE] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("feature_id") %>%
        tidyr::pivot_longer(-feature_id, names_to = "vialLabel", values_to = "value") %>%
        dplyr::inner_join(pheno_data, by = "vialLabel") %>%
        dplyr::mutate(tissue = tissue_name, assay = ome_name)
    }) %>%
      dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows()

  if (nrow(individual_data) == 0) {
    stop("No data corresponds to your requested feature.
         Please double check your input matches something in the feature to gene
         mapping's feature_id column, or, for metabolites, the refmet_name column")
  }

  sex_summary = individual_data %>%
    dplyr::group_by(tissue, assay, feature_id, Sex, randomGroupCode, Timepoint) %>%
    dplyr::summarize(
      Mean = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Count = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SE = SD / sqrt(Count),
      CI_95 = qt((1 + 0.95) / 2, Count - 1),
      CI_low = Mean - CI_95 * SE,
      CI_high = Mean + CI_95 * SE,
      series = Sex
    )

  overall_summary = individual_data %>%
    dplyr::group_by(tissue, assay, feature_id, randomGroupCode, Timepoint) %>%
    dplyr::summarize(
      Mean = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      Count = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SE = SD / sqrt(Count),
      CI_95 = qt((1 + 0.95) / 2, Count - 1),
      CI_low = Mean - CI_95 * SE,
      CI_high = Mean + CI_95 * SE,
      series = "Overall"
    )

  timepoint_levels = c("Pre", "D20M", "D40M", "P10M", "P15-45M", "P3.5/4H", "P24H")
  timepoint_recode = c(
    "pre_exercise"      = "Pre",
    "during_20_min"     = "D20M",
    "during_40_min"     = "D40M",
    "post_10_min"       = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr"     = "P3.5/4H",
    "post_24_hr"        = "P24H"
  )

  plot_data = dplyr::bind_rows(sex_summary, overall_summary) %>%
    dplyr::filter(randomGroupCode %in% c("ADUEndur", "ADUResist", "ADUControl")) %>%
    dplyr::mutate(
      tissue = stringr::str_to_sentence(tissue),
      assay = ifelse(grepl("metab", assay), "metab", assay),
      Timepoint = dplyr::recode(as.character(Timepoint), !!!timepoint_recode),
      Timepoint = factor(Timepoint, levels = timepoint_levels),
      tissue_assay = stringr::str_c(tissue, " ", assay),
      series = factor(series, levels = c("Female", "Male", "Overall"))
    )

  series_colors = c(MotrpacHumanPreSuspensionAnalysis::HUMAN_SEX_COLORS[c("Female", "Male")],
                    "Overall" = "black")
  label_map = c("ADUControl" = "CON", "ADUEndur" = "EE", "ADUResist" = "RE")

  sc = scale_factor * 0.7

  strip_spec = .sex_feature_strip_spec(plot_data, label_map, sc)

  g = ggplot(plot_data, aes(x = Timepoint, color = series, group = series)) +
    geom_line(aes(y = Mean), linewidth = 0.5 * sc) +
    geom_point(aes(y = Mean), size = 1.8 * sc) +
    geom_errorbar(
      aes(ymin = CI_low, ymax = CI_high),
      width = 0.3 * sc,
      linewidth = 0.4 * sc,
      alpha = 0.6
    ) +
    scale_color_manual(values = series_colors) +
    ggh4x::facet_grid2(strip_spec$facet_formula,
                       scales = "free_y",
                       labeller = labeller(randomGroupCode = label_map),
                       strip = strip_spec$strip) +
    ggtitle(strip_spec$title(feature_label)) +
    ylab("log2(normalized value)") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 9 * sc, color = "black"),
      axis.text.y = element_text(size = 10 * sc, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 11 * sc),
      plot.title = element_text(size = 12 * sc, margin = margin(0, 0, 0, 0)),
      strip.text.x = element_text(size = 8 * sc, margin = margin(0.05 * sc, 0, 0.05 * sc, 0, "cm")),
      legend.text = element_text(size = 9 * sc, margin = margin(0, 0, 0, 0)),
      legend.spacing.x = unit(0.01 * sc, "in"),
      legend.spacing.y = unit(0.1 * sc, "in"),
      legend.box.spacing = unit(0.01 * sc, "in"),
      legend.key.size = unit(0.1 * sc, "in"),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      axis.ticks = element_line(linewidth = 0.3 * sc),
      panel.grid = element_line(linewidth = 0.3 * sc),
      panel.border = element_rect(linewidth = 0.3 * sc),
      plot.margin = margin(0.05 * sc, 0.05 * sc, 0.05 * sc, 0.05 * sc, "in")
    )

  if (include_legend == TRUE & !is.null(legend_position)) {
    g = g + theme(legend.position = legend_position)
  }

  if (!is.null(output_file)) {
    if (include_legend == TRUE) {
      ggsave(g, filename = output_file,
             height = 2.5 + 1.5 * sc, width = 5 + 2.45 * sc, dpi = 600, units = "in")
    } else {
      ggsave(g, filename = output_file,
             height = 2.5 + 1.5 * sc, width = 5 + 1.75 * sc, dpi = 600, units = "in")
    }
  }

  return(g)
}


sex_differences_logfc_feature = function(feature,
                                         differential_analysis_results = precawg_and_sex_diff,
                                         selected_tissues = "all",
                                         output_file = NULL,
                                         scale_factor = 1,
                                         include_legend = TRUE,
                                         legend_position = "right",
                                         verbose = TRUE) {

  feature_info = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(tolower(feature_id) == tolower(feature) |
                    tolower(gene_symbol) == tolower(feature) |
                    tolower(refmet_name) == tolower(feature)) %>%
    dplyr::semi_join(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE, by = "gene_symbol") %>%
    dplyr::mutate(feature_id = dplyr::case_when(
      !is.na(refmet_name) ~ refmet_name,
      TRUE ~ feature_id
    ))

  gene_symbol_options = as.character(unique(feature_info$gene_symbol))
  refmet_options = as.character(unique(feature_info$refmet_name))
  label_options = union(gene_symbol_options, refmet_options)
  label_options = label_options[!is.na(label_options)]
  feature_label = label_options[1]

  if (length(label_options) > 1 & verbose) {
    message("if there are multiple gene symbols or refmet names corresponding
            to the input, the first is chosen for the purpose of labeling, use caution")
  }

  data_filtered = differential_analysis_results %>%
    dplyr::filter(feature_id %in% feature_info$feature_id) %>%
    dplyr::filter(if (all(selected_tissues == "all")) TRUE else tissue %in% selected_tissues)

  if (nrow(data_filtered) == 0) {
    stop("No data corresponds to your requested feature.")
  }

  timepoint_levels = c("Pre", "D20M", "D40M", "P10M", "P15-45M", "P3.5/4H", "P24H")
  timepoint_recode = c(
    "pre_exercise"      = "Pre",
    "during_20_min"     = "D20M",
    "during_40_min"     = "D40M",
    "post_10_min"       = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr"     = "P3.5/4H",
    "post_24_hr"        = "P24H"
  )

  plot_data = data_filtered %>%
    dplyr::select(tissue, platform, feature_id, randomGroupCode, Timepoint,
                  series = sex_comparison,
                  logFC = logFC_sex_diff,
                  CI_low = CI.L_sex_diff,
                  CI_high = CI.R_sex_diff) %>%
    dplyr::mutate(
      tissue = stringr::str_to_sentence(tissue),
      Timepoint = dplyr::recode(as.character(Timepoint), !!!timepoint_recode),
      Timepoint = factor(Timepoint, levels = timepoint_levels),
      tissue_assay = stringr::str_c(tissue, " ", platform),
      series = factor(series, levels = c("Female", "Male"))
    )

  series_colors = MotrpacHumanPreSuspensionAnalysis::HUMAN_SEX_COLORS[c("Female", "Male")]
  label_map = c("ADUControl" = "CON", "ADUEndur" = "EE", "ADUResist" = "RE")

  sc = scale_factor * 0.7

  strip_spec = .sex_feature_strip_spec(plot_data, label_map, sc)

  g = ggplot(plot_data, aes(x = Timepoint, color = series, group = series)) +
    geom_hline(yintercept = 0, linewidth = 0.3 * sc, linetype = "dashed", color = "grey60") +
    geom_line(aes(y = logFC), linewidth = 0.5 * sc) +
    geom_point(aes(y = logFC), size = 1.8 * sc) +
    geom_errorbar(
      aes(ymin = CI_low, ymax = CI_high),
      width = 0.2 * sc,
      linewidth = 0.4 * sc,
      alpha = 0.6
    ) +
    scale_color_manual(values = series_colors) +
    ggh4x::facet_grid2(strip_spec$facet_formula,
                       scales = "free_y",
                       labeller = labeller(randomGroupCode = label_map),
                       strip = strip_spec$strip) +
    ggtitle(strip_spec$title(feature_label)) +
    ylab("logFC") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5 * sc, color = "black"),
      axis.text.y = element_text(size = 7 * sc, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8 * sc),
      plot.title = element_text(size = 10 * sc, margin = margin(0, 0, 0, 0)),
      strip.text.x = element_text(size = 8 * sc, margin = margin(0.05 * sc, 0, 0.05 * sc, 0, "cm")),
      legend.text = element_text(size = 7 * sc, margin = margin(0, 0, 0, 0)),
      legend.spacing.x = unit(0.01 * sc, "in"),
      legend.spacing.y = unit(0.1 * sc, "in"),
      legend.box.spacing = unit(0.01 * sc, "in"),
      legend.key.size = unit(0.1 * sc, "in"),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      axis.ticks = element_line(linewidth = 0.3 * sc),
      panel.grid = element_line(linewidth = 0.3 * sc),
      panel.border = element_rect(linewidth = 0.3 * sc),
      plot.margin = margin(0.05 * sc, 0.05 * sc, 0.05 * sc, 0.05 * sc, "in")
    )

  if (include_legend == TRUE & !is.null(legend_position)) {
    g = g + theme(legend.position = legend_position)
  }

  if (!is.null(output_file)) {
    if (include_legend == TRUE) {
      ggsave(g, filename = output_file,
             height = 2.2 * sc, width = 2.45 * sc, dpi = 600, units = "in")
    } else {
      ggsave(g, filename = output_file,
             height = 2.2 * sc, width = 1.75 * sc, dpi = 600, units = "in")
    }
  }

  return(g)
}

# ── PCA panel ─────────────────────────────────────────────────────────────────

make_pca_panel = function(tiss, ome, ome_data, tp_levels, tp_colors_pca) {
  meta = ome_data[["sample_metadata"]] %>%
    dplyr::mutate(vialLabel = as.character(vialLabel))
  if (nrow(meta) < 6) return(NULL)

  mat = ome_data[["qc_norm"]]
  if (is.null(mat)) return(NULL)
  mat = mat[, colnames(mat) %in% meta$vialLabel, drop = FALSE]
  if (ncol(mat) < 6) return(NULL)

  mat[!is.finite(as.matrix(mat))] = NA
  mat = mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
  if (nrow(mat) < 10 || ncol(mat) < 5) return(NULL)

  pca_res = prcomp(t(mat), scale = TRUE, center = TRUE)
  pct_var = pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

  pca_df = as.data.frame(pca_res$x[, 1:2]) %>%
    tibble::rownames_to_column("vialLabel") %>%
    dplyr::left_join(
      meta %>% dplyr::select(vialLabel, Sex, randomGroupCode, Timepoint),
      by = "vialLabel"
    ) %>%
    dplyr::mutate(
      Group = dplyr::recode(randomGroupCode,
                            "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON"),
      Timepoint = factor(Timepoint, levels = tp_levels)
    )

  ome_label = dplyr::case_when(
    grepl("metab", ome)              ~ "Metab",
    ome == "prot-ph"                 ~ "Phos",
    ome %in% c("prot-pr", "prot-ol") ~ "Prot",
    grepl("rna", ome)                ~ "RNA",
    grepl("atac", ome)               ~ "ATAC",
    grepl("methyl", ome)             ~ "Methyl",
    TRUE                             ~ ome
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint, shape = Group)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(
      aes(x = PC1, y = PC2, group = Sex, linetype = Sex),
      inherit.aes = FALSE,
      type = "norm",
      linewidth = 0.7,
      color = "black"
    ) +
    scale_color_manual(values = tp_colors_pca, drop = FALSE, name = "Timepoint") +
    scale_shape_manual(values = c(EE = 16, RE = 17, CON = 15), name = "Group") +
    scale_linetype_manual(values = c(Female = "solid", Male = "dashed"), name = "Sex (ellipse)") +
    labs(
      title = ome_label,
      x = paste0("PC1 (", round(pct_var[1], 1), "%)"),
      y = paste0("PC2 (", round(pct_var[2], 1), "%)")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 9, face = "bold")
    )
}


# ── Sex scatter ───────────────────────────────────────────────────────────────

.make_sex_scatter = function(tiss, ome_type, data, cor_df, highlight_ids = NULL, n_label = 10) {
  df = data %>% dplyr::filter(tissue == tiss, Ome == ome_type)
  if (nrow(df) == 0) return(NULL)

  cor_sub = cor_df %>% dplyr::filter(tissue == tiss, Ome == ome_type)

  lim = max(abs(c(df$logFC_Female, df$logFC_Male)), na.rm = TRUE) * 1.05

  df = df %>%
    dplyr::mutate(highlight = !is.null(highlight_ids) & feature_id %in% highlight_ids)

  lm_labels = df %>%
    dplyr::group_by(Timepoint, Group) %>%
    dplyr::summarize(
      slope_val = tryCatch(coef(lm(logFC_Male ~ logFC_Female))[2], error = function(e) NA_real_),
      intercept_val = tryCatch(coef(lm(logFC_Male ~ logFC_Female))[1], error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      eq_label = dplyr::case_when(
        is.na(slope_val) ~ "",
        intercept_val >= 0 ~ paste0("y=", sprintf("%.2f", slope_val), "x+", sprintf("%.2f", intercept_val)),
        TRUE ~ paste0("y=", sprintf("%.2f", slope_val), "x", sprintf("%.2f", intercept_val))
      )
    )

  label_df = NULL
  if (!is.null(highlight_ids) && length(highlight_ids) > 0) {
    gene_map = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
      dplyr::select(feature_id, gene_symbol, refmet_name) %>%
      dplyr::distinct()

    label_df = df %>%
      dplyr::filter(highlight) %>%
      dplyr::left_join(gene_map, by = "feature_id") %>%
      dplyr::mutate(
        label = dplyr::coalesce(gene_symbol, refmet_name, feature_id),
        dist = sqrt(logFC_Female^2 + logFC_Male^2)
      ) %>%
      dplyr::group_by(Timepoint, Group) %>%
      dplyr::slice_max(dist, n = n_label, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  # Color the facet strips with the canonical timepoint (rows) and exercise-group
  # (columns) palettes used across the other figure scripts. Timepoint labels may be
  # raw ("post_3.5_4_hr") or prettified ("P3.5/4H"); both resolve to the same color.
  tp_pal = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  tp_label_to_raw = c(
    "Pre" = "pre_exercise",
    "D20" = "during_20_min", "D20M" = "during_20_min",
    "D40" = "during_40_min", "D40M" = "during_40_min",
    "P10" = "post_10_min", "P10M" = "post_10_min",
    "P1545" = "post_15_30_45_min", "P15-45M" = "post_15_30_45_min",
    "P35" = "post_3.5_4_hr", "P3.5/4H" = "post_3.5_4_hr",
    "P24" = "post_24_hr", "P24H" = "post_24_hr"
  )
  group_pal = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS
  names(group_pal) = dplyr::recode(names(group_pal),
                                   "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON")

  resolve_tp = function(v) {
    v = as.character(v)
    if (v %in% names(tp_pal)) return(unname(tp_pal[v]))
    if (v %in% names(tp_label_to_raw)) return(unname(tp_pal[tp_label_to_raw[v]]))
    "grey85"
  }
  resolve_grp = function(v) {
    v = as.character(v)
    if (v %in% names(group_pal)) return(unname(group_pal[v])) else "grey85"
  }

  row_vals = levels(droplevels(as.factor(df$Timepoint)))
  col_vals = levels(droplevels(as.factor(df$Group)))
  row_cols = vapply(row_vals, resolve_tp, character(1))
  col_cols = vapply(col_vals, resolve_grp, character(1))

  strip_spec = ggh4x::strip_themed(
    background_x = lapply(col_cols, function(cc) ggplot2::element_rect(fill = cc, colour = "grey20")),
    background_y = lapply(row_cols, function(cc) ggplot2::element_rect(fill = cc, colour = "grey20")),
    text_x = lapply(col_cols, function(cc) ggplot2::element_text(colour = .text_on(cc), face = "bold")),
    text_y = lapply(row_cols, function(cc) ggplot2::element_text(colour = .text_on(cc), face = "bold"))
  )

  p = ggplot(df, aes(x = logFC_Female, y = logFC_Male)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "grey30", linewidth = 0.5) +
    geom_errorbar(
      aes(ymin = CI.L_Male, ymax = CI.R_Male),
      width = 0, alpha = 0.3, linewidth = 0.35, color = "grey40"
    ) +
    geom_errorbarh(
      aes(xmin = CI.L_Female, xmax = CI.R_Female),
      height = 0, alpha = 0.3, linewidth = 0.35, color = "grey40"
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.7) +
    geom_point(alpha = 0.8, size = 1.5, color = "grey30") +
    geom_text(
      data = cor_sub,
      aes(label = paste0("r=", Pearsons_R, " (n=", n_sig, ")")),
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.2,
      size = 5, fontface = "italic",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = lm_labels,
      aes(label = eq_label),
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 2.4,
      size = 5, fontface = "italic",
      inherit.aes = FALSE
    ) +
    coord_fixed(ratio = 1, xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggh4x::facet_grid2(Timepoint ~ Group, strip = strip_spec) +
    labs(
      title = paste0(tiss, " — ", ome_type),
      x     = "logFC (Female participants)",
      y     = "logFC (Male participants)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 9, face = "bold"),
      axis.title       = element_text(size = 15),
      axis.text        = element_text(size = 12),
      legend.position  = "bottom"
    )

  if (!is.null(label_df) && nrow(label_df) > 0) {
    p = p + ggrepel::geom_text_repel(
      data = label_df,
      aes(x = logFC_Female, y = logFC_Male, label = label),
      size = 3.5,
      fontface = "bold",
      bg.color = "white",
      bg.r = 0.15,
      max.overlaps = 30,
      segment.size = 0.3,
      segment.color = "grey40",
      inherit.aes = FALSE
    )
  }

  p
}
