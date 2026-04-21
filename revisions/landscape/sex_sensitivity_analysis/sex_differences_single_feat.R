sex_differences_single_feature = function(feature,
                                          selected_tissues = "all",
                                          output_file = NULL,
                                          scale_factor = 1,
                                          include_legend = TRUE,
                                          legend_position = "right",
                                          verbose = TRUE) {

  if(all(selected_tissues == "all"))
    selected_tissues = MotrpacHumanPreSuspensionAnalysis::tissue_available_list()

  # Match feature to gene symbol and assay
  feature_info = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(tolower(feature_id) == tolower(feature) |
                    tolower(gene_symbol) == tolower(feature) |
                    tolower(refmet_name) == tolower(feature)) %>%
    dplyr::semi_join(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE, by = "gene_symbol") %>%
    dplyr::mutate(feature_id = dplyr::case_when(
      !is.na(refmet_name) ~ refmet_name,
      TRUE ~ feature_id
    ))

  gene_symbol_options <- as.character(unique(feature_info$gene_symbol))
  refmet_options <- as.character(unique(feature_info$refmet_name))
  label_options <- union(gene_symbol_options, refmet_options)
  label_options <- label_options[!is.na(label_options)]
  feature_label <- label_options[1]

  if(length(label_options) > 1 & verbose) {
    message("if there are multiple gene symbols or refmet names corresponding
            to the input, the first is chosen for the purpose of labeling, use caution")
  }

  # Load QC normalized data and extract feature values
  qc_data = MotrpacHumanPreSuspensionData::load_qc(
    selected_tissues = selected_tissues
  )

  pheno_data = MotrpacHumanPreSuspensionData::pheno$data %>%
    dplyr::select(vialLabel, Sex, randomGroupCode, Timepoint)

  individual_data = lapply(names(qc_data), function(tissue_name) {
    lapply(names(qc_data[[tissue_name]]), function(ome_name) {
      mat = qc_data[[tissue_name]][[ome_name]][["qc_norm"]]
      if(is.null(mat)) return(NULL)

      feature_rows = which(rownames(mat) %in% feature_info$feature_id)
      if(length(feature_rows) == 0) return(NULL)

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

  if(nrow(individual_data) == 0) {
    stop("No data corresponds to your requested feature.
         Please double check your input matches something in the feature to gene
         mapping's feature_id column, or, for metabolites, the refmet_name column")
  }

  # Compute sex-stratified summary stats
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

  # Compute overall summary stats
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
    "pre_exercise" = "Pre",
    "during_20_min" = "D20M",
    "during_40_min" = "D40M",
    "post_10_min" = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr" = "P3.5/4H",
    "post_24_hr" = "P24H"
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

  series_colors = c("Female" = "#E87461", "Male" = "#5B9BD5", "Overall" = "black")
  label_map = c("ADUControl" = "CON", "ADUEndur" = "EE", "ADUResist" = "RE")

  sc = scale_factor * 0.7

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
    facet_grid(tissue_assay ~ randomGroupCode,
               scales = "free_y",
               labeller = labeller(randomGroupCode = label_map)) +
    ggtitle(feature_label) +
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

  if(include_legend == TRUE & !is.null(legend_position)) {
    g <- g + theme(legend.position = legend_position)
  }
  #for output file I manually make these for different tissue/combinations for dimensions
  if(!is.null(output_file)) {
    if(include_legend == TRUE) {
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

  gene_symbol_options <- as.character(unique(feature_info$gene_symbol))
  refmet_options <- as.character(unique(feature_info$refmet_name))
  label_options <- union(gene_symbol_options, refmet_options)
  label_options <- label_options[!is.na(label_options)]
  feature_label <- label_options[1]

  if(length(label_options) > 1 & verbose) {
    message("if there are multiple gene symbols or refmet names corresponding
            to the input, the first is chosen for the purpose of labeling, use caution")
  }

  data_filtered = differential_analysis_results %>%
    dplyr::filter(feature_id %in% feature_info$feature_id) %>%
    dplyr::filter(if(all(selected_tissues == "all")) TRUE else tissue %in% selected_tissues)

  if(nrow(data_filtered) == 0) {
    stop("No data corresponds to your requested feature.")
  }

  timepoint_levels = c("Pre", "D20M", "D40M", "P10M", "P15-45M", "P3.5/4H", "P24H")
  timepoint_recode = c(
    "pre_exercise" = "Pre",
    "during_20_min" = "D20M",
    "during_40_min" = "D40M",
    "post_10_min" = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr" = "P3.5/4H",
    "post_24_hr" = "P24H"
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

  series_colors = c("Female" = "#E87461", "Male" = "#5B9BD5")
  label_map = c("ADUControl" = "CON", "ADUEndur" = "EE", "ADUResist" = "RE")

  sc = scale_factor * 0.7

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
    facet_grid(tissue_assay ~ randomGroupCode,
               scales = "free_y",
               labeller = labeller(randomGroupCode = label_map)) +
    ggtitle(feature_label) +
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

  if(include_legend == TRUE & !is.null(legend_position)) {
    g <- g + theme(legend.position = legend_position)
  }

  if(!is.null(output_file)) {
    if(include_legend == TRUE) {
      ggsave(g, filename = output_file,
             height = 2.2 * sc, width = 2.45 * sc, dpi = 600, units = "in")
    } else {
      ggsave(g, filename = output_file,
             height = 2.2 * sc, width = 1.75 * sc, dpi = 600, units = "in")
    }
  }

  return(g)

}
