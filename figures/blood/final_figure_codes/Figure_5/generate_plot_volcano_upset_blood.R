generate_volcano_upset_plot_blood <- function(
    repo_local_dir,
    pdf_output_path,
    omes = "proteomics",
    selected_tissues = "blood",
    GroupCode = "Endur",   # can be "Endur" or "Resist"
    gsutil_cmd = "gsutil",
    nlabels_volcano = 10,
    minsize_upset = 1
) {
  
  valid_omes <- c("proteomics", "metabolomics", "transcriptomics")
  valid_groups <- c("Endur", "Resist", "Control")
  
  if (!(omes %in% valid_omes)) {
    message(glue::glue("⚠️ Skipping: '{omes}' is not a valid omics type. Must be one of {toString(valid_omes)}."))
    return(invisible(NULL))
  }
  
  if (!(GroupCode %in% valid_groups)) {
    message(glue::glue("⚠️ Skipping: '{GroupCode}' is not a valid GroupCode. Must be one of {toString(valid_groups)}."))
    return(invisible(NULL))
  }
  
  # ----------------------------- #
  #   1️⃣ Required Packages
  # ----------------------------- #
  required_libs <- c(
    "dplyr", "ggplot2", "ggrepel", "stringr", "scales",
    "tidyverse", "ComplexUpset",
    "MotrpacHumanPreSuspension", "MotrpacHumanPreSuspensionData"
  )
  
  missing_libs <- required_libs[!sapply(required_libs, requireNamespace, quietly = TRUE)]
  if (length(missing_libs) > 0) {
    stop(paste0("Missing packages: ", paste(missing_libs, collapse = ", "),
                ". Please install them first."))
  }
  suppressPackageStartupMessages(lapply(required_libs, library, character.only = TRUE))
  
  color_exercise <- MotrpacHumanPreSuspension::HUMAN_EXERCISE_GROUP_COLORS
  
  bars_color <- switch(GroupCode,
                       "Resist" = color_exercise[["ADUResist"]],
                       "Endur" = color_exercise[["ADUEndur"]],
                       "Control" = color_exercise[["ADUControl"]],
                       "#000000" # fallback black
  )
  
  # ----------------------------- #
  #   2️⃣ Load differential data
  # ----------------------------- #
  differential_analysis <- load_differential_analysis(
    repo_local_dir = repo_local_dir,
    selected_omes = "all",
    selected_tissues = selected_tissues,
    single_matrix = TRUE,
    combine_with_featgene = TRUE,
    epigen = FALSE
  ) %>%
    dplyr::filter(contrast_type == "exercise_with_controls") %>%
    dplyr::mutate(
      assay = dplyr::case_when(
        stringr::str_detect(assay, "metab") ~ "metab",
        TRUE ~ assay
      ),
      gene_symbol = dplyr::if_else(
        assay == "metab",
        feature_id,     # use feature_id for metabolites
        gene_symbol     # keep original for others
      )
    )
  
  # Map human feature info
  differential_analysis <- differential_analysis %>%
    dplyr::mutate(
      gene = dplyr::case_when(is.na(entrez_gene) ~ feature_id, TRUE ~ entrez_gene),
      assay = dplyr::case_when(
        assay == "prot-ol" ~ "proteomics",
        assay == "metab" ~ "metabolomics",
        assay == "transcript-rna-seq" ~ "transcriptomics",
        TRUE ~ assay
      ),
      randomGroupCode = dplyr::case_when(
        randomGroupCode == "ADUResist" ~ "RE",
        randomGroupCode == "ADUEndur" ~ "EE",
        TRUE ~ randomGroupCode
      )
    )
  
  # ----------------------------- #
  #   3️⃣ Part A: Volcano Plot
  # ----------------------------- #
  selected_group <- ifelse(GroupCode == "Endur", "EE", 
                           ifelse(GroupCode == "Resist", "RE", "Control"))
  
  prot_data <- differential_analysis %>%
    dplyr::filter(assay == omes)
  
  
  ymax <- ceiling(max(-log10(prot_data$p_value), na.rm = TRUE)) + 1
  
  prot_data <- differential_analysis %>%
    dplyr::filter(assay == omes,
                  stringr::str_detect(randomGroupCode, selected_group))
  
  timepoint_levels <- c(
    'during_20_min',
    'during_40_min',
    'post_10_min',
    'post_15_30_45_min',
    'post_3.5_4_hr',
    'post_24_hr'
  )
  
  custom_labels <- c(
    'during_20_min' = "During 20 Min",
    'during_40_min' = "During 40 Min",
    'post_10_min' = "Post 10 Min",
    'post_15_30_45_min' = "Post 30 Min",
    'post_3.5_4_hr' = "Post 3.5 Hour",
    'post_24_hr' = "Post 24 Hour"
  )
  
  prot_data <- prot_data %>% 
    mutate(term = factor(Timepoint, levels = timepoint_levels))
  
  missing_terms <- setdiff(timepoint_levels, unique(prot_data$term))
  
  dummy_rows <- tibble(
    logFC       = NA_real_,
    p_value     = NA_real_,
    adj_p_value = NA_real_,
    fdr_sig2    = NA_character_,
    gene_symbol = NA_character_,  # match factor
    term        = factor(missing_terms, levels = timepoint_levels)
  )
  
  # 4. Combine real + dummy data and add significance label
  prot_data_full <- prot_data %>%
    mutate(fdr_sig2 = case_when(
      adj_p_value <= 0.05 ~ "FDR Significant",
      TRUE ~ "Not Significant"
    )) %>%
    bind_rows(dummy_rows)%>%
    filter(!is.na(term))
  
  
  
  highlighted_points <- prot_data %>%
    dplyr::filter(adj_p_value <= 0.05) %>%
    dplyr::group_by(term = factor(Timepoint,
                                  levels = names(custom_labels))) %>%
    dplyr::slice_min(order_by = adj_p_value, n = nlabels_volcano) %>%
    dplyr::ungroup()
  
  if (omes == "metabolomics" && GroupCode == "Endur") {
    lactic_acid_rows <- prot_data %>%
      dplyr::filter(gene_symbol == "Lactic acid", Timepoint == "post_10_min")
    
    if (nrow(lactic_acid_rows) > 0) {
      highlighted_points <- highlighted_points %>%
        dplyr::bind_rows(lactic_acid_rows) %>%
        dplyr::distinct(term, gene_symbol, .keep_all = TRUE)
    }
  }
  
  gray_overlay <- NULL
  if(length(missing_terms) > 0){
    gray_overlay <- tibble(
      term = factor(missing_terms, levels = timepoint_levels),
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    )
  }
  
  partA_volcano <- ggplot2::ggplot(
    prot_data_full,
    aes(x = logFC, y = -log10(p_value), color = fdr_sig2)
  ) +
    # Gray overlay behind empty panels
    {if(!is.null(gray_overlay)) geom_rect(
      data = gray_overlay,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "gray95"
    )} +
    geom_point(size = 2, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = highlighted_points, aes(label = gene_symbol),
      size = 2.5, color = "black", max.overlaps = Inf
    ) +
    geom_hline(yintercept = -log10(0.05), color = '#C7C7C7', size = 1, alpha = 0.8) +
    scale_color_manual(name = NULL,values = c("FDR Significant" = bars_color, "Not Significant" = "#C7C7C7"),
                       breaks = c("FDR Significant", "Not Significant")) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10),
      limits = c(0, ymax)
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    ggplot2::theme(panel.grid.major = element_blank()) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.background = element_rect(colour = "transparent", fill = "transparent"),
      axis.text.y = element_text(face = 'bold', size = 15),
      axis.text.x = element_text(face = 'bold', size = 15),
      axis.title = element_text(face = "bold", size = 17),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      strip.text = element_text(face = "bold", size = 13)
    ) +
    ggplot2::labs(x = "logFC", y = "-log10(p-value)") +
    ggplot2::guides(colour = guide_legend(override.aes = list(shape = 19, size = 7))) +
    ggplot2::facet_grid(~term, labeller = as_labeller(custom_labels), drop = FALSE)
  
  # ----------------------------- #
  #   4️⃣ Part B: ComplexUpset Plot
  # ----------------------------- #
  da_res_sig <- prot_data %>% dplyr::filter(adj_p_value < 0.05)
  
  desired_order <- c(
    "assay", "gene", "RE_post_24_hr", "RE_post_3.5_4_hr", "RE_post_15_30_45_min", "RE_post_10_min",
    "EE_post_24_hr", "EE_post_3.5_4_hr", "EE_post_15_30_45_min", "EE_post_10_min",
    "EE_during_40_min", "EE_during_20_min"
  )
  
  reshaped_data <- da_res_sig %>%
    tidyr::unite("group_timepoint", randomGroupCode, Timepoint, remove = FALSE) %>%
    dplyr::distinct(assay, gene, group_timepoint) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::spread(group_timepoint, value, fill = 0) %>%
    dplyr::select(any_of(desired_order))
  
  
  rename_map <- c(
    "EE_during_20_min"    = "During 20 Min",
    "EE_during_40_min"    = "During 40 Min",
    "EE_post_10_min"      = "Post 10 Min",
    "EE_post_15_30_45_min"= "Post 30 Min",
    "EE_post_3.5_4_hr"    = "Post 3.5 Hour",
    "RE_post_10_min"      = "Post 10 Min",
    "RE_post_15_30_45_min"= "Post 30 Min",
    "RE_post_3.5_4_hr"    = "Post 3.5 Hour",
    "RE_post_24_hr"       = "Post 24 Hour"
  )
  
  for (old_name in names(rename_map)) {
    if (old_name %in% colnames(reshaped_data)) {
      colnames(reshaped_data)[colnames(reshaped_data) == old_name] <- rename_map[[old_name]]
    }
  }
  
  
  
  group_timepoints <- colnames(reshaped_data)[3:ncol(reshaped_data)]
  reshaped_data[group_timepoints] <- reshaped_data[group_timepoints] == 1
  
  
  partB_upset <- ComplexUpset::upset(
    reshaped_data, group_timepoints, min_size=minsize_upset, name = '',
    base_annotations = list(
      'Intersection size' = ComplexUpset::intersection_size(
        text_mapping = aes(label = !!ComplexUpset::get_size_mode('exclusive_intersection')),
        bar_number_threshold = 1,
        text = list(vjust = -0.5, fontface = "bold", size = 5),
        counts = TRUE,
        mapping = aes(fill = 'bars_color')
      ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
        scale_fill_manual(values = c('bars_color' = bars_color), guide = 'none')+
        theme(axis.text.y = element_text(size = 14, face = "bold"))
    ),
    width_ratio = 0.1,
    set_sizes = (
      ComplexUpset::upset_set_size(
        geom = geom_bar(aes(fill = 'bars_color', x = group), width = 0.8),
        position = 'left'
      ) +
        scale_fill_manual(values = c('bars_color' = bars_color)) +
        theme(
          legend.position = "none",
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
          axis.title.x = element_text(size = 14, face = "bold")
        )
    ),
    sort_sets = FALSE
  ) +
    guides(color = "none") +
    theme(axis.text.y = element_text(size = 14, face = "bold"))
  
  # Derive individual file names
  pdf_dir <- dirname(pdf_output_path)
  if (!dir.exists(pdf_dir)) {
    dir.create(pdf_dir, recursive = TRUE)
  }
  base_name <- tools::file_path_sans_ext(basename(pdf_output_path))
  
  volcano_path <- file.path(pdf_dir, paste0(base_name, "_volcano.pdf"))
  upset_path   <- file.path(pdf_dir, paste0(base_name, "_upset.pdf"))
  combined_path <- pdf_output_path
  
  # --- Volcano-only ---
  grDevices::pdf(volcano_path, width = 10, height = 5)
  print(partA_volcano)
  grDevices::dev.off()
  
  # --- Upset-only ---
  grDevices::pdf(upset_path, width = 10, height = 5)
  print(partB_upset)
  grDevices::dev.off()
  
  # --- Combined PDF ---
  
  library(patchwork)
  
  combined_plot <- partA_volcano + partB_upset + plot_layout(ncol = 2, widths = c(1, 1))
  ggsave(combined_path, combined_plot, width = 20, height = 5)
  
  message("✅ Saved all figures:")
  message(" - Volcano-only: ", volcano_path)
  message(" - Upset-only:   ", upset_path)
  message(" - Combined:     ", combined_path)
  
  message("✅ Combined Volcano + Upset PDF saved to: ", pdf_output_path)
}