generate_heatmap_interactive3 <- function(
    repo_local_dir,
    secprot_path,
    supplemental_excel,
    figure_path,
    tissue_type = c("adipose", "muscle"),  # 👈 choose tissue type
    timepoint = NULL
) {
  tissue_type <- match.arg(tissue_type)
  
  # ---- Load required libraries ----
  required_libs <- c(
    "dplyr", "tidyr", "ComplexHeatmap", "circlize", "RColorBrewer", "grid",
    "openxlsx", "stringr", "purrr"
  )
  suppressPackageStartupMessages(lapply(required_libs, library, character.only = TRUE))
  
  plot_dir  <- dirname(figure_path)
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message(glue::glue("📂 Created Plot output directory: {plot_dir}"))
  }
  
  # ---- Helper ----
  extract_and_process <- function(data, gene_column, tissue_label, gene_data) {
    if (is.null(data) || nrow(data) == 0) return(NULL)
    data %>%
      filter(.data[[gene_column]] %in% gene_data$Genes) %>%
      dplyr::select(all_of(gene_column), adj_p_value, logFC, randomGroupCode, Timepoint) %>%
      dplyr::rename(feature_id = all_of(gene_column)) %>%
      mutate(Tissue_ome = tissue_label)
  }
  
  # ---- Step 1: Load DA data ----
  df_da_full <- load_differential_analysis(
    repo_local_dir = repo_local_dir,
    selected_omes = "all",
    selected_tissues = "all",
    single_matrix = TRUE,
    combine_with_featgene = TRUE,
    epigen = FALSE
  ) %>%
    filter(contrast_type == "exercise_with_controls") %>%
    mutate(
      assay = if_else(str_detect(assay, "metab"), "metab", assay),
      gene_symbol = if_else(assay == "metab", feature_id, gene_symbol)
    )
  
  color_exercise <- MotrpacHumanPreSuspension::HUMAN_EXERCISE_GROUP_COLORS
  color_resist <- color_exercise[["ADUResist"]]
  color_endur  <- color_exercise[["ADUEndur"]]
  color_control <- color_exercise[["ADUControl"]]
  
  # ---- Step 2: Split ----
  da_fil_data_bp <- df_da_full %>% filter(tissue == "blood", assay == "prot-ol")
  da_fil_data_br <- df_da_full %>% filter(tissue == "blood", assay == "transcript-rna-seq")
  da_fil_data_mr <- df_da_full %>% filter(tissue == "muscle", assay == "transcript-rna-seq")
  da_fil_data_ar <- df_da_full %>% filter(tissue == "adipose", assay == "transcript-rna-seq")
  da_fil_data_mp <- df_da_full %>% filter(tissue == "muscle", assay == "prot-pr")
  da_fil_data_ap <- df_da_full %>% filter(tissue == "adipose", assay == "prot-pr")
  
  # ---- Step 3: Secreted protein annotation ----
  secprotcell <- read.csv(secprot_path, check.names = FALSE)
  
  mapping <- MotrpacHumanPreSuspensionData::HUMAN_FEATURE_TO_GENE %>%
    filter(assay == "prot-ol") %>%
    dplyr::select(feature_id, uniprot, gene_symbol) %>%
    mutate(across(everything(), as.character))
  
  da_blood_prot_fil_data <- da_fil_data_bp %>%
    dplyr::select(feature_id, uniprot, gene_symbol, randomGroupCode, Timepoint, logFC, p_value, adj_p_value)
  
  da_sec_prot <- secprotcell %>%
    dplyr::select(name, "tissue label", "uniprot")
  
  da_blood_prot_fil_data_significant <- da_blood_prot_fil_data %>%
    filter(adj_p_value < 0.05) %>%
    mutate(in_sec_prot = uniprot %in% da_sec_prot$uniprot)
  
  da_blood_prot_fil_data_significant_in_sec_prot <- da_blood_prot_fil_data_significant %>%
    filter(in_sec_prot) %>%
    left_join(da_sec_prot, by = "uniprot")
  
  da_sec_prot_sub_da <- secprotcell[secprotcell$uniprot %in% da_blood_prot_fil_data_significant$uniprot,]
  da_sec_prot_sub_da <- merge(mapping, da_sec_prot_sub_da, by = "uniprot")
  
  df_long <- da_sec_prot_sub_da %>%
    dplyr::select(gene_symbol, `Global label`) %>%
    mutate(`Global label` = as.character(`Global label`)) %>%
    separate_rows(`Global label`, sep = "\\.") %>%
    mutate(`Global label` = str_trim(`Global label`)) %>%
    filter(`Global label` != "") %>%
    mutate(present = 1)
  
  gene_data <- read.xlsx(supplemental_excel, sheet = "Sheet3")
  
  df_long_for_assays <- df_long %>% mutate(assay_value_for_cells = gene_symbol)
  
  df_wide_assays <- df_long_for_assays %>%
    pivot_wider(
      id_cols = gene_symbol,
      names_from = `Global label`,
      values_from = assay_value_for_cells,
      values_fill = NA
    )
  
  tissue_columns <- setdiff(names(df_wide_assays), "gene_symbol")
  
  list_of_assays_by_tissue <- lapply(tissue_columns, function(tissue_name) {
    assay_vector <- df_wide_assays[[tissue_name]]
    assays_present <- assay_vector[!is.na(assay_vector)]
    unique(assays_present)
  })
  names(list_of_assays_by_tissue) <- tissue_columns
  
  max_length <- max(sapply(list_of_assays_by_tissue, length))
  padded_list <- lapply(list_of_assays_by_tissue, function(assay_vector) {
    padding_needed <- max_length - length(assay_vector)
    c(assay_vector, rep(NA, padding_needed))
  })
  df_tissues_as_cols <- as.data.frame(padded_list)
  
  gene_data_sheet3_adipose <- drop_na(data.frame(Genes = as.character(df_tissues_as_cols$adiposetissue)))
  gene_data_sheet3_muscle  <- drop_na(data.frame(Genes = as.character(df_tissues_as_cols$muscle)))
  
  # ---- Step 4: Pick tissue ----
  if (tissue_type == "adipose") {
    gene_data_sheet3 <- gene_data_sheet3_adipose
    df_combined <- bind_rows(
      extract_and_process(da_fil_data_bp, "gene_symbol", "blood_protein", gene_data_sheet3),
      extract_and_process(da_fil_data_ar, "gene_symbol", "adipose_rna", gene_data_sheet3),
      extract_and_process(da_fil_data_ap, "gene_symbol", "adipose_protein", gene_data_sheet3)
    )
  } else {
    gene_data_sheet3 <- gene_data_sheet3_muscle
    df_combined <- bind_rows(
      extract_and_process(da_fil_data_bp, "gene_symbol", "blood_protein", gene_data_sheet3),
      extract_and_process(da_fil_data_mr, "gene_symbol", "muscle_rna", gene_data_sheet3),
      extract_and_process(da_fil_data_mp, "gene_symbol", "muscle_protein", gene_data_sheet3)
    )
  }
  
  # ---- Optional timepoint filtering ----
  if (!is.null(timepoint)) {
    df_combined <- df_combined %>% filter(Timepoint == timepoint)
    message(paste("⏱ Filtering to Timepoint:", timepoint))
  }
  
  if (nrow(df_combined) == 0) stop("No data found for the specified filters.")
  
  # ---- Build matrices ----
  df_combined <- df_combined %>% mutate(row_id = paste(feature_id, Tissue_ome, sep = "|"))
  
  all_groups <- unique(df_combined$randomGroupCode)
  all_timepoints <- unique(df_combined$Timepoint)
  all_tissue_omes <- unique(df_combined$Tissue_ome)
  
  all_combos <- expand.grid(
    feature_id = gene_data_sheet3$Genes,
    Tissue_ome = all_tissue_omes,
    randomGroupCode = all_groups,
    Timepoint = all_timepoints,
    stringsAsFactors = FALSE
  ) %>%
    mutate(row_id = paste(feature_id, Tissue_ome, sep = "|"))
  
  # Merge with observed data
  df_combined <- all_combos %>%
    left_join(df_combined, by = c("feature_id","Tissue_ome","randomGroupCode","Timepoint","row_id"))
  
  logFC_mat <- df_combined %>%
    pivot_wider(id_cols = row_id,
                names_from = c(randomGroupCode, Timepoint),
                values_from = logFC) %>%
    tibble::column_to_rownames("row_id") %>%
    as.matrix()
  
  pval_mat <- df_combined %>%
    pivot_wider(id_cols = row_id,
                names_from = c(randomGroupCode, Timepoint),
                values_from = adj_p_value) %>%
    tibble::column_to_rownames("row_id") %>%
    as.matrix()
  
  logFC_mat <- logFC_mat[, colSums(!is.na(logFC_mat)) > 0, drop = FALSE]
  pval_mat  <- pval_mat[, colSums(!is.na(pval_mat)) > 0, drop = FALSE]
  
  # Column metadata
  col_parts <- strsplit(colnames(logFC_mat), "_")
  col_groups <- sapply(col_parts, `[`, 1)
  col_timepoints <- sapply(col_parts, function(x) paste(x[-1], collapse = "_"))
  
  # timepoint_labels <- c(
  #   "during_20_min" = "During 20 Min",
  #   "during_40_min" = "During 40 Min",
  #   "post_10_min" = "Post 10 Min",
  #   "post_15_30_45_min" = "Post 15/30 Min",
  #   "post_3.5_4_hr" = "Post 3.5 Hour",
  #   "post_24_hr" = "Post 24 Hour"
  # )
  
  if (tissue_type == "muscle") {
    
    timepoint_labels <- c(
      "during_20_min"        = "During 20 Min",
      "during_40_min"        = "During 40 Min",
      "post_10_min"          = "Post 10 Min",
      "post_15_30_45_min"    = "Post 15/30 Min",
      "post_3.5_4_hr"        = "Post 3.5 Hour",
      "post_24_hr"           = "Post 24 Hour"
    )
    
  } else if (tissue_type == "adipose") {
    
    timepoint_labels <- c(
      "during_20_min"        = "During 20 Min",
      "during_40_min"        = "During 40 Min",
      "post_10_min"          = "Post 10 Min",
      "post_15_30_45_min"    = "Post 15/45 Min",
      "post_3.5_4_hr"        = "Post 3.5/4 Hour",
      "post_24_hr"           = "Post 24 Hour"
    )
    
  } else {
    
    stop("Unknown tissue_type: must be 'adipose' or 'muscle'")
  }
  
  col_timepoints <- ifelse(col_timepoints %in% names(timepoint_labels),
                           timepoint_labels[col_timepoints],
                           col_timepoints)
  
  # time_colors <- setNames(
  #   c("#FDE725", "#BAD071", "#D1BBD7", "#AE76A3", "#882E72", "#61194F"),
  #   unname(timepoint_labels)
  # )
  
  motrpactimepoint_colors <- MotrpacHumanPreSuspension::HUMAN_ACUTE_TIMEPOINT_COLORS
  
  # Keep only the colors for the timepoints you have
  motrpactimepoint_colors <- motrpactimepoint_colors[names(motrpactimepoint_colors) %in% names(timepoint_labels)]
  
  # Rename color vector to friendly labels
  time_colors <- setNames(motrpactimepoint_colors[names(motrpactimepoint_colors) %in% names(timepoint_labels)],
                          timepoint_labels[names(motrpactimepoint_colors) %in% names(timepoint_labels)])
  
  # Tissue colors
  unique_tissues <- unique(df_combined$Tissue_ome)
  tissue_palette <- c(
    "blood_protein" = "#8dd3c7",
    "blood_rna"     = "#ffffb3",
    "muscle_rna"    = "#bebada",
    "adipose_rna"   = "#fb8072",
    "muscle_protein" = "#80b1d3",
    "adipose_protein" = "#fdb462" 
  )
  tissue_order <- unique_tissues
  tissue_colors <- tissue_palette[tissue_order]
  
  # Feature colors (hidden legend later)
  feature_order <- unique(gsub("\\|.*", "", rownames(logFC_mat)))
  feature_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(feature_order)), feature_order)
  
  # Column ordering
  column_order <- order(col_groups, match(col_timepoints, unname(timepoint_labels), nomatch = length(timepoint_labels) + 1))
  logFC_mat <- logFC_mat[, column_order]
  pval_mat <- pval_mat[, column_order]
  col_groups <- col_groups[column_order]
  col_timepoints <- col_timepoints[column_order]
  
  # Group mapping
  group_display <- c("ADUEndur" = "EE", "ADUResist" = "RE")
  group_colors <- c("EE" = color_endur, "RE" = color_resist)
  
  # Column annotations
  col_annot <- HeatmapAnnotation(
    Group = factor(group_display[col_groups], levels = c("EE", "RE")),
    Timepoint = factor(col_timepoints, levels = unname(timepoint_labels)),
    col = list(Group = group_colors, Timepoint = time_colors),
    # annotation_name_side = "left",
    show_annotation_name = FALSE,
    show_legend = FALSE,
    annotation_legend_param = list(
      Group = list(title = "Group", legend_gp = gpar(fontsize = 18)),
      Timepoint = list(title = "Timepoint", legend_gp = gpar(fontsize = 18))
    )
  )
  
  
  row_annot <- rowAnnotation(
    Tissue = factor(gsub(".*\\|", "", rownames(logFC_mat)), levels = tissue_order),
    col = list(Tissue = tissue_colors),
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    annotation_legend_param = list(
      Tissue = list(
        title = "Tissue",
        title_gp = gpar(fontsize = 14, fontface = "bold"),
        labels_gp = gpar(fontsize = 12)
      )
    ),
    show_legend = FALSE
  )
  
  
  all_na_rows <- rowSums(!is.na(logFC_mat)) == 0
  
  ht <- Heatmap(
    logFC_mat,
    col = colorRamp2(
      c(min(logFC_mat, na.rm = TRUE), 0, max(logFC_mat, na.rm = TRUE)),
      c("blue", "white", "red")
    ),
    column_split = col_groups,
    cluster_columns = FALSE,
    column_order = 1:ncol(logFC_mat),
    column_labels = NULL,
    column_names_rot = 75,
    row_split = factor(gsub("\\|.*", "", rownames(logFC_mat))),
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    row_labels = gsub("\\|.*", "", rownames(logFC_mat)),
    top_annotation = col_annot,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    left_annotation = row_annot,
    column_title = NULL,
    border = TRUE,
    heatmap_legend_param = list(title = "logFC"),
    show_heatmap_legend = FALSE,
    
    # ✅ Combined cell_fun
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.na(logFC_mat[i, j])) {
        if (all_na_rows[i]) {
          grid.rect(x, y, width, height, gp = gpar(fill = "#818589", col = NA))  # Not detected
        } else {
          grid.rect(x, y, width, height, gp = gpar(fill = "lightgray", col = NA)) # Not measured
        }
      }
      if (!is.na(pval_mat[i, j]) && pval_mat[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 12, col = "black"))
      }
    }
  )
  
  na_legend <- Legend(
    labels = c("Not measured", "Not detected"),
    legend_gp = gpar(fill = c("lightgray", "#818589")),
    title = "Missing Data"
  )
  
  lgd_group <- Legend(
    title = "Group",
    at = c("EE", "RE"),
    legend_gp = gpar(fill = c(color_endur, color_resist))
  )
  
  # Timepoint legend
  lgd_timepoint <- Legend(
    title = "Timepoint",
    at = names(time_colors),
    legend_gp = gpar(fill = unname(time_colors))
  )
  
  # Tissue legend
  lgd_tissue <- Legend(
    title = "Tissue",
    at = names(tissue_colors),
    legend_gp = gpar(fill = unname(tissue_colors))
  )
  
  # logFC legend
  lgd_logFC <- Legend(
    title = "logFC",
    col_fun = colorRamp2(
      c(min(logFC_mat, na.rm = TRUE), 0, max(logFC_mat, na.rm = TRUE)),
      c("blue", "white", "red")
    )
  )
  
  # Missing data legend
  lgd_missing <- Legend(
    title = "Missing Data",
    labels = c("Not measured", "Not detected"),
    legend_gp = gpar(fill = c("lightgray", "#818589"))
  )
  
  # Combine in order
  combined_legends <- packLegend(
    lgd_group,
    lgd_timepoint,
    lgd_tissue,
    lgd_logFC,
    lgd_missing,
    direction = "vertical"
  )
  
  
  # ---- Save output ----
  output_file <- file.path(
    figure_path
  )
  
  pdf(output_file, width = 12, height = 15)
  draw(
    ht,
    annotation_legend_list = list(combined_legends),
    heatmap_legend_list = NULL,
    merge_legends = TRUE,              # ✅ merge all legends into one vertical stack
    heatmap_legend_side = "right",     # ✅ legends to right side
    annotation_legend_side = "right"  # ✅ stack legends top-to-bottom
  )
  dev.off()
  
  message(paste("✅ Heatmap saved to:", output_file))
}