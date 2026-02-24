plot_metabolite_clinical_canonical_heatmap_pdf <- function(
    gsutil_cmd = "gsutil",
    permutations = 1000,
    repo_local_dir = NULL,
    cv_number,
    pdf_output_path
) {
  # --- 1. Load required libraries ---
  required_libs <- c(
    "dplyr", "ggplot2", "tidyverse", "glue", "readxl",
    "MotrpacHumanPreSuspension", "MotrpacHumanPreSuspensionData",
    "openxlsx", "PMA", "ComplexHeatmap", "circlize", "RColorBrewer", "tibble"
  )
  missing_libs <- required_libs[!sapply(required_libs, requireNamespace, quietly = TRUE)]
  if (length(missing_libs) > 0) {
    stop(paste0("Missing packages: ", paste(missing_libs, collapse = ", "),
                ". Please install them first."))
  }
  suppressPackageStartupMessages(lapply(required_libs, library, character.only = TRUE))
  
  
  plot_dir  <- dirname(pdf_output_path)
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    message(glue::glue("📂 Created Plot output directory: {plot_dir}"))
  }
  
  # --- 2. Load clinical data ---
  anthro_dataset <- MotrpacHumanPreSuspensionData::cln_curated_anthropometrics_vitals[["data"]] %>%
    filter(visit_code == "ADU_SCP")
  exercise_dataset <- MotrpacHumanPreSuspensionData::cln_curated_ex_performance_testing[["data"]] %>%
    filter(visit_code == "ADU_SCP")
  key_dataset <- MotrpacHumanPreSuspensionData::cln_raw_key[["data"]]
  
  clinical_data <- key_dataset %>%
    merge(anthro_dataset, by = "pid") %>%
    merge(exercise_dataset, by = c("pid", "visit_code")) %>%
    mutate(
      sex = case_when(
        sex_psca == 1 ~ "Female",
        sex_psca == 2 ~ "Male",
        TRUE ~ NA_character_
      ),
      age = age_psca
    )
  
  # --- 3. Load metabolomics omes data ---
  omes_data <- MotrpacHumanPreSuspensionData::load_qc(
    selected_tissues = "blood",
    selected_omes = "all",
    epigen = FALSE,
    gsutil = gsutil_cmd,
    repo_local_dir = repo_local_dir,
    load_acute_only = TRUE,
    remove_unnamed_metab = TRUE,
    remove_redundant_metab = TRUE,
    verbose = TRUE
  )
  
  # Combine metabolomics sample metadata and qc_norm
  metab_keys <- names(omes_data[["blood"]])
  metab_keys <- metab_keys[grepl("^metab", metab_keys)]
  sample_metadata_list <- lapply(metab_keys, function(k) omes_data[["blood"]][[k]][["sample_metadata"]])
  qc_norm_list <- lapply(metab_keys, function(k) omes_data[["blood"]][[k]][["qc_norm"]])
  
  metabolomics_sample_metadata <- dplyr::bind_rows(sample_metadata_list)
  metabolomics_data <- dplyr::bind_rows(qc_norm_list)
  
  metabolomics_sample_metadata <- metabolomics_sample_metadata %>%
    filter(Timepoint == "pre_exercise") %>%
    dplyr::select(vialLabel, pid, randomGroupCode)
  
  metabolomics_vial_labels <- metabolomics_sample_metadata$vialLabel
  metabolomics_subset <- metabolomics_data[, colnames(metabolomics_data) %in% metabolomics_vial_labels]
  metabolomics_subset <- metabolomics_subset[, order(match(colnames(metabolomics_subset), metabolomics_vial_labels))]
  
  metabolomics_bids <- setNames(metabolomics_sample_metadata$pid, metabolomics_sample_metadata$vialLabel)
  colnames(metabolomics_subset) <- metabolomics_bids[colnames(metabolomics_subset)]
  
  metabolomics_pids <- unique(gsub("\\..*", "", colnames(metabolomics_subset)))
  
  # --- 4. Prepare clinical subset ---
  clinical_subset <- clinical_data %>%
    dplyr::select(pid, bmi, peaktorqavg_iske, avg_handgrip, wtkg_cpet, wattend_cpet, hrend_cpet,
                  sysend_cpet, diaend_cpet, sysrecv1_cpet, diarecv1_cpet, hrrecv1_cpet,
                  avg_vecart_cpet, avg_vco2cart_cpet, avg_rercart_cpet,
                  VO2max, VO2max_L, O2_pulse, age) %>%
    filter(complete.cases(VO2max)) %>%
    filter(complete.cases(.))
  
  colnames(metabolomics_subset) <- make.unique(colnames(metabolomics_subset)) 
  
  # Step 2: Convert the wide format into a long format
  metabolomics_long <- metabolomics_subset %>%
    as.data.frame() %>%
    rownames_to_column(var = "metabolite")  %>%
    mutate(metabolite = gsub("\\.\\.\\.[0-9]+$", "", metabolite)) %>%
    distinct(metabolite, .keep_all = TRUE) %>%
    pivot_longer(
      cols = -metabolite,
      names_to = "pid",
      values_to = "value"
    )%>%
    mutate(original_pid = gsub("\\..*", "", pid)) %>%
    dplyr::select(-pid) %>%
    dplyr::rename(pid = original_pid) %>%
    filter(complete.cases(.))
  
  # Step 3: Coalesce data for each unique PID
  metabolomics_coalesced <- metabolomics_long %>%
    group_by(metabolite, pid) %>%
    summarize(
      value = dplyr::first(na.omit(value)),  # Take the first non-NA value for each group
      .groups = "drop"
    )
  
  # Step 4: Convert back to wide format
  metabolomics_wide <- metabolomics_coalesced %>%
    pivot_wider(
      names_from = pid,
      values_from = value
    ) %>%
    column_to_rownames("metabolite")
  
  metabolomics_wide <- metabolomics_wide[, colSums(is.na(metabolomics_wide)) == 0]
  
  metabolomics_subset = metabolomics_wide
  
  # --- 5. Align datasets ---
  common_bids <- intersect(colnames(metabolomics_subset), clinical_subset$pid)
  metabolomics_final <- metabolomics_subset[, common_bids]
  
  metabolomics_final <- as.data.frame(t(apply(metabolomics_final, 1, function(row) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
    row
  })))
  
  clinical_final <- clinical_subset[match(common_bids, clinical_subset$pid), ]
  clinical_final <- clinical_final %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "pid")
  
  metabolomics_final <- t(metabolomics_final)
  clinical_final <- as.matrix(clinical_final)
  
  # --- 6. CCA analysis ---
  set.seed(123)
  
  invisible(capture.output(
    permute_result <- PMA::CCA.permute(metabolomics_final, clinical_final,
                                       typex = "standard", typez = "standard", nperms = permutations)
  ))
  
  cca_result <- PMA::CCA(metabolomics_final, clinical_final,
                         typex = "standard", typez = "standard",
                         penaltyx = permute_result$bestpenaltyx,
                         penaltyz = permute_result$bestpenaltyz,
                         K = 3)
  
  metabolite_weights <- cca_result$u[, cv_number] # Canonical variates for metabolites
  clinical_weights <- cca_result$v[, cv_number]   # Canonical variates for clinical variables
  
  # Create a heatmap of canonical weights
  heatmap_data <- outer(clinical_weights, metabolite_weights)
  
  # Set the row and column names based on clinical and metabolomics data
  rownames(heatmap_data) <- colnames(clinical_final)
  colnames(heatmap_data) <- colnames(metabolomics_final)
  
  rownames(heatmap_data)[rownames(heatmap_data) == "VO2max"] <- "VO2max (ml/kg*min)"
  rownames(heatmap_data)[rownames(heatmap_data) == "VO2max_L"] <- "VO2max (L/min)"
  rownames(heatmap_data)[rownames(heatmap_data) == "O2_pulse"] <- "O2 pulse"
  rownames(heatmap_data)[rownames(heatmap_data) == "wattend_cpet"] <- "End Watts"
  rownames(heatmap_data)[rownames(heatmap_data) == "avg_handgrip"] <- "Handgrip Strength"
  rownames(heatmap_data)[rownames(heatmap_data) == "peaktorqavg_iske"] <- "Peak Torque"
  rownames(heatmap_data)[rownames(heatmap_data) == "avg_vco2cart_cpet"] <- "Avg VCO2"
  rownames(heatmap_data)[rownames(heatmap_data) == "avg_vecart_cpet"] <- "Avg VE"
  rownames(heatmap_data)[rownames(heatmap_data) == "diaend_cpet"] <- "DBP end"
  rownames(heatmap_data)[rownames(heatmap_data) == "diarecv1_cpet"] <- "DBP 1 min"
  rownames(heatmap_data)[rownames(heatmap_data) == "sysend_cpet"] <- "SBP end"
  rownames(heatmap_data)[rownames(heatmap_data) == "sysrecv1_cpet"] <- "SBP 1 min"
  rownames(heatmap_data)[rownames(heatmap_data) == "hrend_cpet"] <- "HR end"
  rownames(heatmap_data)[rownames(heatmap_data) == "hrrecv1_cpet"] <- "HR 1 min"
  rownames(heatmap_data)[rownames(heatmap_data) == "avg_rercart_cpet"] <- "Avg RER"
  rownames(heatmap_data)[rownames(heatmap_data) == "wtkg_cpet"] <- "Weight"
  rownames(heatmap_data)[rownames(heatmap_data) == "age"] <- "Age"
  rownames(heatmap_data)[rownames(heatmap_data) == "bmi"] <- "BMI"
  
  # Filter to include only non-zero rows and columns
  filtered_rows <- heatmap_data[rowSums(heatmap_data != 0) > 0, ]
  filtered_heatmap_data <- filtered_rows[, colSums(filtered_rows != 0) > 0]
  
  # Load necessary libraries
  library(ComplexHeatmap)
  library(circlize)
  
  # Order metabolites and clinical variables based on weights
  # Adjust the order of metabolites and clinical variables
  names(metabolite_weights) <- colnames(heatmap_data)
  metabolite_order <- order(metabolite_weights[colnames(filtered_heatmap_data)])
  filtered_heatmap_data <- filtered_heatmap_data[, metabolite_order]
  metabolite_weights_sorted <- metabolite_weights[colnames(filtered_heatmap_data)]
  
  names(clinical_weights) <- rownames(heatmap_data)
  clinical_order <- order(clinical_weights[rownames(filtered_heatmap_data)])
  filtered_heatmap_data <- filtered_heatmap_data[clinical_order, ]
  clinical_weights_sorted <- clinical_weights[rownames(filtered_heatmap_data)]
  
  top_metabolites <- names(sort(abs(metabolite_weights), decreasing = TRUE))[1:15]
  filtered_metabolite_weights <- metabolite_weights[top_metabolites]
  
  # Filter the heatmap data to include only these metabolites
  filtered_heatmap_data <- filtered_heatmap_data[, top_metabolites]
  
  # Update the metabolite weights to match the filtered data
  metabolite_weights_sorted <- filtered_metabolite_weights
  
  # Order the metabolites based on their weights
  metabolite_order <- order(metabolite_weights_sorted)
  filtered_heatmap_data <- filtered_heatmap_data[, metabolite_order]
  metabolite_weights_sorted <- metabolite_weights_sorted[metabolite_order]
  
  # Reorder the clinical variables (optional; can still use full set)
  names(clinical_weights) <- rownames(heatmap_data)
  clinical_order <- order(clinical_weights[rownames(filtered_heatmap_data)])
  filtered_heatmap_data <- filtered_heatmap_data[clinical_order, ]
  clinical_weights_sorted <- clinical_weights[rownames(filtered_heatmap_data)]
  
  # Create updated annotations for metabolites and clinical variables
  col_annotation <- HeatmapAnnotation(
    " " = anno_barplot(
      metabolite_weights_sorted,
      gp = gpar(fill = ifelse(metabolite_weights_sorted > 0, "red", "blue")),
      border = TRUE,
      axis_param = list(gp = gpar(fontsize = 8)),
      width = unit(2, "cm")
    )
  )
  
  row_annotation <- rowAnnotation(
    " " = anno_barplot(
      clinical_weights_sorted,
      gp = gpar(fill = ifelse(clinical_weights_sorted > 0, "red", "blue")),
      border = TRUE,
      axis_param = list(gp = gpar(fontsize = 8)),
      height = unit(2, "cm")
    )
  )
  
  # Create the heatmap with top 15 metabolites
  canonical_variate_1_heatmap_plot <- Heatmap(
    filtered_heatmap_data,
    name = "Beta", # The name of the heatmap (color legend)
    col = colorRamp2(c(-0.15, 0, 0.15), c("blue", "white", "red")), # Color scale
    bottom_annotation = col_annotation,
    right_annotation = row_annotation,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_heatmap_legend = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8)
  )
  pdf(pdf_output_path, width = 10, height = 8)
  draw(canonical_variate_1_heatmap_plot)
  dev.off()
  
  message("✅ Canonical variate heatmap saved at: ", pdf_output_path)
  return(canonical_variate_1_heatmap_plot)
}