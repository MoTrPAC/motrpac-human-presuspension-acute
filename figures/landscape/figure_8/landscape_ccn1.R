library(ComplexHeatmap)  # For heatmap visualization
library(circlize)        # For color mapping in heatmaps
library(dplyr)           # For data manipulation
library(tibble)          # For handling row names as columns
library(ggplot2)         # For additional plotting (if needed)
library(tidyr)           # For reshaping data
library(gridExtra)       # For arranging multiple heatmaps
library(magrittr)        # For using the pipe (%>%) operator
library(GSEABase)        # Required for gene set enrichment analysis
library(clusterProfiler) # For enrichment analysis (e.g., enrichmap)
library(Cairo)
library(MotrpacHumanPreSuspensionData)
library(MotrpacHumanPreSuspensionAnalysis)
library(TMSig)
library(purrr)


figures_output_folder = file.path(here(), "data", "tmp")

################################################################################
### CCN1 correlation with all genes within tissue (adipose or muscle)
################################################################################
#CCN1 within adipose.----------------------------------------------------------------
motrpac_qc = load_qc(selected_omes = "all",
                     selected_tissues = "adipose")
adi_trans_qc <- motrpac_qc$adipose$`transcript-rna-seq`$qc_norm
adi_trans_meta <- motrpac_qc$adipose$`transcript-rna-seq`$sample_metadata #meta

aid_pre_t <- adi_trans_meta %>%
  filter(Timepoint == "pre_exercise",
         visitcode == "ADU_BAS") %>%
  select(vialLabel)
aid_pre_t_id <- aid_pre_t$vialLabel
adi_trans_qc_pre <- adi_trans_qc[, colnames(adi_trans_qc) %in% aid_pre_t_id] # 172 pre exercise samples. We are missing three? -> 3 didnt get adipose biopsies
adi_trans_qc_pre <- adi_trans_qc_pre %>%
  rownames_to_column(var = "feature_id") %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% select(feature_id, gene_symbol), by = "feature_id") %>%
  filter(!is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>%
  slice(which.max(rowSums(across(where(is.numeric))))) %>% #chooses most highly expressed in cases of conflicts for gene_symbol values
  ungroup() %>%
  select(-feature_id) %>%
  column_to_rownames(var = "gene_symbol")

vial_to_pid <- adi_trans_meta %>%
  select(vialLabel, pid) %>%
  drop_na() %>%
  distinct() %>%
  deframe()  # Converts to named vector


adi_t <- adi_trans_meta %>%
  filter(visitcode == "ADU_BAS") %>%
  select(vialLabel)
aid_t_id <- adi_t$vialLabel
adi_trans_qc <- adi_trans_qc[, colnames(adi_trans_qc) %in% aid_t_id]
ccn1_feature <- HUMAN_FEATURE_TO_GENE %>%
  filter(gene_symbol == "CCN1") %>%
  select(feature_id)
print(ccn1_feature)

adi_trans_meta_filtered <- adi_trans_meta %>%
  filter(vialLabel %in% colnames(adi_trans_qc))


timepoint_results <- list()

# Loop through each time point
for (tp in unique(adi_trans_meta_filtered$Timepoint)) {
  # Subset samples for this time point
  tp_samples <- adi_trans_meta_filtered %>%
    filter(Timepoint == tp) %>%
    pull(vialLabel) %>%
    as.character()  # Ensure tp_samples is character

  # Subset the gene expression matrix
  tp_data <- adi_trans_qc[, tp_samples, drop = FALSE]

  # Extract CCN1 expression
  ccn1_expression <- as.numeric(tp_data["ENSG00000142871.18", , drop = TRUE])  # Ensure numeric

  # Transpose tp_data to match dimensions
  tp_data_transposed <- t(tp_data)

  # Use bicorAndPvalue to compute bicorrelation and p-values
  tp_results <- bicorAndPvalue(tp_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

  # Store bicorrelation and p-values in a list
  timepoint_results[[tp]] <- list(
    bicor = tp_results$bicor,
    pvalue = tp_results$p
  )
}
ccn1_adi_bicor <- do.call(cbind, lapply(timepoint_results, function(res) res$bicor))
ccn1_adi_pval <- do.call(cbind, lapply(timepoint_results, function(res) res$p))

colnames(ccn1_adi_bicor) <- names(timepoint_results)
colnames(ccn1_adi_pval) <- names(timepoint_results)

ccn1_adi_bicor <- ccn1_adi_bicor[, c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")]
save(ccn1_adi_bicor, ccn1_adi_pval, file = "ccn1_adi_bicor.RData")

# Initialize a list to store results
grouped_timepoint_results <- list()

# Loop through each time point
for (tp in unique(adi_trans_meta_filtered$Timepoint)) {
  # Subset samples for this time point
  tp_samples <- adi_trans_meta_filtered %>%
    filter(Timepoint == tp)

  # Loop through each group within the time point
  for (group in unique(tp_samples$randomGroupCode)) {
    # Subset samples for this group
    group_samples <- tp_samples %>%
      filter(randomGroupCode == group) %>%
      pull(vialLabel) %>%
      as.character()  # Ensure group_samples is character

    # Subset the gene expression matrix
    group_data <- adi_trans_qc[, group_samples, drop = FALSE]

    # Extract CCN1 expression
    ccn1_expression <- as.numeric(group_data["ENSG00000142871.18", , drop = TRUE])  # Ensure numeric

    # Transpose group_data to match dimensions
    group_data_transposed <- t(group_data)

    # Use bicorAndPvalue to compute bicorrelation and p-values
    group_results <- bicorAndPvalue(group_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

    # Store results in a nested list
    grouped_timepoint_results[[tp]][[group]] <- list(
      bicor = group_results$bicor,
      pvalue = group_results$p
    )
  }
}

save(grouped_timepoint_results, file = "adipose_ccn1_intracorrelations.RData")
# Output structure:
# grouped_timepoint_results[[timepoint]][[group]]$bicor
# grouped_timepoint_results[[timepoint]][[group]]$pvalue
# Initialize an empty list to store combined data for each group
combined_results <- list()

# List of groups
groups <- c("ADUEndur", "ADUControl", "ADUResist")

# Loop through each group
for (group in groups) {
  # Extract and combine bicor and pvalue across all time points for the current group
  combined_results[[group]] <- list(
    bicor = do.call(cbind, lapply(grouped_timepoint_results, function(tp) tp[[group]]$bicor)),
    pvalue = do.call(cbind, lapply(grouped_timepoint_results, function(tp) tp[[group]]$pvalue))
  )
}

# Access combined bicor and pvalue for each group
adendur_combined <- combined_results[["ADUEndur"]]
aducontrol_combined <- combined_results[["ADUControl"]]
aduresist_combined <- combined_results[["ADUResist"]]

# Example: Access combined bicor matrix for ADUEndur
ee_bicor <- adendur_combined$bicor
re_bicor <- aduresist_combined$bicor
con_bicor <- aducontrol_combined$bicor

#need to rename colnames with correct timepoints
timepoints <- names(grouped_timepoint_results)
colnames(ee_bicor) <- timepoints
colnames(re_bicor) <- timepoints
colnames(con_bicor) <- timepoints

#### Generate bicor heatmap and CAMERA plot for adipose CCN1 intracorrelation ###

generate_heatmap_and_enrichmap <- function(bicor_matrix, output_file, enrichmap_file, title = "Heatmap",
                                           row_title = "Transcripts", column_title = "CCN1 Intercorrelation",
                                           enrichmap_top_terms = 10) {
  # Define desired order for heatmap columns
  ordered_levels <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

  # Reorder columns of bicor_matrix for heatmap
  bicor_matrix <- bicor_matrix[, ordered_levels, drop = FALSE]

  # Timepoint colors
  timepoint_colors <- c(
    "pre_exercise" = "#440154",
    "post_15_30_45_min" = "#4fa78a",
    "post_3.5_4_hr" = "#a9d440",
    "post_24_hr" = "#fde725"
  )

  # Add column annotation with reordered labels
  col_annotation <- columnAnnotation(
    Timepoint = factor(colnames(bicor_matrix), levels = ordered_levels),
    col = list(Timepoint = timepoint_colors),
    annotation_name_gp = gpar(fontsize = 10)
  )
  # Color scale for heatmap
  color_scale <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  # Generate the heatmap
  heatmap_plot <- Heatmap(
    bicor_matrix,
    name = "Correlation",
    col = color_scale,
    row_title = row_title,
    column_title = column_title,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 10),
    top_annotation = col_annotation, # Add the color blocks as annotations
    heatmap_legend_param = list(
      title = "Bicor",
      at = c(-1, 0, 1),
      labels = c("-1", "0", "1")
    )
  )

  # Save the heatmap as a file
  png(output_file, width = 1100, height = 800, res = 300)
  draw(heatmap_plot)
  dev.off()

  # Rename columns in bicor_matrix for enrichmap
  enrichmap_colnames <- c(
    "pre_exercise" = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise",
    "post_15_30_45_min" = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise",
    "post_3.5_4_hr" = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise",
    "post_24_hr" = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise"
  )

  colnames(bicor_matrix) <- enrichmap_colnames[colnames(bicor_matrix)]

  # Prepare data for enrichmap
  cor_data <- bicor_matrix %>%
    as.data.frame() %>%
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id,
                 names_to = "contrast",
                 values_to = "rho") %>%
    mutate(z.std = 0.5 * log((1 + rho) / (1 - rho))) %>%  # Convert to z-scores
    filter(feature_id != "ENSG00000142871.18") %>%        # Exclude CCN1 itself
    list() %>%
    setNames("adipose.transcript-rna-seq")

  # Perform enrichmap analysis
  res <- run_cameraPR(DA_list = cor_data)

  # Rename contrasts for visualization
  res <- res %>%
    mutate(
      contrast_short = case_when(
        contrast == "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise" ~ "Pre",
        contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise" ~ "Post_30min",
        contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise" ~ "Post_4hr",
        contrast == "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise" ~ "Post_24hr",
        TRUE ~ contrast  # Keep other values as is
      ),
      contrast_short = factor(contrast_short, levels = c("Pre", "Post_30min", "Post_4hr", "Post_24hr")) # Set factor levels
    )

  # Select top terms
  top_terms <- res %>%
    filter(database == "REACTOME") %>%
    slice_min(order_by = p_value, by = contrast_short, n = enrichmap_top_terms) %>%
    pull(set_short) %>%
    unique()

  # Create enrichmap with reordered contrasts
  res %>%
    filter(set_short %in% top_terms) %>%
    enrichmap(n_top = Inf,
              set_column = "set_short",
              statistic_column = "z.std",
              contrast_column = "contrast_short",
              padj_column = "adj_p_value",
              padj_legend_title = "BH Adjusted\nP-Value",
              heatmap_args = list(heatmap_legend_param = list(
                title = "Z-Score"
              )),
              filename = enrichmap_file,
              height = 15,
              width = 10)
}
# Example Usage
generate_heatmap_and_enrichmap(
  bicor_matrix = con_bicor,
  output_file = file.path(figures_output_folder, "adi_ccn1_con_heatmap.png"),
  enrichmap_file = file.path(figures_output_folder, "adi_ccn1_con_enrichmap.png"),
  title = "CON CCN1 Heatmap",
  row_title = "Adipose transcripts",
  column_title = "CON CCN1 Intracorrelation",
  enrichmap_top_terms = 10
)

generate_heatmap_and_enrichmap(
  bicor_matrix = ee_bicor,
  output_file = file.path(figures_output_folder, "adi_ccn1_ee_heatmap.png"),
  enrichmap_file = file.path(figures_output_folder, "adi_ccn1_ee_enrichmap.png"),
  title = "EE CCN1 Heatmap",
  row_title = "Adipose transcripts",
  column_title = "EE CCN1 Intracorrelation",
  enrichmap_top_terms = 10
)

generate_heatmap_and_enrichmap(
  bicor_matrix = re_bicor,
  output_file = file.path(figures_output_folder, "adi_ccn1_re_heatmap.png"),
  enrichmap_file = file.path(figures_output_folder, "adi_ccn1_re_enrichmap.png"),
  title = "RE CCN1 Heatmap",
  row_title = "Adipose transcripts",
  column_title = "RE CCN1 Intracorrelation",
  enrichmap_top_terms = 10
)

######################### Muscle CCN1 Intracorrelation ##################################
mus_trans_qc <- motrpac_qc$muscle$`transcript-rna-seq`$qc_norm
mus_trans_meta <- motrpac_qc$muscle$`transcript-rna-seq`$sample_metadata #meta

mus_t <- mus_trans_meta %>%
  filter(visitcode == "ADU_BAS") %>%
  select(vialLabel)
mus_t_id <- mus_t$vialLabel
mus_trans_qc <- mus_trans_qc[, colnames(mus_trans_qc) %in% mus_t_id]

mus_trans_meta_filtered <- mus_trans_meta %>%
  filter(vialLabel %in% colnames(mus_trans_qc))

# Initialize a list to store results
grouped_timepoint_results_m <- list()

# Loop through each time point
for (tp in unique(mus_trans_meta_filtered$Timepoint)) {
  # Subset samples for this time point
  tp_samples <- mus_trans_meta_filtered %>%
    filter(Timepoint == tp)

  # Loop through each group within the time point
  for (group in unique(tp_samples$randomGroupCode)) {
    # Subset samples for this group
    group_samples <- tp_samples %>%
      filter(randomGroupCode == group) %>%
      pull(vialLabel) %>%
      as.character()  # Ensure group_samples is character

    # Subset the gene expression matrix
    group_data <- mus_trans_qc[, group_samples, drop = FALSE]

    # Extract CCN1 expression
    ccn1_expression <- as.numeric(group_data["ENSG00000142871.18", , drop = TRUE])  # Ensure numeric

    # Transpose group_data to match dimensions
    group_data_transposed <- t(group_data)

    # Use bicorAndPvalue to compute bicorrelation and p-values
    group_results <- bicorAndPvalue(group_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

    # Store results in a nested list
    grouped_timepoint_results_m[[tp]][[group]] <- list(
      bicor = group_results$bicor,
      pvalue = group_results$p
    )
  }
}
save(grouped_timepoint_results_m, file = "muscle_ccn1_intracorrelations.RData")

# Output structure:
# grouped_timepoint_results[[timepoint]][[group]]$bicor
# grouped_timepoint_results[[timepoint]][[group]]$pvalue
# Initialize an empty list to store combined data for each group
combined_results_m <- list()

# List of groups
groups <- c("ADUEndur", "ADUControl", "ADUResist")

# Loop through each group
for (group in groups) {
  # Extract and combine bicor and pvalue across all time points for the current group
  combined_results_m[[group]] <- list(
    bicor = do.call(cbind, lapply(grouped_timepoint_results_m, function(tp) tp[[group]]$bicor)),
    pvalue = do.call(cbind, lapply(grouped_timepoint_results_m, function(tp) tp[[group]]$pvalue))
  )
}

# Access combined bicor and pvalue for each group
ee_combined_m <- combined_results_m[["ADUEndur"]]
con_combined_m <- combined_results_m[["ADUControl"]]
re_combined_m <- combined_results_m[["ADUResist"]]

# Example: Access combined bicor matrix for ADUEndur
ee_bicor_m <- ee_combined_m$bicor
re_bicor_m <- re_combined_m$bicor
con_bicor_m <- con_combined_m$bicor

#need to rename colnames with correct timepoints
timepoints <- names(grouped_timepoint_results_m)
colnames(ee_bicor_m) <- timepoints
colnames(re_bicor_m) <- timepoints
colnames(con_bicor_m) <- timepoints


##Musle CCN1 heatmap and CAMERA plot -----------------------------------------
generate_heatmap_and_enrichmap_m <- function(bicor_matrix, output_file, enrichmap_file, title = "Heatmap",
                                           row_title = "Transcripts", column_title = "CCN1 Intracorrelation",
                                           enrichmap_top_terms = 10) {

  # Define desired order for heatmap columns
  ordered_levels <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

  # Reorder columns of bicor_matrix for heatmap
  bicor_matrix <- bicor_matrix[, ordered_levels, drop = FALSE]

  # Timepoint colors
  timepoint_colors <- c(
    "pre_exercise" = "#440154",
    "post_15_30_45_min" = "#4fa78a",
    "post_3.5_4_hr" = "#a9d440",
    "post_24_hr" = "#fde725"
  )

  # Add column annotation with reordered labels
  col_annotation <- columnAnnotation(
    Timepoint = factor(colnames(bicor_matrix), levels = ordered_levels),
    col = list(Timepoint = timepoint_colors),
    annotation_name_gp = gpar(fontsize = 10)
  )
  # Color scale for heatmap
  color_scale <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  # Generate the heatmap
  heatmap_plot <- Heatmap(
    bicor_matrix,
    name = "Correlation",
    col = color_scale,
    row_title = row_title,
    column_title = column_title,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 10),
    top_annotation = col_annotation, # Add the color blocks as annotations
    heatmap_legend_param = list(
      title = "Bicor",
      at = c(-1, 0, 1),
      labels = c("-1", "0", "1")
    )
  )

  # Save the heatmap as a file
  png(output_file, width = 1100, height = 800, res = 300)
  draw(heatmap_plot)
  dev.off()

  # Rename columns in bicor_matrix for enrichmap
  enrichmap_colnames <- c(
    "pre_exercise" = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise",
    "post_15_30_45_min" = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise",
    "post_3.5_4_hr" = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise",
    "post_24_hr" = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise"
  )

  colnames(bicor_matrix) <- enrichmap_colnames[colnames(bicor_matrix)]

  # Prepare data for enrichmap
  cor_data <- bicor_matrix %>%
    as.data.frame() %>%
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id,
                 names_to = "contrast",
                 values_to = "rho") %>%
    mutate(z.std = 0.5 * log((1 + rho) / (1 - rho))) %>%  # Convert to z-scores
    filter(feature_id != "ENSG00000142871.18") %>%        # Exclude CCN1 itself
    list() %>%
    setNames("muscle.transcript-rna-seq")

  # Perform enrichmap analysis
  res <- run_cameraPR(DA_list = cor_data)

  # Rename contrasts for visualization
  res <- res %>%
    mutate(
      contrast_short = case_when(
        contrast == "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise" ~ "Pre",
        contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise" ~ "Post_30min",
        contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise" ~ "Post_4hr",
        contrast == "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise" ~ "Post_24hr",
        TRUE ~ contrast  # Keep other values as is
      ),
      contrast_short = factor(contrast_short, levels = c("Pre", "Post_30min", "Post_4hr", "Post_24hr")) # Set factor levels
    )

  # Select top terms
  top_terms <- res %>%
    filter(database == "REACTOME") %>%
    slice_min(order_by = p_value, by = contrast_short, n = enrichmap_top_terms) %>%
    pull(set_short) %>%
    unique()

  # Create enrichmap with reordered contrasts
  res %>%
    filter(set_short %in% top_terms) %>%
    enrichmap(n_top = Inf,
              set_column = "set_short",
              statistic_column = "z.std",
              contrast_column = "contrast_short",
              padj_column = "adj_p_value",
              padj_legend_title = "BH Adjusted\nP-Value",
              heatmap_args = list(heatmap_legend_param = list(
                title = "Z-Score"
              )),
              filename = enrichmap_file,
              height = 15,
              width = 10)
}
generate_heatmap_and_enrichmap_m(
  bicor_matrix = con_bicor_m,
  output_file = file.path(figures_output_folder, "mus_ccn1_con_heatmap.png"),
  enrichmap_file = file.path(figures_output_folder, "mus_ccn1_con_enrichmap.png"),
  title = "CON CCN1 Heatmap",
  row_title = "Muscle transcripts",
  column_title = "CON CCN1 Intracorrelation",
  enrichmap_top_terms = 10
)

generate_heatmap_and_enrichmap_m(
  bicor_matrix = ee_bicor_m,
  output_file = file.path(figures_output_folder, "mus_ccn1_ee_heatmap.png"),
  enrichmap_file = file.path(figures_output_folder, "mus_ccn1_ee_enrichmap.png"),
  title = "EE CCN1 Heatmap",
  row_title = "Muscle transcripts",
  column_title = "EE CCN1 Intracorrelation",
  enrichmap_top_terms = 10
)

generate_heatmap_and_enrichmap_m(
  bicor_matrix = re_bicor_m,
  output_file = file.path(figures_output_folder, "mus_ccn1_re_heatmap.png"),
  enrichmap_file = file.path(figures_output_folder, "mus_ccn1_re_enrichmap.png"),
  title = "RE CCN1 Heatmap",
  row_title = "Muscle transcripts",
  column_title = "RE CCN1 Intracorrelation",
  enrichmap_top_terms = 10
)


############### Adipose CCN1 Intracorrelation - merged version ################
con_bicor_df <- as.data.frame(con_bicor) %>%
  rename_with(~ paste0("CON_", .)) %>%
  tibble::rownames_to_column("GeneID")  # Keep rownames for merging

ee_bicor_df <- as.data.frame(ee_bicor) %>%
  rename_with(~ paste0("EE_", .)) %>%
  tibble::rownames_to_column("GeneID")

re_bicor_df <- as.data.frame(re_bicor) %>%
  rename_with(~ paste0("RE_", .)) %>%
  tibble::rownames_to_column("GeneID")

# Merge by GeneID (rownames)
merged_bicor <- reduce(list(con_bicor_df, ee_bicor_df, re_bicor_df), full_join, by = "GeneID")

# Convert back to rownames format if needed
merged_bicor <- column_to_rownames(merged_bicor, var = "GeneID")
merged_bicor <- as.matrix(merged_bicor)

generate_merged_heatmap_and_enrichmap <- function(merged_bicor, output_file, enrichmap_file, title = "Merged Heatmap",
                                                  row_title = "Transcripts", column_title = "CCN1 Intercorrelation",
                                                  enrichmap_top_terms = 10) {

  # Define desired order for each group (CON, EE, RE)
  ordered_levels <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

  # Extract and reorder columns for each group
  con_cols <- paste0("CON_", ordered_levels)
  ee_cols <- paste0("EE_", ordered_levels)
  re_cols <- paste0("RE_", ordered_levels)

  # Ensure the columns exist
  valid_cols <- intersect(c(con_cols, ee_cols, re_cols), colnames(merged_bicor))
  bicor_matrix <- merged_bicor[, valid_cols, drop = FALSE]

  # Timepoint colors
  timepoint_colors <- c(
    "pre_exercise" = "#440154",
    "post_15_30_45_min" = "#4fa78a",
    "post_3.5_4_hr" = "#a9d440",
    "post_24_hr" = "#fde725"
  )

  # Group colors
  group_colors <- c("RE" = "#1b9e77", "EE" = "#d95f02", "CON" = "#7570b3")

  # Create column annotation for Timepoints and Exercise Groups
  col_annotation <- columnAnnotation(
    Group = factor(gsub("_.*", "", colnames(bicor_matrix)), levels = c("CON", "EE", "RE")),
    Timepoint = factor(sub("^(CON|EE|RE)_", "", colnames(bicor_matrix)), levels = ordered_levels),
    col = list(Group = group_colors, Timepoint = timepoint_colors),
    annotation_legend_param = list(title = "Exercise Group & Timepoint"),
    annotation_name_gp = gpar(fontsize = 10)
  )

  # Color scale for heatmap
  color_scale <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  # Generate the heatmap
  heatmap_plot <- Heatmap(
    bicor_matrix,
    name = "Correlation",
    col = color_scale,
    row_title = row_title,
    column_title = column_title,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 10),
    top_annotation = col_annotation, # Now includes Exercise Group colors
    heatmap_legend_param = list(
      title = "Bicor",
      at = c(-1, 0, 1),
      labels = c("-1", "0", "1")
    )
  )

  # Save the heatmap as a file
  pdf(output_file, width = 12, height = 8)
  draw(heatmap_plot)
  dev.off()

  # Rename columns in bicor_matrix for enrichmap
  enrichmap_colnames <- c(
    "CON_pre_exercise" = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise",
    "CON_post_15_30_45_min" = "group_timepointADUControl.post_15_30_45_min - group_timepointADUControl.pre_exercise",
    "CON_post_3.5_4_hr" = "group_timepointADUControl.post_3.5_4_hr - group_timepointADUControl.pre_exercise",
    "CON_post_24_hr" = "group_timepointADUControl.post_24_hr - group_timepointADUControl.pre_exercise",
    "EE_pre_exercise" = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise",
    "EE_post_15_30_45_min" = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
    "EE_post_3.5_4_hr" = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
    "EE_post_24_hr" = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
    "RE_pre_exercise" = "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise",
    "RE_post_15_30_45_min" = "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
    "RE_post_3.5_4_hr" = "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
    "RE_post_24_hr" = "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
  )

  colnames(bicor_matrix) <- enrichmap_colnames[colnames(bicor_matrix)]

  # Prepare data for enrichmap
  cor_data <- bicor_matrix %>%
    as.data.frame() %>%
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id,
                 names_to = "contrast",
                 values_to = "rho") %>%
    mutate(z.std = 0.5 * log((1 + rho) / (1 - rho))) %>%  # Convert to z-scores
    filter(feature_id != "ENSG00000142871.18") %>%        # Exclude CCN1 itself
    list() %>%
    setNames("adipose.transcript-rna-seq")

  # Perform enrichmap analysis
  res <- run_cameraPR(DA_list = cor_data)

  # Rename contrasts for visualization
  res <- res %>%
    mutate(
      contrast_short = case_when(
        contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise" ~ "CON Pre",
        contrast == "group_timepointADUControl.post_15_30_45_min - group_timepointADUControl.pre_exercise" ~ "CON 30min",
        contrast == "group_timepointADUControl.post_3.5_4_hr - group_timepointADUControl.pre_exercise" ~ "CON 4hr",
        contrast == "group_timepointADUControl.post_24_hr - group_timepointADUControl.pre_exercise" ~ "CON 24hr",
        contrast == "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise" ~ "EE Pre",
        contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise" ~ "EE 30min",
        contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise" ~ "EE 4hr",
        contrast == "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise" ~ "EE 24hr",
        contrast == "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise" ~ "RE Pre",
        contrast == "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise" ~ "RE 30min",
        contrast == "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise" ~ "RE 4hr",
        contrast == "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise" ~ "RE 24hr",
        TRUE ~ contrast
      ),
      contrast_short = factor(contrast_short, levels = c("CON Pre", "CON 30min", "CON 4hr", "CON 24hr",
                                                         "EE Pre", "EE 30min", "EE 4hr", "EE 24hr",
                                                         "RE Pre", "RE 30min", "RE 4hr", "RE 24hr")) # Set factor levels
    )

  # Select top terms
  top_terms <- res %>%
    filter(database == "REACTOME") %>%
    slice_min(order_by = p_value, by = contrast_short, n = enrichmap_top_terms) %>%
    pull(set_short) %>%
    unique()

  # Create enrichmap with reordered contrasts
  res %>%
    filter(set_short %in% top_terms) %>%
    enrichmap(n_top = Inf,
              set_column = "set_short",
              statistic_column = "z.std",
              contrast_column = "contrast_short",
              padj_column = "adj_p_value",
              padj_legend_title = "BH Adjusted\nP-Value",
              heatmap_args = list(heatmap_legend_param = list(
                title = "Z-Score"
              )),
              filename = enrichmap_file,
              height = 15,
              width = 10)
}

# Run the function with the merged bicor matrix
generate_merged_heatmap_and_enrichmap(
  merged_bicor = merged_bicor,
  output_file = file.path(figures_output_folder, "adipose_ccn1_intra_heatmap.pdf"),
  enrichmap_file = file.path(figures_output_folder, "adipose_ccn1_intra_camera.png"),
  title = "Adipose CCN1 Intracorrelation",
  row_title = "Adipose transcripts",
  column_title = "Adipose CCN1 Intracorrelation",
  enrichmap_top_terms = 6
)
############### Muscle CCN1 Intracorrelation - merged version ################

con_bicor_df_m <- as.data.frame(con_bicor_m) %>%
  rename_with(~ paste0("CON_", .)) %>%
  tibble::rownames_to_column("GeneID")  # Keep rownames for merging

ee_bicor_df_m <- as.data.frame(ee_bicor_m) %>%
  rename_with(~ paste0("EE_", .)) %>%
  tibble::rownames_to_column("GeneID")

re_bicor_df_m <- as.data.frame(re_bicor_m) %>%
  rename_with(~ paste0("RE_", .)) %>%
  tibble::rownames_to_column("GeneID")

# Merge by GeneID (rownames)
merged_bicor_m <- reduce(list(con_bicor_df_m, ee_bicor_df_m, re_bicor_df_m), full_join, by = "GeneID")

# Convert back to rownames format if needed
merged_bicor_m <- column_to_rownames(merged_bicor_m, var = "GeneID")
merged_bicor_m <- as.matrix(merged_bicor_m)

generate_merged_heatmap_and_enrichmap_m <- function(merged_bicor, output_file, enrichmap_file, title = "Merged Heatmap",
                                                  row_title = "Transcripts", column_title = "CCN1 Intercorrelation",
                                                  enrichmap_top_terms = 10) {

  # Define desired order for each group (CON, EE, RE)
  ordered_levels <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

  # Extract and reorder columns for each group
  con_cols <- paste0("CON_", ordered_levels)
  ee_cols <- paste0("EE_", ordered_levels)
  re_cols <- paste0("RE_", ordered_levels)

  # Ensure the columns exist
  valid_cols <- intersect(c(con_cols, ee_cols, re_cols), colnames(merged_bicor))
  bicor_matrix <- merged_bicor[, valid_cols, drop = FALSE]

  # Timepoint colors
  timepoint_colors <- c(
    "pre_exercise" = "#440154",
    "post_15_30_45_min" = "#4fa78a",
    "post_3.5_4_hr" = "#a9d440",
    "post_24_hr" = "#fde725"
  )

  # Group colors
  group_colors <- c("RE" = "#1b9e77", "EE" = "#d95f02", "CON" = "#7570b3")

  # Create column annotation for Timepoints and Exercise Groups
  col_annotation <- columnAnnotation(
    Group = factor(gsub("_.*", "", colnames(bicor_matrix)), levels = c("CON", "EE", "RE")),
    Timepoint = factor(sub("^(CON|EE|RE)_", "", colnames(bicor_matrix)), levels = ordered_levels),
    col = list(Group = group_colors, Timepoint = timepoint_colors),
    annotation_legend_param = list(title = "Exercise Group & Timepoint"),
    annotation_name_gp = gpar(fontsize = 10)
  )

  # Color scale for heatmap
  color_scale <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

  # Generate the heatmap
  heatmap_plot <- Heatmap(
    bicor_matrix,
    name = "Correlation",
    col = color_scale,
    row_title = row_title,
    column_title = column_title,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 10),
    top_annotation = col_annotation, # Now includes Exercise Group colors
    heatmap_legend_param = list(
      title = "Bicor",
      at = c(-1, 0, 1),
      labels = c("-1", "0", "1")
    )
  )

  # Save the heatmap as a file
  pdf(output_file, width = 12, height = 8)
  draw(heatmap_plot)
  dev.off()

  # Rename columns in bicor_matrix for enrichmap
  enrichmap_colnames <- c(
    "CON_pre_exercise" = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise",
    "CON_post_15_30_45_min" = "group_timepointADUControl.post_15_30_45_min - group_timepointADUControl.pre_exercise",
    "CON_post_3.5_4_hr" = "group_timepointADUControl.post_3.5_4_hr - group_timepointADUControl.pre_exercise",
    "CON_post_24_hr" = "group_timepointADUControl.post_24_hr - group_timepointADUControl.pre_exercise",
    "EE_pre_exercise" = "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise",
    "EE_post_15_30_45_min" = "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
    "EE_post_3.5_4_hr" = "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
    "EE_post_24_hr" = "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
    "RE_pre_exercise" = "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise",
    "RE_post_15_30_45_min" = "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
    "RE_post_3.5_4_hr" = "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
    "RE_post_24_hr" = "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
  )

  colnames(bicor_matrix) <- enrichmap_colnames[colnames(bicor_matrix)]

  # Prepare data for enrichmap
  cor_data <- bicor_matrix %>%
    as.data.frame() %>%
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id,
                 names_to = "contrast",
                 values_to = "rho") %>%
    mutate(z.std = 0.5 * log((1 + rho) / (1 - rho))) %>%  # Convert to z-scores
    filter(feature_id != "ENSG00000142871.18") %>%        # Exclude CCN1 itself
    list() %>%
    setNames("muscle.transcript-rna-seq")

  # Perform enrichmap analysis
  res <- run_cameraPR(DA_list = cor_data)

  # Rename contrasts for visualization
  res <- res %>%
    mutate(
      contrast_short = case_when(
        contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUResist.post_15_30_45_min + group_timepointADUResist.pre_exercise" ~ "CON Pre",
        contrast == "group_timepointADUControl.post_15_30_45_min - group_timepointADUControl.pre_exercise" ~ "CON 30min",
        contrast == "group_timepointADUControl.post_3.5_4_hr - group_timepointADUControl.pre_exercise" ~ "CON 4hr",
        contrast == "group_timepointADUControl.post_24_hr - group_timepointADUControl.pre_exercise" ~ "CON 24hr",
        contrast == "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise" ~ "EE Pre",
        contrast == "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise" ~ "EE 30min",
        contrast == "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise" ~ "EE 4hr",
        contrast == "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise" ~ "EE 24hr",
        contrast == "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise" ~ "RE Pre",
        contrast == "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise" ~ "RE 30min",
        contrast == "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise" ~ "RE 4hr",
        contrast == "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise" ~ "RE 24hr",
        TRUE ~ contrast
      ),
      contrast_short = factor(contrast_short, levels = c("CON Pre", "CON 30min", "CON 4hr", "CON 24hr",
                                                         "EE Pre", "EE 30min", "EE 4hr", "EE 24hr",
                                                         "RE Pre", "RE 30min", "RE 4hr", "RE 24hr")) # Set factor levels
    )

  # Select top terms
  top_terms <- res %>%
    filter(database == "REACTOME") %>%
    slice_min(order_by = p_value, by = contrast_short, n = enrichmap_top_terms) %>%
    pull(set_short) %>%
    unique()

  # Create enrichmap with reordered contrasts
  res %>%
    filter(set_short %in% top_terms) %>%
    enrichmap(n_top = Inf,
              set_column = "set_short",
              statistic_column = "z.std",
              contrast_column = "contrast_short",
              padj_column = "adj_p_value",
              padj_legend_title = "BH Adjusted\nP-Value",
              heatmap_args = list(heatmap_legend_param = list(
                title = "Z-Score"
              )),
              filename = enrichmap_file,
              height = 15,
              width = 10)
}

# Run the function with the merged bicor matrix
generate_merged_heatmap_and_enrichmap_m(
  merged_bicor = merged_bicor_m,
  output_file = file.path(figures_output_folder, "muscle_ccn1_intra_heatmap.pdf"),
  enrichmap_file = file.path(figures_output_folder, "muscle_ccn1_intra_camera.png"),
  title = "Muscle CCN1 Intracorrelation",
  row_title = "Muscle transcripts",
  column_title = "Muscle CCN1 Intracorrelation",
  enrichmap_top_terms = 6
)

##################### CCN1 cross tissue ####################################
#Adipose CCN1 & muscle all------------------------------------------------------
common_samples <- intersect(colnames(adi_trans_qc), colnames(mus_trans_qc))
adi_trans <- adi_trans_qc[, common_samples, drop = FALSE]
mus_trans <- mus_trans_qc[, common_samples, drop = FALSE]
# Define timepoints and groups
timepoints <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

# Extract CCN1 expression (ENSG00000142871.18) from adipose for each group
adi_con_ccn1 <- adi_trans %>% select(matches("^.*_ADUControl_")) %>% filter(rownames(.) == "ENSG00000142871.18")
adi_ee_ccn1  <- adi_trans %>% select(matches("^.*_ADUEndur_")) %>% filter(rownames(.) == "ENSG00000142871.18")
adi_re_ccn1  <- adi_trans %>% select(matches("^.*_ADUResist_")) %>% filter(rownames(.) == "ENSG00000142871.18")

# Match muscle columns for each group
mus_con <- mus_trans %>% select(matches("^.*_ADUControl_"))
mus_ee  <- mus_trans %>% select(matches("^.*_ADUEndur_"))
mus_re  <- mus_trans %>% select(matches("^.*_ADUResist_"))

# Create an empty list to store correlations by Group and Timepoint
group_timepoint_correlations <- list()

# Define groups and time points
groups <- c("ADUControl", "ADUEndur", "ADUResist")
timepoints <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

# Loop through each exercise group (CON, EE, RE)
for (group in groups) {

  # Initialize a list for the current group
  group_correlations <- list()

  # Loop through each time point
  for (tp in timepoints) {

    # Get adipose samples for the current group & time point
    adi_samples <- colnames(adi_trans)[grepl(paste0("_", group, "_", tp), colnames(adi_trans))]

    # Get muscle samples for the current group & time point
    mus_samples <- colnames(mus_trans)[grepl(paste0("_", group, "_", tp), colnames(mus_trans))]

    # Ensure shared PIDs between adipose and muscle
    shared_samples <- intersect(adi_samples, mus_samples)

    if (length(shared_samples) > 0) {
      # Subset CCN1 expression for adipose
      ccn1_expression <- adi_trans["ENSG00000142871.18", shared_samples, drop = TRUE]

      # Subset muscle data for the time point
      mus_data <- mus_trans[, shared_samples, drop = FALSE]

      # Compute bicorrelation between CCN1 and each muscle gene
      tp_cor <- apply(mus_data, 1, function(x) bicor(x, ccn1_expression, use = "pairwise.complete.obs"))

      # Store correlations with proper column name
      group_correlations[[tp]] <- tp_cor
    }
  }

  # Convert the correlations for the current group into a matrix
  group_cor_matrix <- do.call(cbind, group_correlations)

  # Add row and column names
  rownames(group_cor_matrix) <- rownames(mus_trans_qc)  # Muscle genes
  colnames(group_cor_matrix) <- paste0(group, "_", names(group_correlations))  # Label by Group & Timepoint

  # Store in the main list
  group_timepoint_correlations[[group]] <- group_cor_matrix
}

# Combine all group matrices into one dataframe
con_bicor_a2m <- group_timepoint_correlations[["ADUControl"]]
ee_bicor_a2m <- group_timepoint_correlations[["ADUEndur"]]
re_bicor_a2m <- group_timepoint_correlations[["ADUResist"]]

# Merge them together
a2m_ccn1 <- cbind(con_bicor_a2m, ee_bicor_a2m, re_bicor_a2m)
colnames(a2m_ccn1) <- gsub("ADUControl", "CON", colnames(a2m_ccn1))
colnames(a2m_ccn1) <- gsub("ADUEndur", "EE", colnames(a2m_ccn1))
colnames(a2m_ccn1) <- gsub("ADUResist", "RE", colnames(a2m_ccn1))

generate_merged_heatmap_and_enrichmap_m(
  merged_bicor = a2m_ccn1,
  output_file = file.path(figures_output_folder, "ccn1_a2m_heatmap.pdf"),
  enrichmap_file = file.path(figures_output_folder, "ccn1_a2m_camera.png"),
  title = "Adipose CCN1 - Muscle Correlation",
  row_title = "Muscle transcripts",
  column_title = "Adipose CCN1 - Muscle Correlation",
  enrichmap_top_terms = 6
)


# Muscle CCN1 & Adipose all ----------------------------------------------------

# Extract CCN1 expression (ENSG00000142871.18) from muscle for each group
mus_con_ccn1 <- mus_trans %>% select(matches("^.*_ADUControl_")) %>% filter(rownames(.) == "ENSG00000142871.18")
mus_ee_ccn1  <- mus_trans %>% select(matches("^.*_ADUEndur_")) %>% filter(rownames(.) == "ENSG00000142871.18")
mus_re_ccn1  <- mus_trans %>% select(matches("^.*_ADUResist_")) %>% filter(rownames(.) == "ENSG00000142871.18")

# Match muscle columns for each group
adi_con <- adi_trans %>% select(matches("^.*_ADUControl_"))
adi_ee  <- adi_trans %>% select(matches("^.*_ADUEndur_"))
adi_re  <- adi_trans %>% select(matches("^.*_ADUResist_"))

# Create an empty list to store correlations by Group and Timepoint
group_timepoint_correlations_m2a <- list()

# Define groups and time points
groups <- c("ADUControl", "ADUEndur", "ADUResist")
timepoints <- c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")

# Loop through each exercise group (CON, EE, RE)
for (group in groups) {

  # Initialize a list for the current group
  group_correlations <- list()

  # Loop through each time point
  for (tp in timepoints) {

    # Get adipose samples for the current group & time point
    adi_samples <- colnames(adi_trans)[grepl(paste0("_", group, "_", tp), colnames(adi_trans))]

    # Get muscle samples for the current group & time point
    mus_samples <- colnames(mus_trans)[grepl(paste0("_", group, "_", tp), colnames(mus_trans))]

    # Ensure shared PIDs between adipose and muscle
    shared_samples <- intersect(adi_samples, mus_samples)

    if (length(shared_samples) > 0) {
      # Subset CCN1 expression for adipose
      ccn1_expression <- mus_trans["ENSG00000142871.18", shared_samples, drop = TRUE]

      # Subset adipose data for the time point
      adi_data <- adi_trans[, shared_samples, drop = FALSE]

      # Compute bicorrelation between CCN1 and each muscle gene
      tp_cor <- apply(adi_data, 1, function(x) bicor(x, ccn1_expression, use = "pairwise.complete.obs"))

      # Store correlations with proper column name
      group_correlations[[tp]] <- tp_cor
    }
  }

  # Convert the correlations for the current group into a matrix
  group_cor_matrix <- do.call(cbind, group_correlations)

  # Add row and column names
  rownames(group_cor_matrix) <- rownames(adi_trans_qc)  # Adipose genes
  colnames(group_cor_matrix) <- paste0(group, "_", names(group_correlations))  # Label by Group & Timepoint

  # Store in the main list
  group_timepoint_correlations_m2a[[group]] <- group_cor_matrix
}

# Combine all group matrices into one dataframe
con_bicor_m2a <- group_timepoint_correlations_m2a[["ADUControl"]]
ee_bicor_m2a <- group_timepoint_correlations_m2a[["ADUEndur"]]
re_bicor_m2a <- group_timepoint_correlations_m2a[["ADUResist"]]

# Merge them together
m2a_ccn1 <- cbind(con_bicor_m2a, ee_bicor_m2a, re_bicor_m2a)
colnames(m2a_ccn1) <- gsub("ADUControl", "CON", colnames(m2a_ccn1))
colnames(m2a_ccn1) <- gsub("ADUEndur", "EE", colnames(m2a_ccn1))
colnames(m2a_ccn1) <- gsub("ADUResist", "RE", colnames(m2a_ccn1))

generate_merged_heatmap_and_enrichmap(
  merged_bicor = m2a_ccn1,
  output_file = file.path(figures_output_folder, "ccn1_m2a_heatmap.pdf"),
  enrichmap_file = file.path(figures_output_folder, "ccn1_m2a_camera.png"),
  title = "Muscle CCN1 - Adipose Correlation",
  row_title = "Adipose transcripts",
  column_title = "Muscle CCN1 - Adipose Correlation",
  enrichmap_top_terms = 6
)
