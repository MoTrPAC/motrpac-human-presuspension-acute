library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(stringr)
library(readr)
library(tibble)
library(tidyr)
library(circlize)
library(WGCNA)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(ComplexUpset)
library(ggplot2)
library(tidyverse)

motrpac_da <- load_differential_analysis()
motrpac_qc <- load_qc()


# Identify blood protein DA
blood_da_pr <- motrpac_da$blood$`prot-ol` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))

#Map Olink id with gene symbol using motrpac_qc
blood_prot_map <- motrpac_qc$blood$`prot-ol`$feature_metadata
blood_da_pr <- blood_da_pr %>%
  left_join(blood_prot_map %>% dplyr::select(feature_id, assay), by = "feature_id") %>%
  dplyr::rename(gene_symbol = assay.y) %>%
  filter(adj_p_value<0.1)

blood_da_pr_list <- unique(blood_da_pr$gene_symbol) #278 unique blood DA proteins at all time points, in EE & RE. ## These are the blood proteins we will use to screen for secretomes.

############## Pull out DA secretomes from each tissue (muscle, adipose, blood)
# 1. Muscle - transcripts
muscle_da_t <- motrpac_da$muscle$`transcript-rna-seq` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
muscle_da_t <- muscle_da_t %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
muscle_da_t2 <- muscle_da_t %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)

# 2. Muscle protein
muscle_da_pr <- motrpac_da$muscle$`prot-pr` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
muscle_da_pr <- muscle_da_pr %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
muscle_da_pr2 <- muscle_da_pr %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)

# 3. Muscle phosphosite. For phospho, we ignore the sites. Just pulling out the protein names from DA phosphosites.
muscle_da_ph <- motrpac_da$muscle$`prot-ph` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
muscle_da_ph <- muscle_da_ph %>%
  mutate(feature_id = sub("_.*", "", feature_id))
muscle_da_ph <- muscle_da_ph %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
muscle_da_ph2 <- muscle_da_ph %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)

# 4. Adipose - transcripts
adipose_da_t <- motrpac_da$adipose$`transcript-rna-seq` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
adipose_da_t <- adipose_da_t %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
adipose_da_t2 <- adipose_da_t %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)

# 5. Adipose - proteins
adipose_da_pr <- motrpac_da$adipose$`prot-pr` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
adipose_da_pr <- adipose_da_pr %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
adipose_da_pr2 <- adipose_da_pr %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)

# 6. Adipose - phosphosites. For phospho, we ignore the sites. Just pulling out the protein names from DA phosphosites.
adipose_da_ph <- motrpac_da$adipose$`prot-ph` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
adipose_da_ph <- adipose_da_ph %>%
  mutate(feature_id = sub("_.*", "", feature_id))
adipose_da_ph <- adipose_da_ph %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
adipose_da_ph2 <- adipose_da_ph %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)

# 7. Blood - transcript
blood_da_t <- motrpac_da$blood$`transcript-rna-seq` %>%
  filter(contrast_type == "exercise_with_controls" & contrast_category %in% c("EE-CON", "RE-CON"))
blood_da_t <- blood_da_t %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")
blood_da_t2 <- blood_da_t %>%
  mutate(
    blood_da = ifelse(gene_symbol %in% blood_da_pr_list, "Yes", "No"),
    Exercise = str_extract(contrast_category, "^(EE|RE)"),
    Direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  filter(adj_p_value<0.1)
### End of secretome profile from each tissue and ome.

## Combine them all and move forward with secondary analyses
combined_da <- bind_rows(
  muscle_da_t2,
  muscle_da_pr2,
  muscle_da_ph2,
  adipose_da_t2,
  adipose_da_pr2,
  adipose_da_ph2,
  blood_da_t2
)

### Call in COMPARTMENTS data.
compartments <-read_tsv(
  "https://download.jensenlab.org/human_compartment_integrated_full.tsv",
  col_names = c("feature_id", "gene_symbol", "GO", "location", "score"),
  show_col_types = FALSE
)
## This is the COMPARTMENT database that has extracellular scores.

target_locations <- c("Extracellular region", "Extracellular space", "Extracellular exosome", "Extracellular vesicle") # These are the 4 locations that has 'extracellular' regions.

# Filter the dataset for the selected genes and locations
compartments_filtered <- compartments %>%
  filter(gene_symbol %in% blood_da_pr$gene_symbol & location %in% target_locations)

# Pick the highest extracellular score among 4 locations.
compartments_filtered_max <- compartments_filtered %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup()

# Now, attach COMPARMENT information to our combined DA secretomes.
combined_da_fil <- combined_da %>%
  filter(gene_symbol %in% blood_da_pr$gene_symbol)
combined_da_fil <- combined_da_fil %>%
  group_by(gene_symbol) %>%
  mutate(Exercise2 = ifelse(n_distinct(Exercise) > 1, "Both", Exercise)) %>%
  ungroup()
combined_da_fil <- combined_da_fil %>%
  left_join(dplyr::select(compartments_filtered_max, gene_symbol, score), by = "gene_symbol")
blood_logFC_info <- blood_da_pr %>%
  dplyr::select(gene_symbol, logFC) %>%
  group_by(gene_symbol) %>%
  summarise(
    blood_direction = case_when(
      all(logFC > 0, na.rm = TRUE) ~ "Up",
      all(logFC < 0, na.rm = TRUE) ~ "Down",
      TRUE ~ "Mixed"
    ),
    blood_logFC = logFC[which.max(abs(logFC))],  # Select the highest absolute logFC
    .groups = "drop"
  )

# Join with combined_da_fil
combined_da_fil <- combined_da_fil %>%
  left_join(blood_logFC_info, by = "gene_symbol")

## Converting names for better readability in the figure
convert_column_names <- function(colname) {
  case_when(
    grepl("muscle-transcript", colname) ~ "Muscle - Transcript",
    grepl("muscle-prot-pr", colname) ~ "Muscle - Protein",
    grepl("muscle-prot-ph", colname) ~ "Muscle - Phospho",
    grepl("adipose-transcript", colname) ~ "Adipose - Transcript",
    grepl("adipose-prot-pr", colname) ~ "Adipose - Protein",
    grepl("adipose-prot-ph", colname) ~ "Adipose - Phospho",
    grepl("blood-transcript", colname) ~ "Blood - Transcript",
    TRUE ~ colname  # Default: Keep as is if unmatched
  )
}



################### Figure 8A and S8A ################
assign_heatmap_label <- function(tissue_logFC, blood_direction) {
  case_when(
    tissue_logFC > 0 & blood_direction == "Up"   ~ "Tissue Up & Blood Up",
    tissue_logFC < 0 & blood_direction == "Up"   ~ "Tissue Down & Blood Up",
    tissue_logFC > 0 & blood_direction == "Down" ~ "Tissue Up & Blood Down",
    tissue_logFC < 0 & blood_direction == "Down" ~ "Tissue Down & Blood Down",
    TRUE ~ "Mixed"  # If inconsistent cases exist
  )
}

heatmap_data_ext <- combined_da_fil %>%
  filter(score >= 4, blood_direction == "Up") %>%
  distinct(gene_symbol, assay, tissue, Exercise2, score, logFC, blood_direction) %>%
  mutate(
    column_name = paste(tissue, assay, sep = "-"),
    heatmap_label = assign_heatmap_label(logFC, blood_direction)
  ) %>%
  select(gene_symbol, column_name, heatmap_label) %>%
  pivot_wider(
    names_from = column_name,
    values_from = heatmap_label,
    values_fill = "Absent",
    values_fn = list(heatmap_label = ~ if (length(unique(.x)) > 1) "Mixed" else unique(.x))
  ) %>%
  column_to_rownames(var = "gene_symbol")

heatmap_matrix <- as.matrix(heatmap_data_ext)

# ✅ Define the color mapping **at plotting stage**
custom_heatmap_colors <- c(
  "Tissue Up & Blood Up"   = "#C84C05",
  "Tissue Down & Blood Up" = "#FFB200",
  "Mixed"                  = "purple",
  "Absent"                 = "white"
)

# ✅ Define row annotation for Exercise and Score
gene_metadata <- combined_da_fil %>%
  filter(score >= 4, blood_direction == "Up") %>%
  distinct(gene_symbol, Exercise2, score) %>%
  column_to_rownames(var = "gene_symbol")

score_colors <- colorRamp2(
  c(3.5, 5),  # ✅ Match to score range
  c("white", "green")
)

row_anno <- rowAnnotation(
  Exercise = gene_metadata$Exercise2,
  Score = anno_barplot(
    gene_metadata$score,
    gp = gpar(fill = score_colors(gene_metadata$score)),
    ylim = c(3.5, 5),
    baseline = 3.5,  # ← start from minimum shown value
    border = FALSE,  # ✅ no need to use NA here
    axis_param = list(
      at = c(3.5, 4, 4.5, 5),
      labels = c("3.5", "4", "4.5", "5")
    )
  ),
  col = list(Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77", "Both" = "black"))
)

convert_column_names <- function(colname) {
  case_when(
    grepl("muscle-transcript", colname) ~ "Muscle - Transcript",
    grepl("muscle-prot-pr", colname) ~ "Muscle - Protein",
    grepl("muscle-prot-ph", colname) ~ "Muscle - Phospho",
    grepl("adipose-transcript", colname) ~ "Adipose - Transcript",
    grepl("adipose-prot-pr", colname) ~ "Adipose - Protein",
    grepl("adipose-prot-ph", colname) ~ "Adipose - Phospho",
    grepl("blood-transcript", colname) ~ "Blood - Transcript",
    TRUE ~ colname  # Default: Keep as is if unmatched
  )
}

# Apply function to rename column names
colnames(heatmap_matrix) <- sapply(colnames(heatmap_matrix), convert_column_names)

# Generate heatmap
Heatmap(
  heatmap_matrix,
  name = "LogFC-Blood Direction",
  col = custom_heatmap_colors,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Blood & Tissue Direction",
    at = names(custom_heatmap_colors),  # Keep categorical values
    labels = c("Blood Up + Tissue Up", "Blood Up + Tissue Down",
               "Blood Up + Tissue Up/Down", "Absent"),
    legend_gp = gpar(fill = custom_heatmap_colors)  # Keep legend unchanged
  ),
  right_annotation = row_anno,  # Ensure modified Score annotation is applied
  border = TRUE
)


################# Figure S8 ##############
main_features <- rownames(heatmap_matrix)
heatmap_data_ext_s <- combined_da_fil %>%
  filter(!gene_symbol %in% main_features) %>% # Use >4 for the main figure. Also filtering for only blood up
  distinct(gene_symbol, assay, tissue, Exercise2, score, logFC, blood_direction) %>%
  mutate(
    column_name = paste(tissue, assay, sep = "-"),
    heatmap_label = assign_heatmap_label(logFC, blood_direction)  # Assign descriptive labels
  ) %>%
  select(gene_symbol, column_name, heatmap_label) %>%
  pivot_wider(
    names_from = column_name,
    values_from = heatmap_label,
    values_fill = "Absent",  # Default missing values to "Absent"
    values_fn = list(heatmap_label = ~ if (length(unique(.x)) > 1) "Mixed" else unique(.x))
  ) %>%
  column_to_rownames(var = "gene_symbol")
heatmap_matrix_s <- as.matrix(heatmap_data_ext_s)

custom_heatmap_colors_s <- c(
  "Tissue Up & Blood Up"   = "#C84C05",  # Red
  "Tissue Down & Blood Up" = "#FFB200",
  "Tissue Up & Blood Down" = "#DDEB9D",
  "Tissue Down & Blood Down" = "#638C6D",
  "Mixed"                  = "purple",
  "Absent"                 = "white"
)

# define row annotation for Exercise and Score
gene_metadata_s <- combined_da_fil %>%
  filter(!gene_symbol %in% main_features) %>%
  distinct(gene_symbol, Exercise2, score) %>%
  column_to_rownames(var = "gene_symbol")

score_colors_s <- colorRamp2(
  c(0, 5),  # Match to score range
  c("white", "green")
)

row_anno_s <- rowAnnotation(
  Exercise = gene_metadata_s$Exercise2,
  Score = anno_barplot(
    gene_metadata_s$score,
    gp = gpar(fill = score_colors_s(gene_metadata_s$score)),
    ylim = c(0, 5)  # Removed the incorrect comma
  ),  # Closed `anno_barplot()`
  col = list(Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77", "Both" = "black"))
)  # Closed `rowAnnotation()`


# Apply function to rename column names
colnames(heatmap_matrix_s) <- sapply(colnames(heatmap_matrix_s), convert_column_names)

Heatmap(
  heatmap_matrix_s,
  name = "LogFC-Blood Direction",
  col = custom_heatmap_colors_s,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Blood & Tissue Direction",
    at = names(custom_heatmap_colors_s),  # ✅ Keep categorical values
    labels = c("Blood Up + Tissue Up", "Blood Up + Tissue Down",
               "Blood Down + Tissue Up", "Blood Down + Tissue Down",
               "Mixed", "Absent"),
    legend_gp = gpar(fill = custom_heatmap_colors_s)  # ✅ Keep legend unchanged
  ),
  right_annotation = row_anno_s,  # ✅ Ensure modified Score annotation is applied
  border = TRUE
)



#-------------------------------------------------------------------------------
# Figure 8D and 8E

####################  Adipose CCN1 intra-correlation ###########################
adi_trans_qc <- motrpac_qc$adipose$`transcript-rna-seq`$qc_norm
adi_trans_meta <- motrpac_qc$adipose$`transcript-rna-seq`$sample_metadata #meta

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

timepoint_group_results <- list()

# Loop through each time point
for (tp in unique(adi_trans_meta_filtered$Timepoint)) {

  # If timepoint is "pre_exercise", compute correlation for all samples together
  if (tp == "pre_exercise") {
    tp_samples <- adi_trans_meta_filtered %>%
      filter(Timepoint == tp) %>%
      pull(vialLabel) %>%
      as.character()

    # Subset the gene expression matrix
    tp_data <- adi_trans_qc[, tp_samples, drop = FALSE]

    # Extract CCN1 expression
    ccn1_expression <- as.numeric(tp_data["ENSG00000142871.18", , drop = TRUE])  # Ensure numeric

    # Transpose tp_data to match dimensions
    tp_data_transposed <- t(tp_data)

    # Use bicorAndPvalue to compute bicorrelation and p-values
    tp_results <- bicorAndPvalue(tp_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

    # Store bicorrelation and p-values in a list
    timepoint_group_results[[tp]] <- list(
      bicor = tp_results$bicor,
      pvalue = tp_results$p
    )

  } else {
    # For other timepoints, compute correlation separately for each group
    for (grp in unique(adi_trans_meta_filtered$randomGroupCode)) {

      # Subset metadata for the current timepoint and group
      tp_grp_samples <- adi_trans_meta_filtered %>%
        filter(Timepoint == tp, randomGroupCode == grp) %>%
        pull(vialLabel) %>%
        as.character()

      # Subset the gene expression matrix
      tp_grp_data <- adi_trans_qc[, tp_grp_samples, drop = FALSE]

      # Skip if there are not enough samples
      if (ncol(tp_grp_data) < 3) next

      # Extract CCN1 expression
      ccn1_expression <- as.numeric(tp_grp_data["ENSG00000142871.18", , drop = TRUE])

      # Transpose tp_grp_data to match dimensions
      tp_grp_data_transposed <- t(tp_grp_data)

      # Compute bicorrelation and p-values
      tp_grp_results <- bicorAndPvalue(tp_grp_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

      # Store results using timepoint-group format
      tp_grp_name <- paste(tp, grp, sep = "-")
      timepoint_group_results[[tp_grp_name]] <- list(
        bicor = tp_grp_results$bicor,
        pvalue = tp_grp_results$p
      )
    }
  }
}

# Convert results into matrices
ccn1_adi_bicor <- do.call(cbind, lapply(timepoint_group_results, function(res) res$bicor))
ccn1_adi_pval <- do.call(cbind, lapply(timepoint_group_results, function(res) res$p))

# Assign column names based on timepoint-group
colnames(ccn1_adi_bicor) <- names(timepoint_group_results)
colnames(ccn1_adi_pval) <- names(timepoint_group_results)

# Check results
head(ccn1_adi_bicor)
head(ccn1_adi_pval)

##### For Figure 8D Adipose bicor heatmap, follow below.
ccn1_adi_hm <- ccn1_adi_bicor[, c("pre_exercise",
                                  "post_15_30_45_min-ADUEndur", "post_3.5_4_hr-ADUEndur", "post_24_hr-ADUEndur",
                                  "post_15_30_45_min-ADUResist", "post_3.5_4_hr-ADUResist", "post_24_hr-ADUResist")]

colnames_ccn1 <- colnames(ccn1_adi_hm)

# 2. Parse annotations
tissue_colors <- c("adipose" = "#ffffbf", "muscle" = "#abd9e9")
exercise_colors <- c("Pre" = "#bebebe", "EE" = "#d95f02", "RE" = "#1b9e77")
timepoint_colors <- c(
  "Pre" = "#bebebe",
  "30minPost" = "#AE76A3",
  "4hrPost" = "#882E72",
  "24hrPost" = "#61194F"
)
exercise <- ifelse(grepl("ADUEndur", colnames_ccn1), "EE",
                   ifelse(grepl("ADUResist", colnames_ccn1), "RE", "Pre"))

timepoint <- ifelse(grepl("30|15|45", colnames_ccn1), "30minPost",
                    ifelse(grepl("3.5_4_hr", colnames_ccn1), "4hrPost",
                           ifelse(grepl("24_hr", colnames_ccn1), "24hrPost", "Pre")))

# 3. Convert to factors (optional but useful)
exercise <- factor(exercise, levels = c("Pre", "EE", "RE"))
timepoint <- factor(timepoint, levels = c("Pre", "30minPost", "4hrPost", "24hrPost"))

# 4. Create column annotations
col_anno <- columnAnnotation(
  Exercise = exercise,
  Timepoint = timepoint,
  col = list(
    Exercise = exercise_colors,
    Timepoint = timepoint_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Exercise"),
    Timepoint = list(title = "Timepoint")
  )
)

# 5. Create heatmap
ccn1_adi_bicor_hm <- Heatmap(
  ccn1_adi_hm,
  name = "bicor",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = col_anno,
  column_split = exercise,
  column_title = NULL,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "Bicor")
) # FIgure 8D-adipose heatmap


##### For Figure 8E CAMERA plot, follow below.
# prepare genesets
genesets_list <- c(
  MOLECULAR_SIGNATURES$GOBP,
  MOLECULAR_SIGNATURES$GOCC,
  MOLECULAR_SIGNATURES$GOMF,
  MOLECULAR_SIGNATURES$BIOCARTA,
  MOLECULAR_SIGNATURES$PID,
  MOLECULAR_SIGNATURES$REACTOME,
  MOLECULAR_SIGNATURES$WP
)
genesets_list <- c(
  MOLECULAR_SIGNATURES$GOBP
)
mapping <- HUMAN_FEATURE_TO_GENE[
  HUMAN_FEATURE_TO_GENE$assay == "transcript-rna-seq",
  c("feature_id", "gene_symbol", "ensembl_gene")
]
# Keep only rows that are in mapping
keep_rows <- rownames(ccn1_adi_bicor) %in% mapping$feature_id
ccn1_adi_bicor <- ccn1_adi_bicor[keep_rows, ]

# Add gene symbols as new rownames
rownames(ccn1_adi_bicor) <- mapping$gene_symbol[
  match(rownames(ccn1_adi_bicor), mapping$feature_id)
]
ccn1_adi_z <- as.data.frame(apply(ccn1_adi_bicor, 2, function(x) {
  z <- 0.5 * log((1 + x) / (1 - x))
  z[!is.finite(z)] <- NA  # remove +/- Inf
  return(z)
}))

res <- cameraPR.matrix(
  statistic = as.matrix(ccn1_adi_z),
  index = genesets_list,
  use.ranks = FALSE,        # Use numeric values (not ranks)
  inter.gene.cor = 0.01,    # Default correlation adjustment
  alternative = "two.sided",
  min.size = 10             # Ignore very small sets
)
res <- res %>%
  mutate(tissue = "adi")

####################  Muscle CCN1 intra-correlation ###########################
mus_trans_qc <- motrpac_qc$muscle$`transcript-rna-seq`$qc_norm
mus_trans_meta <- motrpac_qc$muscle$`transcript-rna-seq`$sample_metadata #meta

mus_t <- mus_trans_meta %>%
  filter(visitcode == "ADU_BAS") %>%
  select(vialLabel)
mus_t_id <- mus_t$vialLabel
mus_trans_qc <- mus_trans_qc[, colnames(mus_trans_qc) %in% mus_t_id]

mus_trans_meta_filtered <- mus_trans_meta %>%
  filter(vialLabel %in% colnames(mus_trans_qc))

timepoint_group_results_mus <- list()

# Loop through each time point
for (tp in unique(mus_trans_meta_filtered$Timepoint)) {

  # If timepoint is "pre_exercise", compute correlation for all samples together
  if (tp == "pre_exercise") {
    tp_samples <- mus_trans_meta_filtered %>%
      filter(Timepoint == tp) %>%
      pull(vialLabel) %>%
      as.character()

    # Subset the gene expression matrix
    tp_data <- mus_trans_qc[, tp_samples, drop = FALSE]

    # Extract CCN1 expression
    ccn1_expression <- as.numeric(tp_data["ENSG00000142871.18", , drop = TRUE])  # Ensure numeric

    # Transpose tp_data to match dimensions
    tp_data_transposed <- t(tp_data)

    # Use bicorAndPvalue to compute bicorrelation and p-values
    tp_results <- bicorAndPvalue(tp_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

    # Store bicorrelation and p-values in a list
    timepoint_group_results_mus[[tp]] <- list(
      bicor = tp_results$bicor,
      pvalue = tp_results$p
    )

  } else {
    # For other timepoints, compute correlation separately for each group
    for (grp in unique(mus_trans_meta_filtered$randomGroupCode)) {

      # Subset metadata for the current timepoint and group
      tp_grp_samples <- mus_trans_meta_filtered %>%
        filter(Timepoint == tp, randomGroupCode == grp) %>%
        pull(vialLabel) %>%
        as.character()

      # Subset the gene expression matrix
      tp_grp_data <- mus_trans_qc[, tp_grp_samples, drop = FALSE]

      # Skip if there are not enough samples
      if (ncol(tp_grp_data) < 3) next

      # Extract CCN1 expression
      ccn1_expression <- as.numeric(tp_grp_data["ENSG00000142871.18", , drop = TRUE])

      # Transpose tp_grp_data to match dimensions
      tp_grp_data_transposed <- t(tp_grp_data)

      # Compute bicorrelation and p-values
      tp_grp_results <- bicorAndPvalue(tp_grp_data_transposed, ccn1_expression, use = "pairwise.complete.obs")

      # Store results using timepoint-group format
      tp_grp_name <- paste(tp, grp, sep = "-")
      timepoint_group_results_mus[[tp_grp_name]] <- list(
        bicor = tp_grp_results$bicor,
        pvalue = tp_grp_results$p
      )
    }
  }
}

# Convert results into matrices
ccn1_mus_bicor <- do.call(cbind, lapply(timepoint_group_results_mus, function(res) res$bicor))
ccn1_mus_pval <- do.call(cbind, lapply(timepoint_group_results_mus, function(res) res$p))

# Assign column names based on timepoint-group
colnames(ccn1_mus_bicor) <- names(timepoint_group_results_mus)
colnames(ccn1_mus_pval) <- names(timepoint_group_results_mus)

##### For Figure 7D Muscle bicor heatmap, follow below.
ccn1_mus_hm <- ccn1_mus_bicor[, c("pre_exercise",
                                  "post_15_30_45_min-ADUEndur", "post_3.5_4_hr-ADUEndur", "post_24_hr-ADUEndur",
                                  "post_15_30_45_min-ADUResist", "post_3.5_4_hr-ADUResist", "post_24_hr-ADUResist")]

colnames_ccn1 <- colnames(ccn1_mus_hm)

# 2. Parse annotations
exercise <- ifelse(grepl("ADUEndur", colnames_ccn1), "EE",
                   ifelse(grepl("ADUResist", colnames_ccn1), "RE", "Pre"))

timepoint <- ifelse(grepl("30|15|45", colnames_ccn1), "30minPost",
                    ifelse(grepl("3.5_4_hr", colnames_ccn1), "4hrPost",
                           ifelse(grepl("24_hr", colnames_ccn1), "24hrPost", "Pre")))

# 3. Convert to factors (optional but useful)
exercise <- factor(exercise, levels = c("Pre", "EE", "RE"))
timepoint <- factor(timepoint, levels = c("Pre", "30minPost", "4hrPost", "24hrPost"))

# 4. Create column annotations
col_anno <- columnAnnotation(
  Exercise = exercise,
  Timepoint = timepoint,
  col = list(
    Exercise = exercise_colors,
    Timepoint = timepoint_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Exercise"),
    Timepoint = list(title = "Timepoint")
  )
)

# 5. Create heatmap
ccn1_mus_bicor_hm <- Heatmap(
  ccn1_mus_hm,
  name = "bicor",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = col_anno,
  column_split = exercise,
  column_title = NULL,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "Bicor")
)

ccn1_adi_bicor_hm@column_title <- "CCN1 Adipose Intracorrelation"
ccn1_mus_bicor_hm@column_title <- "CCN1 Muscle Intracorrelation"

#note: this may take a while because this is a heatmap of all features.
draw(
  ccn1_adi_bicor_hm %v% ccn1_mus_bicor_hm,
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

### Run CAMERA-PR on muscle genes for Figure 8E
# Keep only rows that are in mapping
keep_rows2 <- rownames(ccn1_mus_bicor) %in% mapping$feature_id
ccn1_mus_bicor <- ccn1_mus_bicor[keep_rows2, ]

# Add gene symbols as new rownames
rownames(ccn1_mus_bicor) <- mapping$gene_symbol[
  match(rownames(ccn1_mus_bicor), mapping$feature_id)
]
ccn1_mus_z <- as.data.frame(apply(ccn1_mus_bicor, 2, function(x) {
  z <- 0.5 * log((1 + x) / (1 - x))
  z[!is.finite(z)] <- NA  # remove +/- Inf
  return(z)
}))

res_mus <- cameraPR.matrix(
  statistic = as.matrix(ccn1_mus_z),
  index = genesets_list,
  use.ranks = FALSE,        # Use numeric values (not ranks)
  inter.gene.cor = 0.01,    # Default correlation adjustment
  alternative = "two.sided",
  min.size = 10             # Ignore very small sets
)
res_mus <- res_mus %>%
  mutate(tissue = "mus")

######## CCN1 intracorrelation
cam_res_ccn1 <- rbind(res, res_mus)
cam_res_ccn1$Contrast <- factor(cam_res_ccn1$Contrast, levels =c(
  "pre_exercise", "post_15_30_45_min-ADUEndur", "post_3.5_4_hr-ADUEndur", "post_24_hr-ADUEndur",
  "post_15_30_45_min-ADUResist", "post_3.5_4_hr-ADUResist", "post_24_hr-ADUResist"
))

# remove any unused levels
cam_res_ccn1 <- droplevels(cam_res_ccn1)
cam_res_ccn1 <- cam_res_ccn1 %>%
  filter(!is.na(Contrast)) %>%  # Remove NA contrast values
  mutate(tissue_contrast = paste(tissue, Contrast, sep = "_"))

top_upregulated <- cam_res_ccn1 %>%
  group_by(GeneSet) %>%
  filter(all(c("adi", "mus") %in% unique(tissue))) %>%
  ungroup() %>%
  group_by(tissue_contrast) %>%
  slice_max(order_by = ZScore, n = 1, with_ties = FALSE) %>%  # Select 2 most upregulated pathways
  ungroup()

top_downregulated <- cam_res_ccn1 %>%
  group_by(GeneSet) %>%
  filter(all(c("adi", "mus") %in% unique(tissue))) %>%
  ungroup() %>%
  group_by(tissue_contrast) %>%
  slice_min(order_by = ZScore, n = 1, with_ties = FALSE) %>%
  ungroup()

top_terms <- bind_rows(top_upregulated, top_downregulated) %>%
  distinct(GeneSet) %>%  # eep only unique pathway names
  pull(GeneSet)  # Extract as a vector
# Define Column Order for Heatmap
heatmap_columns <- c(
  "adi_pre_exercise", "adi_post_15_30_45_min-ADUEndur", "adi_post_3.5_4_hr-ADUEndur", "adi_post_24_hr-ADUEndur",
  "adi_post_15_30_45_min-ADUResist", "adi_post_3.5_4_hr-ADUResist", "adi_post_24_hr-ADUResist",
  "mus_pre_exercise", "mus_post_15_30_45_min-ADUEndur", "mus_post_3.5_4_hr-ADUEndur", "mus_post_24_hr-ADUEndur",
  "mus_post_15_30_45_min-ADUResist", "mus_post_3.5_4_hr-ADUResist", "mus_post_24_hr-ADUResist"
)
tissue <- ifelse(grepl("^adi", heatmap_columns), "Adipose", "Muscle")

exercise <- ifelse(grepl("pre_exercise", heatmap_columns), "Pre",
                   ifelse(grepl("ADUEndur", heatmap_columns), "Endurance",
                          ifelse(grepl("ADUResist", heatmap_columns), "Resistance", NA)))

timepoint <- ifelse(grepl("pre_exercise", heatmap_columns), "Pre",
                    ifelse(grepl("post_15_30_45_min", heatmap_columns), "Post_30min",
                           ifelse(grepl("post_3.5_4_hr", heatmap_columns), "Post_4hr",
                                  ifelse(grepl("post_24_hr", heatmap_columns), "Post_24hr", NA))))

# Colors
tissue_colors <- c("Adipose" = "#ffffbf", "Muscle" = "#abd9e9")
exercise_colors <- c("Pre" = "#bebebe", "Endurance" = "#d95f02", "Resistance" = "#1b9e77")
timepoint_colors <- c("Pre" = "#bebebe", "Post_30min" = "#AE76A3", "Post_4hr" = "#882E72", "Post_24hr" = "#61194F")

# Build annotation
col_anno <- columnAnnotation(
  Tissue = factor(tissue, levels = c("Adipose", "Muscle")),
  Exercise = factor(exercise, levels = c("Pre", "Endurance", "Resistance")),
  Timepoint = factor(timepoint, levels = c("Pre", "Post_30min", "Post_4hr", "Post_24hr")),
  col = list(
    Tissue = tissue_colors,
    Exercise = exercise_colors,
    Timepoint = timepoint_colors
  )
)


# Generate Heatmap
#pdf("CCN1_camera_heatmap.pdf", width = 18, height = 18)
cam_res_ccn1 %>%
  filter(GeneSet %in% top_terms) %>%
  enrichmap(n_top = Inf,
            set_column = "GeneSet",
            statistic_column = "ZScore",
            contrast_column = "tissue_contrast",
            padj_column = "FDR",
            padj_legend_title = "BH Adjusted\nP-Value",
            heatmap_args = list(
              heatmap_legend_param = list(title = "Z-Score"),
              column_split = factor(sub("_.*", "", heatmap_columns), levels = c("adi", "mus")),
              top_annotation = col_anno,
              column_order = heatmap_columns,
              show_column_names = TRUE
            ),
            height = 18,
            width = 18) ### Somehow the column annotation colors don't match with the column names. The colors were subsequently manually altered in the Illustrator.
#dev.off()

### Figure S8A. Upset plots
# Assay to Omics Layer Mapping
assay_map <- c(
  "transcript-rna-seq" = "Transcript",
  "prot-pr" = "Protein",
  "prot-ph" = "Phospho"
)

# Assay Priority
assay_priority <- c("Protein" = 1, "Transcript" = 2, "Phospho" = 3)

# Prepare Input and Apply Assay Hierarchy
upset_input <- combined_da_fil %>%
  mutate(
    Omics = recode(assay, !!!assay_map),
    Tissue = case_when(
      grepl("blood", tissue, ignore.case = TRUE) ~ "Blood",
      grepl("muscle", tissue, ignore.case = TRUE) ~ "Muscle",
      grepl("adipose", tissue, ignore.case = TRUE) ~ "Adipose",
      TRUE ~ tissue
    )
  ) %>%
  select(gene_symbol, Tissue, Omics, Exercise2) %>%
  distinct() %>%
  mutate(priority = assay_priority[Omics]) %>%
  group_by(gene_symbol, Exercise2, Tissue) %>%
  slice_min(order_by = priority, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = Tissue, values_from = value, values_fill = FALSE) %>%
  mutate(
    Omics = factor(Omics, levels = c("Transcript", "Protein", "Phospho")),
    Exercise2 = factor(Exercise2, levels = c("EE", "RE", "Both"))
  )
# Set Exercise Levels
exercise_levels <- levels(upset_input$Exercise2)

# Create a named list to store each plot
upset_plots <- list()

# Loop and save each plot in the list
for (exercise in exercise_levels) {
  df_subset <- upset_input %>% filter(Exercise2 == exercise)

  p <- ComplexUpset::upset(
    data = df_subset,
    intersect = c("Muscle", "Adipose", "Blood"),
    min_size = 1,
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        mapping = aes(fill = Omics)
      ) +
        scale_fill_manual(
          values = c("Transcript" = "#4477AA", "Protein" = "#228833", "Phospho" = "#F3A02B"),
          guide = guide_legend(title = "Omics Layer")
        )
    ),
    set_sizes = upset_set_size() +
      theme(axis.text.y = element_text(size = 12, face = "bold")),
    width_ratio = 0.2,
    stripes = c("#ffffbf","#d7191c", "#abd9e9"),
    name = "Tissue Overlap"
  ) + ggtitle(paste("Exercise:", exercise))

  upset_plots[[exercise]] <- p
}

# Combine into one figure (3 rows, 1 column layout)
wrap_plots(upset_plots, ncol = 1) ## Note that the color stripes don't always match with tissue. These were manually fixed in the Illustrator.

ggsave("Figure_S8A_upset.pdf",
       plot = wrap_plots(upset_plots, ncol = 1),
       width = 8,
       height = 12,
       units = "in")
#### Figure S8C
####################### Adipose CCN1 vs Muslce CCN1 ############################
plot_correlation_by_timepoint <- function(mus_trans_qc, mus_trans_meta_filtered, adi_trans_qc, adi_trans_meta_filtered, gene_id, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Define group colors
  group_colors <- c(
    "RE" = "#1b9e77",
    "EE" = "#d95f02",
    "CON" = "#7570b3"
  )

  # Get unique time points
  timepoints <- unique(mus_trans_meta_filtered$Timepoint)

  # Loop through each time point
  for (tp in timepoints) {
    # Filter metadata for the current time point
    mus_tp_meta <- mus_trans_meta_filtered %>% filter(Timepoint == tp)
    adi_tp_meta <- adi_trans_meta_filtered %>% filter(Timepoint == tp)

    # Extract matching sample IDs
    mus_samples <- mus_tp_meta %>% pull(vialLabel) %>% as.character()
    adi_samples <- adi_tp_meta %>% pull(vialLabel) %>% as.character()

    # Subset expression data
    mus_tp_data <- mus_trans_qc[, mus_samples, drop = FALSE]
    adi_tp_data <- adi_trans_qc[, adi_samples, drop = FALSE]

    # Rename columns from vialLabel to pid
    colnames(mus_tp_data) <- mus_tp_meta %>% arrange(match(vialLabel, colnames(mus_tp_data))) %>% pull(pid)
    colnames(adi_tp_data) <- adi_tp_meta %>% arrange(match(vialLabel, colnames(adi_tp_data))) %>% pull(pid)

    # Get common samples
    common_samples <- intersect(colnames(mus_tp_data), colnames(adi_tp_data))

    # Subset to common samples
    mus_tp_data <- mus_tp_data[, common_samples]
    adi_tp_data <- adi_tp_data[, common_samples]

    # Extract gene expression
    muscle_gene_expr <- as.numeric(mus_tp_data[gene_id, ])
    adipose_gene_expr <- as.numeric(adi_tp_data[gene_id, ])

    # Create dataframe for plotting
    plot_data <- mus_tp_meta %>%
      filter(pid %in% common_samples) %>%
      mutate(
        Adipose = adipose_gene_expr,
        Muscle = muscle_gene_expr,

        # Recode group names
        randomGroupCode = recode(randomGroupCode,
                                 "ADUResist" = "RE",
                                 "ADUEndur" = "EE",
                                 "ADUControl" = "CON")
      )

    # Calculate correlation per group
    cor_results <- plot_data %>%
      group_by(randomGroupCode) %>%
      summarise(
        r = round(cor(Adipose, Muscle, use = "pairwise.complete.obs"), 3),
        p = signif(cor.test(Adipose, Muscle)$p.value, 3)
      ) %>%
      mutate(label = paste0(randomGroupCode, ": r = ", r, ", p = ", p)) %>%
      pull(label) %>%
      paste(collapse = "\n")

    # Generate scatter plot with group-specific trend lines
    p <- ggplot(plot_data, aes(x = Adipose, y = Muscle, color = randomGroupCode)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE) +  # Separate trend lines per group
      labs(
        title = paste("Correlation at", tp),
        x = "Adipose CCN1",
        y = "Muscle CCN1",
        caption = cor_results
      ) +
      scale_color_manual(name = "Group", values = group_colors) +
      theme_minimal(base_size = 15) +
      theme(plot.caption = element_text(hjust = 0.5, size = 12))

    # Save the plot
    plot_file <- file.path(output_dir, paste0(tp, "_correlation_scatter.pdf"))
    ggsave(plot_file, plot = p, width = 5, height = 5)
  }
}

# Example usage:
plot_correlation_by_timepoint(
  mus_trans_qc = mus_trans_qc,
  mus_trans_meta_filtered = mus_trans_meta_filtered,
  adi_trans_qc = adi_trans_qc,
  adi_trans_meta_filtered = adi_trans_meta_filtered,
  gene_id = "ENSG00000142871.18", # CCN1 gene
  output_dir = "correlation_plots"
)

#------------------------------------------------------------------------------
# Supplementary table. Includes features, tissue origin, directionality
combined_da_fil2 <- combined_da %>%
  filter(gene_symbol %in% blood_da_pr$gene_symbol) %>%
  left_join(
    blood_da_pr %>% dplyr::select(gene_symbol, plasma_contrast = contrast),
    by = "gene_symbol"
  )
sec_30min <- combined_da_fil2 %>%
  filter(str_detect(contrast, "30min")) %>%
  filter(str_detect(plasma_contrast, "30min|4hr|24hr")) %>%
  pull(gene_symbol) %>%
  unique() # 40

sec_4hr <- combined_da_fil2 %>%
  filter(str_detect(contrast, "4hr")) %>%
  filter(str_detect(plasma_contrast, "4hr|24hr")) %>%
  pull(gene_symbol) %>%
  unique() # 35

sec_24hr <- combined_da_fil2 %>%
  filter(str_detect(contrast, "24hr")) %>%
  filter(str_detect(plasma_contrast, "24hr")) %>%
  pull(gene_symbol) %>%
  unique() # 0

opposing_genes <- blood_da_pr %>%
  group_by(gene_symbol) %>%
  filter(n() > 1) %>%  # Keep only duplicated gene symbols
  summarise(has_opposing = any(logFC > 0) & any(logFC < 0)) %>%  # Check if both + and - logFC exist
  filter(has_opposing) %>%
  pull(gene_symbol)  # Extract only gene names with opposing logFC

# Print results
print(opposing_genes) #"CD300E" "GH1"    "SPOCK1" "STC1"  These are blood da features that have opposing logFC directions...


# Subset and rename columns
subset_da_fil <- combined_da_fil %>%
  select(
    gene_symbol, assay, tissue, contrast, Exercise2, logFC, score, blood_logFC
  ) %>%
  rename(
    Extracellular_score = score,
    Exercise = Exercise2,
    Tissue_logFC = logFC
  )

# Define mapping for contrast values
contrast_map <- c(
  "EE_30min_vs_Control"  = "EE_Post_15_30_45min",
  "RE_30min_vs_Control"  = "RE_Post_15_30_45min",
  "RE_4hr_vs_Control"    = "RE_Post_3.5_4hr",
  "EE_4hr_vs_Control"    = "EE_Post_3.5_4hr",
  "RE_24hr_vs_Control"   = "RE_Post_24hr",
  "EE_24hr_vs_Control"   = "EE_Post_24hr",
  "EE_d40min_vs_Control" = "EE_During_40min",
  "EE_d20min_vs_Control" = "EE_During_20min",
  "RE_10min_vs_Control"  = "RE_Post_10min",
  "EE_10min_vs_Control"  = "EE_Post_10min"
)

# Apply mapping
subset_da_pretty <- subset_da_fil %>%
  mutate(
    Tissue_direction = ifelse(Tissue_logFC > 0, "Up", "Down"),
    Plasma_direction = ifelse(blood_logFC > 0, "Up", "Down"),
    Directionality_concordant = ifelse(Tissue_direction == Plasma_direction, "Yes", "No"),
    Temporally_concordant = ifelse(gene_symbol %in% union(sec_30min, sec_4hr), "Yes", "No"),
    contrast = recode(contrast, !!!contrast_map)
  )
subset_da_pretty <- subset_da_pretty %>%
  rename(Modality = Exercise) %>%
  select(-Tissue_logFC, -blood_logFC) # Supplementary Table 8 ready

