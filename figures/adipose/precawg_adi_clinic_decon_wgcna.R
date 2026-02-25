# PreCWAG-Adipose: Clinical. These codes cover Figures 1C, S1A-B, S4A, 5B, 5C, S5A-C
# Load library
library(dplyr)
library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(readr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(TMSig)
library(purrr)
library(tibble)
library(circlize)
library(RColorBrewer)
library(WGCNA)
library(ComplexHeatmap)
library(reshape2)
library(lmerTest)
library(scales)

# qc data and DA results

#the `load_qc` function requires access to sample level data (via "MotrpacHumanPreSuspensionData"). To repeat this specific analysis, you will
#need to request access to sample level data via a request to the consortium.
motrpac_qc <- load_qc()
motrpac_da <- load_differential_analysis()

# extract pid.
adi_trans_meta <- motrpac_qc$adipose$`transcript-rna-seq`$sample_metadata #meta
adi_pre_t <- adi_trans_meta %>%
  filter(Timepoint == "pre_exercise",
         visitcode == "ADU_BAS") %>%
  dplyr::select(pid) %>%
  pull() # vecctor of 172 PIDs.
adi_pre_t <- as.character(adi_pre_t)
# Construct clinical data matrix
# Age, sex
cli_pheno <- pheno$data %>%
  filter(pid %in% adi_pre_t & visitcode == "ADU_BAS" & Timepoint == "pre_exercise") %>%
  select(pid, Age_groups, Sex) %>%
  dplyr::distinct(pid, .keep_all = TRUE)

# vitals: BMI, WC, SBP, DBP, HR
cli_vital <- cln_curated_anthropometrics_vitals$data %>%
  filter(pid %in% adi_pre_t & visit_code == "ADU_SCP") %>%
  select(pid, bmi, wccm, sys_bp, dias_bp, hr_avg)

# daily total steps, vector magnitude, energy expenditure
cli_accel <- cln_curated_accel_derived_variables_baseline$data %>%
  filter(pid %in% adi_pre_t) %>%
  select(pid, daily_total_steps, daily_total_vm, daily_total_eneg_exp)

# vo2peak, relative vo2peak, O2 pulse, grip strength, knee torque
cli_ex <- cln_curated_ex_performance_testing$data %>%
  filter(pid %in% adi_pre_t & visit_code == "ADU_SCP") %>%
  select(pid, peaktorq_iske, peak_handgrip, VO2max_L, VO2max, O2_pulse)

# HbA1c, total cholesterol, HDL, LDL, Triglyceride, glucose
cli_lab <- cln_raw_local_lab_results$data %>%
  filter(pid %in% adi_pre_t & visit_code == "ADU_SCP") %>%
  select(pid, hba1c_labr, totalc_labr, hdlc_labr, ldlc_labr, trig_labr, gluc_labr)

cli_pheno <- cli_pheno %>% mutate(pid = as.character(pid))
cli_vital <- cli_vital %>% mutate(pid = as.character(pid))
cli_accel <- cli_accel %>% mutate(pid = as.character(pid))
cli_ex    <- cli_ex %>% mutate(pid = as.character(pid))
cli_lab   <- cli_lab %>% mutate(pid = as.character(pid))
# merge all clinical tables
cli_all <- reduce(
  list(cli_pheno, cli_vital, cli_accel, cli_ex, cli_lab),
  full_join,
  by = "pid"
)

## Blood metab data
filter_metabolite_matrix <- function(matrix, metadata,
                                     timepoint = "pre_exercise",
                                     visitcode = "ADU_BAS",
                                     pid_filter = NULL) {
  # Filter metadata
  filtered_meta <- metadata %>%
    dplyr::filter(Timepoint == timepoint, visitcode == visitcode)

  # Vial labels to keep (must exist in matrix colnames)
  keep_vials <- intersect(filtered_meta$vialLabel, colnames(matrix))

  if (length(keep_vials) == 0) {
    stop("No matching vialLabels found in matrix for given filters.")
  }

  # Subset matrix
  filtered_matrix <- matrix[, keep_vials, drop = FALSE]

  # Map vialLabel → pid
  vial_to_pid <- filtered_meta %>%
    dplyr::filter(vialLabel %in% keep_vials) %>%
    dplyr::select(vialLabel, pid) %>%
    dplyr::distinct()

  # Reorder pids to match columns
  pid_ordered <- vial_to_pid$pid[match(colnames(filtered_matrix), vial_to_pid$vialLabel)]

  # Replace colnames with pid
  colnames(filtered_matrix) <- pid_ordered

  # Optional: filter by pid_filter
  if (!is.null(pid_filter)) {
    pid_filter <- as.character(pid_filter)  # ensure character
    keep_pids <- intersect(colnames(filtered_matrix), pid_filter)
    filtered_matrix <- filtered_matrix[, keep_pids, drop = FALSE]
  }

  return(filtered_matrix)
}

#update: Feb 25, 2026. The clinical analytes are also now available through `MotrpacHumanPreSuspensionData` (still only available via request).
#this implementation needs to be revamped.
# setwd("/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/6 - MoTrPAC Pre cawg/R/Files/clinical_metab/")
# 1. insulin
precawg_ins <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-imm-ins_qc-norm_log2_v1.2.txt",
  row.names = 1,
  check.names = FALSE
)
precawg_ins_meta <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-imm-ins_metadata_samples_v1.2.txt",
  check.names = FALSE
)
insulin_pre <- filter_metabolite_matrix(
  matrix = precawg_ins,
  metadata = precawg_ins_meta,
  pid_filter = adi_pre_t
)
# 2. t-conv
precawg_conv <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-conv_qc-norm_log2_v1.2.txt",
  row.names = 1,
  check.names = FALSE
)
precawg_conv_meta <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-conv_metadata_samples_v1.2.txt",
  check.names = FALSE
)
conv_pre <- filter_metabolite_matrix(
  matrix = precawg_conv,
  metadata = precawg_conv_meta,
  pid_filter = adi_pre_t
)
conv_pre <- conv_pre[!rownames(conv_pre) %in% "Glucose", ]

# 3. BCAAs
precawg_amine <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-amines_qc-norm_log2_v1.2.txt",
  row.names = 1,
  check.names = FALSE
)
precawg_amine_meta <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-amines_metadata_samples_v1.2.txt",
  check.names = FALSE
)
amine_pre <- filter_metabolite_matrix(
  matrix = precawg_amine,
  metadata = precawg_amine_meta,
  pid_filter = adi_pre_t
)
amine_pre <- amine_pre[c("Leucine", "Isoleucine", "Valine"),]
# 4. cortisol
precawg_cor <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-imm-crt_qc-norm_log2_v1.2.txt",
  row.names = 1,
  check.names = FALSE
)
precawg_cor_meta <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-imm-crt_metadata_samples_v1.2.txt",
  check.names = FALSE
)
cortisol_pre <- filter_metabolite_matrix(
  matrix = precawg_cor,
  metadata = precawg_cor_meta,
  pid_filter = adi_pre_t
)
# 5. glucagon
precawg_glc <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-imm-glc_qc-norm_log2_v1.2.txt",
  row.names = 1,
  check.names = FALSE
)
precawg_glc_meta <- read.delim(
  "human-precovid-sed-adu_t02-plasma_metab-t-imm-glc_metadata_samples_v1.2.txt",
  check.names = FALSE
)
glucagon_pre <- filter_metabolite_matrix(
  matrix = precawg_glc,
  metadata = precawg_glc_meta,
  pid_filter = adi_pre_t
)

metab_list <- list(
  Insulin   = insulin_pre,
  Conv      = conv_pre,
  Amines    = amine_pre,
  Cortisol  = cortisol_pre,
  Glucagon  = glucagon_pre
)

metab_dfs <- metab_list %>%
  imap(function(df, name) {
    df_t <- as.data.frame(t(df)) %>%        # transpose (samples as rows)
      rownames_to_column("pid") %>%         # move colnames → pid column
      mutate(pid = as.numeric(pid)) %>%     # ensure pid is numeric
      rename_with(~ paste(name, ., sep = "_"), -pid)  # prefix with list name
  })
metab_all <- reduce(metab_dfs, full_join, by = "pid")
cli_full <- full_join(
  cli_all %>% mutate(pid = as.character(pid)),
  metab_all %>% mutate(pid = as.character(pid)),
  by = "pid"
)

# renmae variables for publication purpose
cli_full <- cli_full %>%
  rename(
    Age = Age_groups,
    BMI = bmi,
    WC = wccm,
    SBP = sys_bp,
    DBP = dias_bp,
    HR = hr_avg,
    Glucose = gluc_labr,
    HbA1c = hba1c_labr,
    Cholesterol = totalc_labr,
    HDL = hdlc_labr,
    LDL = ldlc_labr,
    Grip_strength = peak_handgrip,
    VO2max_rel = VO2max,
    VO2max = VO2max_L,
    O2_Pulse = O2_pulse,
    Knee_torque = peaktorq_iske,
    Insulin = Insulin_Insulin,
    Trig = trig_labr,
    Lactate = Conv_Lactate,
    Glycerol = Conv_Glycerol,
    NEFA = Conv_NEFA,
    KET = Conv_KET,
    Leucine = Amines_Leucine,
    Isoleucine = Amines_Isoleucine,
    Valine = Amines_Valine,
    Cortisol = Cortisol_Cortisol,
    Glucagon = Glucagon_Glucagon,
    Steps = daily_total_steps,
    Movement = daily_total_vm,
    Energy_exp = daily_total_eneg_exp
  )

# Convert metabolite data back to original
metab_vars <- c(
  "Insulin", "Lactate", "Glycerol", "NEFA",
  "KET", "Glucagon", "Cortisol", "Isoleucine", "Leucine", "Valine"
)
cli_full <- cli_full %>%
  mutate(across(all_of(metab_vars), ~ 2^.x - 1))
# Convert insulin to uU/ml
cli_full <- cli_full %>%
  mutate(Insulin = Insulin * 0.023)
# Calculate HOMA_IR and Adipo_IR
cli_full <- cli_full %>%
  mutate(
    HOMA_IR = (Glucose * Insulin) / 405,
    Adipo_IR = NEFA * Insulin
  )

# Encode Sex and Age
cli_full <- cli_full %>%
  mutate(
    Sex = ifelse(Sex == "Female", 1, 0),
    Age = case_when(
      Age == "10-20" ~ 1,
      Age == "20-30" ~ 2,
      Age == "30-40" ~ 3,
      Age == "40-50" ~ 4,
      Age == "50-60" ~ 5,
      Age == "60-70" ~ 6,
      Age == "70-80" ~ 7,
      TRUE ~ NA_real_
    )
  )
cli_full <- cli_full %>%
  column_to_rownames("pid")
# Convert any Inf to NA
cli_full[] <- lapply(cli_full, function(col) {
  col[is.infinite(col)] <- NA
  col
})
# z scale
cli_full_z <- as.data.frame(scale(cli_full))


## Figure S1A: Intercorrelation of clinical variables
cc1 = bicorAndPvalue(cli_full_z, cli_full_z, use = 'p')
dim(cc1$p) <- dim(cc1$bicor)
dimnames(cc1$p) <- dimnames(cc1$bicor)

ht_colors <- brewer.pal(3, "RdBu")

cor_mat <- cc1$bicor
p_mat <- cc1$p
col_fun <- colorRamp2(
  c(1, 0, -1),
  ht_colors
)

column_groups_named <- c(
  "Sex" = "Anthropometrics",
  "Age" = "Anthropometrics",
  "BMI" = "Anthropometrics",
  "WC" = "Anthropometrics",
  "SBP" = "Anthropometrics",
  "DBP" = "Anthropometrics",
  "HR" = "Anthropometrics",
  "HbA1c" = "Blood",
  "Cholesterol" = "Blood",
  "HDL" = "Blood",
  "LDL" = "Blood",
  "Grip_strength" = "Muscle strength",
  "VO2max_rel" = "Cardiorespiration",
  "VO2max" = "Cardiorespiration",
  "O2_Pulse" = "Cardiorespiration",
  "Knee_torque" = "Muscle strength",
  "Steps" = "Accelerometry",
  "Movement" = "Accelerometry",
  "Energy_exp" = "Accelerometry",
  "Insulin" = "Blood",
  "Glucose" = "Blood",
  "Lactate" = "Blood",
  "Glycerol" = "Blood",
  "NEFA" = "Blood",
  "Trig" = "Blood",
  "KET" = "Blood",
  "Glucagon" = "Blood",
  "Cortisol" = "Blood",
  "Isoleucine" = "Plasma BCAA",
  "Leucine" = "Plasma BCAA",
  "Valine" = "Plasma BCAA",
  "HOMA_IR" = "IR_index",
  "Adipo_IR" = "IR_index"
)

group_colors <- c(
  "Anthropometrics" = "#FCC737",
  "Accelerometry"   = "#BDE8CA",
  "Cardiorespiration" = "#FF8000",
  "Muscle strength" = "#0A5EB0",
  "IR_index"        = "#41B3A2",
  "Blood"           = "#F95454",
  "Plasma BCAA"     = "#D7C3F1"
)

# Step 5: Map group names to row/col correctly
row_group_vec <- column_groups_named[rownames(cor_mat)]
col_group_vec <- column_groups_named[colnames(cor_mat)]

# Step 6: Build annotations
column_annotation <- HeatmapAnnotation(
  Clinical_traits = col_group_vec,
  col = list(Clinical_traits = group_colors),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

row_annotation <- rowAnnotation(
  Clinical_traits = row_group_vec,
  col = list(Clinical_traits = group_colors),
  show_annotation_name = FALSE
)
#pdf("clinical_intercorrelation.pdf", width = 11, height = 8)
Heatmap(
  cor_mat,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = column_annotation,
  right_annotation = row_annotation,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  na_col = "white",
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    pval <- p_mat[i, j]
    if (!is.na(pval)) {
      label <- if (pval <= 0.001) "***"
      else if (pval <= 0.01) "**"
      else if (pval <= 0.05) "*"
      else ""
      if (label != "") {
        grid.text(label, x, y, gp = gpar(fontsize = 9, fontface = "bold"))
      }
    }
  }
) # Figure S1A

#dev.off()

## Correlation of Deconvolution outcome and clinical variables
decon_pre <- read.csv("../precawg_adi_baseline_decon.csv", row.names = 1)

decon_pre_z <- decon_pre %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

# align clinical matrix and decon matrix to have same order of pid
common_samples <- intersect(rownames(cli_full_z), rownames(decon_pre_z))

# Subset both data frames to keep only the common samples and ensure they are in the same order
cli_full_z2 <- cli_full_z[common_samples, ]
decon_pre_z2 <- decon_pre_z[common_samples, ]
stopifnot(identical(rownames(cli_full_z2), rownames(decon_pre_z2))) # confirm the order


## rearrage the order of columns
new_order <- c(
  "Sex", "Age", "BMI", "WC", "SBP", "DBP", "HR",
  "Steps", "Movement", "Energy_exp",
  "VO2max_rel", "VO2max", "O2_Pulse",
  "Grip_strength", "Knee_torque",
  "HOMA_IR", "Adipo_IR",
  "HbA1c", "Insulin", "Glucose", "Lactate", "Glycerol", "NEFA", "Trig", "Cholesterol", "HDL", "LDL", "KET", "Glucagon", "Cortisol",
  "Isoleucine", "Leucine", "Valine"
)
cli_full_z2 <- cli_full_z2[, new_order]

cor_decon_clin = bicorAndPvalue(decon_pre_z2, cli_full_z2) # alternative biweight midcor2
df_decon_clin = melt(cor_decon_clin$bicor) %>% dplyr::rename(bicor = value)
df_decon_clin$pval = melt(cor_decon_clin$p)$value
df_decon_clin$obs = melt(cor_decon_clin$nObs)$value
df_decon_clin <- df_decon_clin %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

# set up correlation dataframe for heatmap
deconXclin = dcast(df_decon_clin, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(deconXclin) = deconXclin$Var1
deconXclin$Var1 = NULL

# set up signif overlay
sig_table_decon_clin = dcast(df_decon_clin, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_decon_clin) = sig_table_decon_clin$Var1
sig_table_decon_clin$Var1 = NULL

## Complex Heatmap
color_palette <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))
col_fun <- color_palette(10)
col_group_vec <- column_groups_named[colnames(deconXclin)]

# Build annotations
column_annotation <- HeatmapAnnotation(
  Clinical_traits = col_group_vec,
  col = list(Clinical_traits = group_colors),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

# Figure 1C
#pdf(file = 'deconXtrait_complexheat.pdf', width = 12, height = 3)
Heatmap(
  as.matrix(deconXclin),
  name = "Correlation",
  col = col_fun,
  top_annotation = column_annotation,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, col = "black"),
  column_names_side = "top",
  show_heatmap_legend = TRUE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    sig_marker <- sig_table_decon_clin[i, j]

    if (!is.na(sig_marker) && sig_marker != "") {
      # Compute luminance to choose text color
      fill_rgb <- col2rgb(fill)
      luminance <- (0.299 * fill_rgb[1] + 0.587 * fill_rgb[2] + 0.114 * fill_rgb[3]) / 255
      text_color <- if (luminance < 0.5) "white" else "black"

      grid.text(sig_marker, x, y, gp = gpar(fontsize = 10, col = text_color))
    }
  }
)
#dev.off()

### Additional deconvolution-related figures: Figure S1B and S4A
# S1B is baseline sex comparison. S4A is deconvolution in all acute exercise time points.
# Import expanded deconvolution outcome that includes data from all time points.
dtangle_props <- read.csv("../precawg_decon_BA_6000hvg.csv", row.names = 1)

#Adding meta data
adi_pre_t_vial <- adi_trans_meta %>%
  filter(Timepoint == "pre_exercise",
         visitcode == "ADU_BAS") %>%
  dplyr::select(vialLabel) %>%
  pull()

dtangle_props <- dtangle_props %>%
  rownames_to_column(var = "vialLabel")
adi_trans_meta <- adi_trans_meta %>%
  mutate(vialLabel = as.character(vialLabel))
adi_trans_meta_subset <- adi_trans_meta %>%
  dplyr::select(vialLabel, pid, Timepoint, Sex, BMI, randomGroupCode)
dtangle_props <- dtangle_props %>%
  left_join(adi_trans_meta_subset, by = "vialLabel")

#Remove NA rows
dtangle_props <- dtangle_props %>%
  filter(!is.na(pid))

#visualize
dtangle_long <- dtangle_props %>%
  pivot_longer(cols = Adip.1:Vascular, names_to = "CellType", values_to = "Proportion")
dtangle_long$Timepoint <- factor(dtangle_long$Timepoint, levels = c("pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"))
dtangle_filtered <- dtangle_long %>%
  filter(randomGroupCode %in% c("ADUEndur", "ADUResist"))

dtangle_filtered <- dtangle_filtered %>%
  mutate(
    # 1) Clean timepoint labels
    tp_label = case_when(
      Timepoint == "pre_exercise"          ~ "Pre",
      Timepoint == "post_15_30_45_min"     ~ "Post 30min",
      Timepoint == "post_3.5_4_hr"         ~ "Post 4hr",
      Timepoint == "post_24_hr"            ~ "Post 24hr",
      TRUE ~ NA_character_
    ),

    # 2) Clean group labels
    grp_label = case_when(
      randomGroupCode == "ADUEndur"  ~ "EE",
      randomGroupCode == "ADUResist" ~ "RE",
      TRUE ~ randomGroupCode
    ),

    # 3) Combine them
    group_time = factor(
      paste(grp_label, tp_label),
      levels = c(
        "EE Pre", "EE Post 30min", "EE Post 4hr", "EE Post 24hr",
        "RE Pre", "RE Post 30min", "RE Post 4hr", "RE Post 24hr"
      )
    )
  )
# Figure S4A
ggplot(dtangle_filtered, aes(x = group_time, y = Proportion, fill = randomGroupCode)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8),
           width = 0.7, color = "black") +  # Bar plot with mean and black outline
  stat_summary(fun = "mean", fun.min = function(x) mean(x), fun.max = function(x) mean(x) + sd(x),
               geom = "errorbar", position = position_dodge(width = 0.8), width = 0.3) +  # Error bars only on top
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.3, color = "black") +  # Overlay individual points
  facet_wrap(~ CellType, scales = "free_y", nrow = 1) +  # One row of plots
  scale_fill_manual(values = c("#d95f02", "#1b9e77"), labels = c("Endurance", "Resistance")) +  # Colors and legend labels
  labs(x = "", y = "Proportion (%)", fill = "Modality") +  # Change legend title here
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(size = 12),
    panel.grid = element_blank()  # Remove gridlines
  )



# Baseline sex comparison (Figure S1B)
pre_exercise_data <- subset(dtangle_props, Timepoint == "pre_exercise")
cell_types <- colnames(pre_exercise_data)[2:8]  # Assuming 'Adip.1' to 'Vascular' are in columns 2 to 8

# Function to perform t-test for each cell type
t_test_results <- lapply(cell_types, function(cell_type) {
  t.test(pre_exercise_data[[cell_type]] ~ pre_exercise_data$Sex, data = pre_exercise_data)
})

# Name the list elements with the corresponding cell types
names(t_test_results) <- cell_types

# Print t-test results
t_test_results

dtangle_long <- pre_exercise_data %>%
  pivot_longer(cols = Adip.1:Vascular, names_to = "Cell_Type", values_to = "Value")

# Figure S1B
ggplot(dtangle_long, aes(x = Sex, y = Value, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8),
           width = 0.7, color = "black") +
  stat_summary(fun = "mean", fun.min = function(x) mean(x), fun.max = function(x) mean(x) + sd(x),
               geom = "errorbar", position = position_dodge(width = 0.8), width = 0.3) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.3, color = "black") +
  facet_wrap(~ Cell_Type, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("#f95c6f", "#5555ff"), labels = c("Female", "Male")) +
  theme_minimal() +
  labs(x = "", y = "Proportion (%)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12),
    panel.grid = element_blank()
  )

#### WGCNA
# Build networks
# Transcriptomics
trans_qc <- motrpac_qc$adipose$`transcript-rna-seq`$qc_norm
trans_meta <- motrpac_qc$adipose$`transcript-rna-seq`$sample_metadata #meta
adi_pre_t <- trans_meta %>%
  filter(Timepoint == "pre_exercise",
         visitcode == "ADU_BAS") %>%
  dplyr::select(vialLabel)
adi_pre_t_id <- adi_pre_t$vialLabel

trans_qc_pre <- trans_qc[, colnames(trans_qc) %in% adi_pre_t_id]
trans_qc_pre <- trans_qc_pre %>%
  rownames_to_column(var = "feature_id") %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% select(feature_id, gene_symbol), by = "feature_id") %>%
  filter(!is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>%
  slice(which.max(rowSums(across(where(is.numeric))))) %>%
  ungroup() %>%
  select(-feature_id) %>%
  column_to_rownames(var = "gene_symbol") # This code converts to gene symbols and remove duplicates. Prioritize higher expression. As a result, gene count becomes 17287

vial_to_pid <- trans_meta %>%
  select(vialLabel, pid) %>%
  drop_na() %>%
  distinct() %>%
  deframe()
colnames(trans_qc_pre) <- vial_to_pid[colnames(trans_qc_pre)]

#### filter for ENTS detected in >25% of samples and prepare counts matrix for WGCNA ####
#filter for ENTS that are detected in over 25%
trans_qc_pre$Zero_count = rowSums(trans_qc_pre == 0)
hist(trans_qc_pre$Zero_count)
filtered_matrix_trans1 = trans_qc_pre[trans_qc_pre$Zero_count<0.25*172,] # filtered out low count genes
wgcna_input_trans <- filtered_matrix_trans1 %>%
  dplyr::select(-Zero_count)### 2.0 - Check for outlier/bad genes and samples ####

datExpr0 = data.frame(t(wgcna_input_trans))

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# check if all genes/samples are good
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK #[1] TRUE

# Soft thresholding and network topology analysis ####

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)


# Module Generation
pwr = 6
ngenes = 100

# takes several hours
net = blockwiseModules(datExpr0,
                       power = pwr, #determined above
                       maxBlockSize = 20000, # we can set this high because we have a lot of ram. We have ~17000 genes so this will do them all at onec
                       TOMType = "signed", # unsigned was default
                       minModuleSize = ngenes,
                       reassignThreshold = 0.3,
                       mergeCutHeight = 0.3,
                       numericLabels = TRUE, # names modules as numbers isntead of colors
                       pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       nThreads = 48,
                       verbose = 3)

MEs = net$MEs

# ME0 is error module, but not removing for now.
module_membership = as.data.frame(net$colors)
colnames(module_membership) = 'module'
module_membership$ID = row.names(module_membership)
module_membership$module <- paste0("T", module_membership$module)
module_membership$module <- factor(module_membership$module, levels = paste0("T", 0:12))
module_membership <- module_membership[order(module_membership$module), ]

## Proteomics
prot_qc <- motrpac_qc$adipose$`prot-pr`$qc_norm
prot_meta <- motrpac_qc$adipose$`prot-pr`$sample_metadata

gene_mapping <- HUMAN_FEATURE_TO_GENE %>% # load Humanpresuspension library
  ungroup() %>%  # Ungroup to prevent automatic grouping
  dplyr::select(uniprot, gene_symbol) %>%
  filter(!is.na(uniprot), !is.na(gene_symbol)) %>%
  distinct()  # Remove duplicates if any
#Prcoess prot df
colnames(prot_qc) <- sub("^X", "", colnames(prot_qc))
prot_qc_long <- prot_qc %>%
  rownames_to_column("uniprot") %>%
  pivot_longer(-uniprot, names_to = "vialLabel", values_to = "expression")
prot_meta <- prot_meta %>%
  mutate(vialLabel = as.character(vialLabel))
prot_qc_long <- prot_qc_long %>%
  left_join(prot_meta %>% dplyr::select(vialLabel, pid, timepointDescription), by = "vialLabel")

prot_pre <- prot_qc_long %>%
  filter(timepointDescription %in% c("Rest 1", "Pre ex", "")) %>%

  # Summarize duplicates by averaging expression values
  group_by(uniprot, pid) %>%
  summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%

  # Pivot to wide format with uniprot as rows and pid as columns
  pivot_wider(names_from = pid, values_from = expression) %>%

  # Set uniprot as row names
  column_to_rownames("uniprot")

#### WGCNA pipeline for proteomics
prot_pre$Zero_count <- rowSums(is.na(prot_pre))
hist(prot_pre$Zero_count)
prot_pre1 = prot_pre[prot_pre$Zero_count<0.25*22,] # filtered out low count genes. From 7758 to 7098 proteins
wgcna_input_prot <- prot_pre1 %>%
  dplyr::select(-Zero_count)### 2.0 - Check for outlier/bad genes and samples ####

datExpr1 = data.frame(t(wgcna_input_prot))

options(stringsAsFactors = FALSE)

# check if all genes/samples are good
gsg = goodSamplesGenes(datExpr1, verbose = 3)
gsg$allOK #[1] TRUE

# Soft thresholding and network topology analysis ####

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr1, powerVector = powers, verbose = 5)

# Module construction
pwr = 12
ngenes = 50

# takes several hours
net1 = blockwiseModules(datExpr1,
                        power = pwr, #determined above
                        maxBlockSize = 20000, # we can set this high because we have a lot of ram. We have ~17000 genes so this will do them all at onec
                        TOMType = "signed", # unsigned was default
                        minModuleSize = ngenes,
                        reassignThreshold = 0.3,
                        mergeCutHeight = 0.3,
                        numericLabels = TRUE, # names modules as numbers isntead of colors
                        pamRespectsDendro = FALSE,
                        saveTOMs = F,
                        nThreads = 48,
                        verbose = 3)

MEs1 = net1$MEs

module_membership1 = as.data.frame(net1$colors)
colnames(module_membership1) = 'module'
module_membership1$ID = row.names(module_membership1)
module_membership1$module <- paste0("Pr", module_membership1$module)
module_membership1$module <- factor(module_membership1$module, levels = paste0("Pr", 0:14))
module_membership1 <- module_membership1[order(module_membership1$module), ]

## Phosphoproteomics
phos_qc <- motrpac_qc$adipose$`prot-ph`$qc_norm
phos_meta <- motrpac_qc$adipose$`prot-ph`$sample_metadata
phos_feature_meta <- motrpac_qc$adipose$`prot-ph`$feature_metadata

#Prcoess prot df
phos_qc_long <- phos_qc %>%
  rownames_to_column("uniprot") %>%
  pivot_longer(-uniprot, names_to = "vialLabel", values_to = "expression")
phos_meta <- phos_meta %>%
  mutate(vialLabel = as.character(vialLabel))
phos_qc_long <- phos_qc_long %>%
  left_join(phos_meta %>% select(vialLabel, pid, timepointDescription), by = "vialLabel")
unique(phos_qc_long$timepointDescription) #"4 hr post"   "Rest 1"      "Pre ex"      "4 hr Rest 3" ""

phos_pre <- phos_qc_long %>%
  filter(timepointDescription %in% c("Rest 1", "Pre ex", "")) %>%

  # Summarize duplicates by averaging expression values
  group_by(uniprot, pid) %>%
  summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%

  # Pivot to wide format with uniphos as rows and pid as columns
  pivot_wider(names_from = pid, values_from = expression) %>%

  # Set uniphos as row names
  column_to_rownames("uniprot")

#### WGCNA pipeline for phoseomics
phos_pre$Zero_count <- rowSums(is.na(phos_pre))
hist(phos_pre$Zero_count)
phos_pre1 = phos_pre[phos_pre$Zero_count<0.25*22,] # filtered out low count genes. From 21022 to 8998 phoseins
wgcna_input_phos <- phos_pre1 %>%
  dplyr::select(-Zero_count)

datExpr2 = data.frame(t(wgcna_input_phos), check.names = FALSE)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# check if all genes/samples are good
gsg = goodSamplesGenes(datExpr2, verbose = 3)
gsg$allOK #[1] TRUE

#Soft thresholding and network topology analysis ####

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr2, powerVector = powers, verbose = 5)

# Module construction
pwr = 12
ngenes = 50

# takes several hours
net2 = blockwiseModules(datExpr2,
                        power = pwr, #determined above
                        maxBlockSize = 20000, # we can set this high because we have a lot of ram. We have ~17000 genes so this will do them all at onec
                        TOMType = "signed", # unsigned was default
                        minModuleSize = ngenes,
                        reassignThreshold = 0.3,
                        mergeCutHeight = 0.3,
                        numericLabels = TRUE, # names modules as numbers isntead of colors
                        pamRespectsDendro = FALSE,
                        saveTOMs = F,
                        nThreads = 48,
                        verbose = 3)

MEs2 = net2$MEs

module_membership2 = as.data.frame(net2$colors)
colnames(module_membership2) = 'module'
module_membership2$ID = row.names(module_membership2)
module_membership2$module <- paste0("Ph", module_membership2$module)
module_membership2$module <- factor(module_membership2$module, levels = paste0("Ph", 0:18))
module_membership2 <- module_membership2[order(module_membership2$module), ]

## Metabolomics
qc_norm_list <- list()

# Loop through all objects starting with "metab"
for (metab_name in names(motrpac_qc$adipose)) {
  if (startsWith(metab_name, "metab")) {
    # Extract qc_norm and sample_metadata for the current metab object
    qc_norm <- motrpac_qc$adipose[[metab_name]]$qc_norm
    sample_metadata <- motrpac_qc$adipose[[metab_name]]$sample_metadata

    # Filter sample_metadata to keep only pre-exercise samples
    pre_exercise_samples <- sample_metadata$Timepoint == "pre_exercise" & sample_metadata$visitcode == "ADU_BAS"
    filtered_metadata <- sample_metadata[pre_exercise_samples, ]

    # Subset qc_norm to include only pre-exercise vialLabels
    pre_exercise_qc_norm <- qc_norm[, colnames(qc_norm) %in% filtered_metadata$vialLabel, drop = FALSE]

    # Rename vialLabels (columns) to their corresponding pid from sample_metadata
    colnames(pre_exercise_qc_norm) <- filtered_metadata$pid[match(colnames(pre_exercise_qc_norm), filtered_metadata$vialLabel)]

    # Add the filtered qc_norm to the list
    qc_norm_list[[metab_name]] <- pre_exercise_qc_norm
  }
}
# Combine all dataframes by rows, allowing for differing columns
all_columns <- unique(unlist(lapply(qc_norm_list, colnames)))
qc_norm_list <- lapply(qc_norm_list, function(df) {
  missing_columns <- setdiff(all_columns, colnames(df))
  df[, missing_columns] <- NA  # Add missing columns filled with NA
  df <- df[, all_columns]  # Reorder columns to match the full set
  return(df)
})

# Combine all dataframes by rows
final_qc_norm <- do.call(rbind, qc_norm_list)

# View the final dataframe
head(final_qc_norm)
cleaned_rownames <- sub("^[^.]+\\.", "", rownames(final_qc_norm))
final_qc_norm$feature_id <- cleaned_rownames

# Aggregate by averaging rows with the same feature_id
final_qc_norm_aggregated <- final_qc_norm %>%
  group_by(feature_id) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    across(where(negate(is.numeric)), first)
  ) %>%
  ungroup()

# Convert back to a dataframe with rownames
final_qc_norm_aggregated <- as.data.frame(final_qc_norm_aggregated)
rownames(final_qc_norm_aggregated) <- final_qc_norm_aggregated$feature_id
final_qc_norm_aggregated$feature_id <- NULL


#### WGCNA pipeline for metabeomics
final_qc_norm_aggregated$Zero_count <- rowSums(is.na(final_qc_norm_aggregated))
hist(final_qc_norm_aggregated$Zero_count)
final_qc_norm_aggregated1 = final_qc_norm_aggregated[final_qc_norm_aggregated$Zero_count<0.25*172,] # filtered out low count metabs. All 676 survived.
wgcna_input_metab <- final_qc_norm_aggregated1 %>%
  dplyr::select(-Zero_count)### 2.0 - Check for outlier/bad genes and samples ####

datExpr3 <- data.frame(t(wgcna_input_metab), check.names = FALSE)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# check if all genes/samples are good
gsg = goodSamplesGenes(datExpr3, verbose = 3)
gsg$allOK #[1] TRUE

# Soft thresholding and network topology analysis ####

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr3, powerVector = powers, verbose = 5)

# Module Construction
pwr = 4
ngenes = 20

# takes several hours
net3 = blockwiseModules(datExpr3,
                        power = pwr, #determined above
                        maxBlockSize = 20000, # we can set this high because we have a lot of ram. We have ~17000 genes so this will do them all at onec
                        TOMType = "signed", # unsigned was default
                        minModuleSize = ngenes,
                        reassignThreshold = 0.3,
                        mergeCutHeight = 0.3,
                        numericLabels = TRUE, # names modules as numbers isntead of colors
                        pamRespectsDendro = FALSE,
                        saveTOMs = F,
                        nThreads = 48,
                        verbose = 3)

MEs3 = net3$MEs

module_membership3 = as.data.frame(net3$colors)
colnames(module_membership3) = 'module'
module_membership3$ID = row.names(module_membership3)
module_membership3$module <- paste0("M", module_membership3$module)
module_membership3$module <- factor(module_membership3$module, levels = paste0("M", 0:7))
module_membership3 <- module_membership3[order(module_membership3$module), ]
## Completion of module construction

### WGCNA downstream analyses
# 1. Transcriptomics - ORA
mm = module_membership %>%
  filter(module != "T0" ) %>%
  droplevels()

signatures = split(mm$ID, mm$module)
bckg <- colnames(datExpr0)
ora_wgcna <- lapply(signatures, function(input_i) {
  run_ORA(input = input_i,
          background = bckg)
}) %>%
  bind_rows(.id = "Module") %>%
  mutate(log10p = -log10(p_value),
         Module = factor(Module, levels = unique(Module)))
ora_wgcna_fil <- ora_wgcna %>%
  filter(adj_p_value<0.05)
top_pathways_t <- ora_wgcna %>%
  filter(database == "GOBP") %>%
  group_by(Module) %>%
  slice_min(p_value, n = 1) %>%  # Select most significant pathway per module
  ungroup() %>%
  pull(set_short)  # Extract pathway names

# Filter the original dataset to keep all rows associated with selected pathways
ora_top_per_module_t <- ora_wgcna %>%
  filter(set_short %in% top_pathways_t) %>%
  mutate(set_short = str_remove(set_short, "^GOBP_"))

heatmap_mat_ora_t <- ora_top_per_module_t %>%
  select(Module, set_short, log10p) %>%
  pivot_wider(names_from = set_short, values_from = log10p, values_fill = NA) %>%
  column_to_rownames("Module") %>%
  as.matrix()
col_fun_ora <- colorRamp2(c(0, max(heatmap_mat_ora_t, na.rm = TRUE)), c("white", "#377EB8"))

adjp_threshold <- 0.05

heat_ora <- Heatmap(
  heatmap_mat_ora_t,
  name = "-log10(P-value)",
  col = col_fun_ora,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "-log10(P-value)"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))
    if (heatmap_mat_ora_t[i, j] >= -log10(adjp_threshold)) {
      grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
    }
  }
)
draw(heat_ora, heatmap_legend_side = "left")

# 2. trait-module (T) correlation
colnames(MEs) <- gsub("ME", "T", colnames(MEs))  # Replace "ME" with "T"
MEs <- MEs[, colnames(MEs) != "T0"]             # Remove column T0

# Reorder columns from T1 to T13
desired_order <- paste0("T", 1:12)              # Generate order from T1 to T12
MEs <- MEs[, desired_order]
# make sure clinical and ME rows are aligned
common_samples <- intersect(rownames(cli_full_z2), rownames(MEs))

# Subset both data frames to keep only the common samples and ensure they are in the same order
cli_full_z2_t <- cli_full_z2[common_samples, ]
MEs <- MEs[common_samples, ]
stopifnot(identical(rownames(cli_full_z2_t), rownames(MEs)))

# correlation
cor = bicorAndPvalue(cli_full_z2_t, MEs)
df = melt(cor$bicor) %>% dplyr::rename(bicor = value)
df$pval = melt(cor$p)$value
df$obs = melt(cor$nObs)$value
df <- df %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

traitXmod = dcast(df, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod) = traitXmod$Var1
traitXmod$Var1 = NULL
traitXmod <- traitXmod %>%
  select(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12)
# set up signif overlay
sig_table = dcast(df, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table) = sig_table$Var1
sig_table$Var1 = NULL
sig_table <- sig_table %>%
  select(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12)

traitXmod_mat_t <- t(traitXmod)
sig_table_mat_t <- t(sig_table)

col_fun_t <- colorRamp2(seq(min(traitXmod_mat_t), max(traitXmod_mat_t), length.out = 10),
                        rev(brewer.pal(10, "RdBu")))
col_group_vec <- column_groups_named[colnames(traitXmod_mat_t)]

column_annotation <- HeatmapAnnotation(
  Clinical_traits = factor(col_group_vec, levels = names(group_colors)),
  col = list(Clinical_traits = group_colors),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

T_count <- module_membership %>%
  filter(module %in% rownames(traitXmod_mat_t)) %>%  # Ensure matching modules
  group_by(module) %>%
  summarise(ngenes = n(), .groups = "drop") %>%
  mutate(module = factor(module, levels = rownames(traitXmod_mat_t))) %>%  # Order modules correctly
  arrange(module)

# Step 3: Create Right Annotation (Barplot)
right_annot_t <- rowAnnotation(
  "Feature Count" = anno_barplot(
    T_count$ngenes,  # Gene counts per module
    gp = gpar(fill = "#377EB8", col = "black"),  # Bar color and outline
    width = unit(3, "cm"),
    axis_param = list(
      at = seq(0, 3000, by = 1000),       # Ticks at 0, 500, ..., 3000
      labels = seq(0, 3000, by = 1000),   # Matching labels
      labels_rot = 0                     # Horizontal labels
    )
  ),
  annotation_name_gp = gpar(fontsize = 10)
)
heat_trait_mod_t <- Heatmap(
  traitXmod_mat_t,
  name = "Correlation",
  col = col_fun_t,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Correlation"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  top_annotation = column_annotation,
  right_annotation = right_annot_t,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))

    if (!is.na(sig_table_mat_t[i, j])) {
      # Get the fill color and convert it to RGB
      fill_rgb <- col2rgb(fill)
      luminance <- (0.299 * fill_rgb[1] + 0.587 * fill_rgb[2] + 0.114 * fill_rgb[3]) / 255

      # Choose black or white based on brightness
      text_color <- if (luminance < 0.5) "white" else "black"

      grid.text(sig_table_mat_t[i, j], x, y, gp = gpar(fontsize = 9, col = text_color))
    }
  }
)

# 3. DA transcripts on modules.
mm_t <- data.frame(ID = names(net$colors), module = net$colors) %>%
  mutate(
    module = factor(paste0("T", module), levels = paste0("T", 0:12))
  ) %>%
  arrange(module)

# Call in motrpac da.
motrpac_da <- load_differential_analysis()
precawg_trans_da <- as.data.frame(motrpac_da$adipose$`transcript-rna-seq`)
precawg_trans_da <- precawg_trans_da %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol), by = "feature_id")

precawg_trans_da_sig_ee <- precawg_trans_da %>%
  filter(contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                         "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                         "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise")) %>%
  filter(adj_p_value<0.05)%>%
  pull(gene_symbol)
mm_t_ee <- mm_t %>%
  mutate(da = ifelse(ID %in% precawg_trans_da_sig_ee, "yes", "no"))

mm_t_summary_ee <- mm_t_ee %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Total count of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0))


#Trans-RE
precawg_trans_da_sig_re <- precawg_trans_da %>%
  filter(contrast %in% c("group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                         "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                         "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise")) %>%
  filter(adj_p_value<0.05)%>%
  pull(gene_symbol)
mm_t_re <- mm_t %>%
  mutate(da = ifelse(ID %in% precawg_trans_da_sig_re, "yes", "no"))

mm_t_summary_re <- mm_t_re %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Total count of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0))

# integrated plot
mm_t_combined <- bind_rows(
  mm_t_summary_ee %>% mutate(Modality = "EE"),
  mm_t_summary_re %>% mutate(Modality = "RE")
)
mm_t_matrix <- mm_t_combined %>%
  filter(module != "T0") %>%  # Exclude "T0", keep T1, T2, ...
  group_by(module, Modality) %>%  # Group by module & Modality
  summarise(yes_count = sum(yes_count), .groups = "drop") %>%  # Ensure unique pairs
  pivot_wider(names_from = Modality, values_from = yes_count, values_fill = 0) %>%
  column_to_rownames("module") %>%
  as.matrix()

# Define color scale
col_fun_t_da <- colorRamp2(
  c(0, max(mm_t_matrix, na.rm = TRUE)),
  c("white", "#377EB8")  # White for low, Blue for high
)
colanno_ex <- columnAnnotation(
  Modality = colnames(mm_t_matrix),  # Assign EE and RE labels
  col = list(Modality = c("EE" = "#d95f02", "RE" = "#1b9e77")),  # Correct color mapping
  annotation_legend_param = list(title = "Modality")
)
# Create heatmap
heat_t_da <- Heatmap(
  mm_t_matrix,
  name = "Yes Count",
  col = col_fun_t_da,
  cluster_rows = FALSE,  # Keep module order fixed
  cluster_columns = FALSE,  # Keep EE and RE order
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "DA count"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),  # Grid lines
  top_annotation = colanno_ex,  # Add column annotation
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))  # Subtle gridlines
    grid.text(mm_t_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))  # Display Yes Count
  }
)

# Draw the heatmap
draw(heat_t_da)

# stitching
wgcna_figure_t <- heat_ora + heat_trait_mod_t + heat_t_da
draw(wgcna_figure_t)

## Proteomics
# 1. ORA
datExpr1_long <- datExpr1 %>%
  as.data.frame() %>%
  rownames_to_column("pid") %>%
  pivot_longer(-pid, names_to = "uniprot", values_to = "expression")

# Step 2: Join with gene_mapping to replace `uniprot` with `gene_symbol`
datExpr1_long <- datExpr1_long %>%
  left_join(gene_mapping, by = "uniprot") %>%
  drop_na(gene_symbol)  # Drop rows where `gene_symbol` is missing
datExpr1_long <- datExpr1_long %>%
  filter(!is.na(expression)) %>%
  group_by(pid, gene_symbol) %>%
  summarise(expression = max(expression, na.rm = TRUE), .groups = "drop")
datExpr1_transformed <- datExpr1_long %>%
  pivot_wider(names_from = gene_symbol, values_from = expression) %>%
  column_to_rownames("pid")
#Removing duplicates left 6618 gene-symbols

mm1 = module_membership1 %>%
  filter(module != "Pr0" ) %>%
  droplevels()
mm1 <- mm1 %>% # converting uniprot to gene_symbol
  left_join(gene_mapping, by = c("ID" = "uniprot")) %>%
  # Replace ID with gene_symbol and drop unneeded columns
  mutate(ID = gene_symbol) %>%
  select(-gene_symbol)
signatures1 = split(mm1$ID, mm1$module)
bckg_pr <- as.character(unique(datExpr1_long$gene_symbol))

ora_wgcna_pr <- lapply(signatures1, function(input_i) {
  run_ORA(input = as.character(input_i),
          background = bckg_pr)
}) %>%
  bind_rows(.id = "Module") %>%
  mutate(log10p = -log10(p_value),
         Module = factor(Module, levels = unique(Module)))

top_pathways_pr <- ora_wgcna_pr %>%
  filter(database == "GOBP") %>%
  group_by(Module) %>%
  slice_min(p_value, n = 1) %>%  # Select most significant pathway per module
  ungroup() %>%
  pull(set_short)  # Extract pathway names

# Filter the original dataset to keep all rows associated with selected pathways
ora_top_per_module_pr <- ora_wgcna_pr %>%
  filter(set_short %in% top_pathways_pr) %>%
  mutate(set_short = str_remove(set_short, "^GOBP_"))

heatmap_mat_ora_pr <- ora_top_per_module_pr %>%
  select(Module, set_short, log10p) %>%
  pivot_wider(names_from = set_short, values_from = log10p, values_fill = NA) %>%
  column_to_rownames("Module") %>%
  as.matrix()
col_fun_ora_pr <- colorRamp2(c(0, max(heatmap_mat_ora_pr, na.rm = TRUE)), c("white", "#228833"))

adjp_threshold <- 0.05

heat_ora_pr <- Heatmap(
  heatmap_mat_ora_pr,
  name = "-log10(P-value)",
  col = col_fun_ora_pr,
  cluster_rows = TRUE,  # Cluster modules
  cluster_columns = TRUE,  # Cluster pathways
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,  # Ensure text rotates properly
  column_names_gp = gpar(fontsize = 9),  # Keep font size, but remove rotation from here
  heatmap_legend_param = list(title = "-log10(P-value)"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))
    if (heatmap_mat_ora_pr[i, j] >= -log10(adjp_threshold)) {
      grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
    }
  }
)
draw(heat_ora_pr)

## trait-module (Pr) correlation
colnames(MEs1) <- gsub("ME", "Pr", colnames(MEs1))  # Replace "ME" with "T"
MEs1 <- MEs1[, colnames(MEs1) != "Pr0"]             # Remove column T0

# Reorder columns from T1 to T13
desired_order <- paste0("Pr", 1:14)              # Generate order from T1 to T12
MEs1 <- MEs1[, desired_order]
# make sure clinical and ME rows are aligned
common_samples <- intersect(rownames(cli_full_z2), rownames(MEs1))

# Subset both data frames to keep only the common samples and ensure they are in the same order
cli_full_z2_pr <- cli_full_z2[common_samples, ]
MEs1 <- MEs1[common_samples, ]
stopifnot(identical(rownames(cli_full_z2_pr), rownames(MEs1)))

# correlation
cor1 = bicorAndPvalue(cli_full_z2_pr, MEs1) # alternative biweight midcor
df1 = melt(cor1$bicor) %>% dplyr::rename(bicor = value)
df1$pval = melt(cor1$p)$value
df1$obs = melt(cor1$nObs)$value
df1 <- df1 %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
traitXmod_pr = dcast(df1, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod_pr) = traitXmod_pr$Var1
traitXmod_pr$Var1 = NULL
traitXmod_pr <- traitXmod_pr %>%
  select(Pr1, Pr2, Pr3, Pr4, Pr5, Pr6, Pr7, Pr8, Pr9, Pr10, Pr11, Pr12,  Pr13,  Pr14)

# set up signif overlay
sig_table_pr = dcast(df1, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_pr) = sig_table_pr$Var1
sig_table_pr$Var1 = NULL
sig_table_pr <- sig_table_pr %>%
  select(Pr1, Pr2, Pr3, Pr4, Pr5, Pr6, Pr7, Pr8, Pr9, Pr10, Pr11, Pr12,  Pr13,  Pr14)

traitXmod_mat_pr <- t(traitXmod_pr)
sig_table_mat_pr <- t(sig_table_pr)

col_fun_pr <- colorRamp2(seq(min(traitXmod_mat_pr), max(traitXmod_mat_pr), length.out = 10),
                         rev(brewer.pal(10, "RdBu")))

Pr_count <- module_membership1 %>%
  filter(module %in% rownames(traitXmod_mat_pr)) %>%  # Ensure matching modules
  group_by(module) %>%
  summarise(ngenes = n(), .groups = "drop") %>%
  mutate(module = factor(module, levels = rownames(traitXmod_mat_pr))) %>%  # Order modules correctly
  arrange(module)

# Step 3: Create Right Annotation (Barplot)
right_annot_pr <- rowAnnotation(
  "Feature Count" = anno_barplot(
    Pr_count$ngenes,  # Gene counts per module
    gp = gpar(fill = "#228833", col = "black"),  # Bar color and outline
    width = unit(3, "cm"),
    axis_param = list(
      at = seq(0, 1200, by = 400),     # Ticks at 0, 400, 800, 1200
      labels = seq(0, 1200, by = 400),
      labels_rot = 0                   # Horizontal labels
    )
  ),
  annotation_name_gp = gpar(fontsize = 10)
)
heat_trait_mod_pr <- Heatmap(
  traitXmod_mat_pr,
  name = "Correlation",
  col = col_fun_pr,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Correlation"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  top_annotation = column_annotation,
  right_annotation = right_annot_pr,  # Attach gene count barplot to the right
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))

    if (!is.na(sig_table_mat_pr[i, j])) {
      fill_rgb <- col2rgb(fill)
      luminance <- (0.299 * fill_rgb[1] + 0.587 * fill_rgb[2] + 0.114 * fill_rgb[3]) / 255
      text_color <- if (luminance < 0.5) "white" else "black"

      grid.text(sig_table_mat_pr[i, j], x, y, gp = gpar(fontsize = 9, col = text_color))
    }
  }
)

## Protein DA
mm_pr <- data.frame(feature_id = names(net1$colors), module = net1$colors) %>%
  mutate(
    module = factor(paste0("Pr", module), levels = paste0("Pr", 0:14))
  ) %>%
  arrange(module)
mm_pr$feature_id <- gsub("\\.", "-", mm_pr$feature_id)
mm_pr <- mm_pr %>%
  left_join(
    HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol),
    by = "feature_id"
  )
#EE
precawg_prot_da <- as.data.frame(motrpac_da$adipose$`prot-pr`)
precawg_prot_da <- precawg_prot_da %>%
  left_join(
    HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, uniprot, gene_symbol),
    by = "feature_id"
  )
precawg_prot_da_sig_ee <- precawg_prot_da %>%
  filter(contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise") %>%
  filter(adj_p_value<0.05)%>%
  pull(gene_symbol)
mm_pr_ee <- mm_pr %>%
  mutate(da = ifelse(gene_symbol %in% precawg_prot_da_sig_ee, "yes", "no"))

mm_pr_summary_ee <- mm_pr_ee %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Calculate total number of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0)) # Ensure all modules are included

non_yes_count <- mm_pr_ee %>%
  filter(module == "Pr5") %>%  # Filter for module Pr5
  summarize(non_yes_count = sum(da == "yes", na.rm = TRUE))  # Count "yes" values
mm_pr_ee_filtered <- mm_pr_ee %>%
  filter(module == "Pr5", da == "yes")


#RE
precawg_prot_da_sig_re <- precawg_prot_da %>%
  filter(contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise") %>%
  filter(adj_p_value<0.05)%>%
  pull(gene_symbol)
mm_pr_re <- mm_pr %>%
  mutate(da = ifelse(gene_symbol %in% precawg_prot_da_sig_re, "yes", "no"))

mm_pr_summary_re <- mm_pr_re %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Calculate total number of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0)) # Ensure all modules are included

# integrated plot
mm_pr_combined <- bind_rows(
  mm_pr_summary_ee %>% mutate(Modality = "EE"),
  mm_pr_summary_re %>% mutate(Modality = "RE")
)
mm_pr_matrix <- mm_pr_combined %>% #This is from wgcna_downstream.R
  filter(module != "Pr0") %>%  # Exclude "T0", keep T1, T2, ...
  group_by(module, Modality) %>%  # Group by module & Modality
  summarise(yes_count = sum(yes_count), .groups = "drop") %>%  # Ensure unique pairs
  pivot_wider(names_from = Modality, values_from = yes_count, values_fill = 0) %>%
  column_to_rownames("module") %>%
  as.matrix()

# Define color scale
col_fun_pr_da <- colorRamp2(
  c(0, max(mm_pr_matrix, na.rm = TRUE)),
  c("white", "#228833")  # White for low, Blue for high
)
colanno_ex <- columnAnnotation(
  Modality = colnames(mm_pr_matrix),  # Assign EE and RE labels
  col = list(Modality = c("EE" = "#d95f02", "RE" = "#1b9e77")),  # Correct color mapping
  annotation_legend_param = list(title = "Modality")
)
# Create heatmap
heat_pr_da <- Heatmap(
  mm_pr_matrix,
  name = "Yes Count",
  col = col_fun_pr_da,
  cluster_rows = FALSE,  # Keep module order fixed
  cluster_columns = FALSE,  # Keep EE and RE order
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "DA count"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),  # Grid lines
  top_annotation = colanno_ex,  # Add column annotation
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))  # Subtle gridlines
    grid.text(mm_pr_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))  # Display Yes Count
  }
)

# Draw the heatmap
draw(heat_pr_da)

# stitching
wgcna_figure_pr <- heat_ora_pr + heat_trait_mod_pr + heat_pr_da
draw(wgcna_figure_pr)

###Phospho
# ORA
module_membership2 <- module_membership2 %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% select(feature_id, gene_symbol),
            by = c("ID" = "feature_id"))
mm2 = module_membership2 %>%
  filter(module != "Ph0" ) %>%
  droplevels()
signatures2 = split(mm2$gene_symbol, mm2$module)


bckg_ph <- as.character(unique(module_membership2$gene_symbol))
ora_wgcna_ph <- lapply(signatures2, function(input_i) {
  run_ORA(input = as.character(input_i),
          background = bckg_ph,
          overlap_cutoff = 0.3)
}) %>%
  bind_rows(.id = "Module") %>%
  mutate(log10p = -log10(p_value),
         Module = factor(Module, levels = unique(Module)))
top_pathways_ph <- ora_wgcna_ph %>%
  filter(database == "GOBP") %>%
  group_by(Module) %>%
  slice_min(p_value, n = 1) %>%  # Select most significant pathway per module
  ungroup() %>%
  pull(set_short)  # Extract pathway names

# Filter the original dataset to keep all rows associated with selected pathways
ora_top_per_module_ph <- ora_wgcna_ph %>%
  filter(set_short %in% top_pathways_ph) %>%
  mutate(set_short = str_remove(set_short, "^GOBP_"))

heatmap_mat_ora_ph <- ora_top_per_module_ph %>%
  select(Module, set_short, log10p) %>%
  pivot_wider(names_from = set_short, values_from = log10p, values_fill = NA) %>%
  column_to_rownames("Module") %>%
  as.matrix()
col_fun_ora_ph <- colorRamp2(c(0, max(heatmap_mat_ora_ph, na.rm = TRUE)), c("white", "#F3A02B"))

adjp_threshold <- 0.05

heat_ora_ph <- Heatmap(
  heatmap_mat_ora_ph,
  name = "-log10(P-value)",
  col = col_fun_ora_ph,
  cluster_rows = TRUE,  # Cluster modules
  cluster_columns = TRUE,  # Cluster pathways
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,  # Ensure text rotates properly
  column_names_gp = gpar(fontsize = 9),  # Keep font size, but remove rotation from here
  heatmap_legend_param = list(title = "-log10(P-value)"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))
    if (heatmap_mat_ora_ph[i, j] >= -log10(adjp_threshold)) {
      grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
    }
  }
)

# module (ph)-trait correlation
#make sure clinical and ME rows are aligned
colnames(MEs2) <- gsub("ME", "Ph", colnames(MEs2))  # Replace "ME" with "T"
MEs2 <- MEs2[, colnames(MEs2) != "Ph0"]             # Remove column T0

# Reorder columns from T1 to T13
desired_order <- paste0("Ph", 1:18)              # Generate order from T1 to T13
MEs2 <- MEs2[, desired_order] #proteomics ME ready

common_samples <- intersect(rownames(cli_full_z2), rownames(MEs2))

# Subset both data frames to keep only the common samples and ensure they are in the same order
cli_full_z2_ph <- (cli_full_z2)[common_samples, ]
MEs2 <- MEs2[common_samples, ]
stopifnot(identical(rownames(cli_full_z2_ph), rownames(MEs2)))

# correlation
cor2 = bicorAndPvalue(cli_full_z2_ph, MEs2) # alternative biweight midcor
df2 = melt(cor2$bicor) %>% dplyr::rename(bicor = value)
df2$pval = melt(cor2$p)$value
df2$obs = melt(cor2$nObs)$value
df2 <- df2 %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

traitXmod_ph = dcast(df2, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod_ph) = traitXmod_ph$Var1
traitXmod_ph$Var1 = NULL
traitXmod_ph <- traitXmod_ph %>%
  select(Ph1, Ph2, Ph3, Ph4, Ph5, Ph6, Ph7, Ph8, Ph9, Ph10, Ph11, Ph12,  Ph13,  Ph14, Ph15,Ph16,Ph17,Ph18)
# set up signif overlay
sig_table_ph = dcast(df2, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_ph) = sig_table_ph$Var1
sig_table_ph$Var1 = NULL
sig_table_ph <- sig_table_ph %>%
  select(Ph1, Ph2, Ph3, Ph4, Ph5, Ph6, Ph7, Ph8, Ph9, Ph10, Ph11, Ph12,  Ph13,  Ph14, Ph15,Ph16,Ph17,Ph18)

traitXmod_mat_ph <- t(traitXmod_ph)
sig_table_mat_ph <- t(sig_table_ph)

col_fun_ph <- colorRamp2(seq(min(traitXmod_mat_ph), max(traitXmod_mat_ph), length.out = 10),
                         rev(brewer.pal(10, "RdBu")))

Ph_count <- module_membership2%>%
  filter(module %in% rownames(traitXmod_mat_ph)) %>%  # Ensure matching modules
  group_by(module) %>%
  summarise(ngenes = n(), .groups = "drop") %>%
  mutate(module = factor(module, levels = rownames(traitXmod_mat_ph))) %>%  # Order modules correctly
  arrange(module)

# Step 3: Create Right Annotation (Barplot)
right_annot_ph <- rowAnnotation(
  "Feature Count" = anno_barplot(
    Ph_count$ngenes,  # Gene counts per module
    gp = gpar(fill = "#F3A02B", col = "black"),  # Bar color and outline
    width = unit(3, "cm"),
    axis_param = list(
      at = seq(0, 1200, by = 400),  # Tick marks at 0, 400, 800, 1200
      labels = seq(0, 1200, by = 400),
      labels_rot = 0  # Horizontal labels
    )
  ),
  annotation_name_gp = gpar(fontsize = 10)
)
heat_trait_mod_ph <- Heatmap(
  traitXmod_mat_ph,
  name = "Correlation",
  col = col_fun_ph,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Correlation"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  top_annotation = column_annotation,
  right_annotation = right_annot_ph,  # Attach gene count barplot to the right
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))

    if (!is.na(sig_table_mat_ph[i, j])) {
      fill_rgb <- col2rgb(fill)
      luminance <- (0.299 * fill_rgb[1] + 0.587 * fill_rgb[2] + 0.114 * fill_rgb[3]) / 255
      text_color <- if (luminance < 0.5) "white" else "black"

      grid.text(sig_table_mat_ph[i, j], x, y, gp = gpar(fontsize = 9, col = text_color))
    }
  }
)
# phospho DA
precawg_phos_da <- as.data.frame(motrpac_da$adipose$`prot-ph`)
precawg_phos_da <- precawg_phos_da %>%
  # Join with HUMAN_FEATURE_TO_GENE to map feature_id to gene_symbol
  left_join(
    HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol),
    by = "feature_id"
  ) %>%
  # Process the feature_id after mapping
  mutate(
    # Extract base_feature_id (protein ID) before the first underscore
    base_feature_id = sub("_.*$", "", feature_id),

    # Extract phosphosite(s) from everything after the first underscore
    phosphosites = sub("^.*_", "", feature_id) %>%
      gsub("s$", "", .) %>%  # Remove trailing 's'
      gsub("([STY]\\d+)[a-zA-Z]", "\\1;", .) %>%  # Replace any extra letters after phosphosites with ';'
      gsub(";$", "", .),  # Remove trailing semicolon

    # Create gene_symbol_with_phosphosite
    gene_symbol_with_phosphosite = ifelse(
      !is.na(gene_symbol),
      paste0(gene_symbol, "-", phosphosites),  # Combine gene symbol and cleaned phosphosites
      feature_id  # Fallback to feature_id if gene_symbol is NA
    )
  )

mm_ph <- data.frame(feature_id = names(net2$colors), module = net2$colors) %>%
  mutate(
    module = factor(paste0("Ph", module), levels = paste0("Ph", 0:18))
  ) %>%
  arrange(module)
mm_ph <- mm_ph %>%
  # Join with HUMAN_FEATURE_TO_GENE to map feature_id to gene_symbol
  left_join(
    HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol),
    by = "feature_id"
  ) %>%
  # Process the feature_id after mapping
  mutate(
    # Split feature_id into base_feature_id and phosphosites based on the last `_`
    base_feature_id = sub("_.*$", "", feature_id),  # Extract part before the last underscore
    phosphosites = sub("^.*_", "", feature_id) %>%  # Extract part after the last underscore
      gsub("s$", "", .) %>%                         # Remove trailing 's'
      gsub("([STY]\\d+)[a-zA-Z]", "\\1;", .) %>%   # Replace any lowercase letters after phosphosites with ';'
      gsub(";$", "", .)                            # Remove trailing semicolon
  ) %>%
  # Create gene_symbol_with_phosphosite
  mutate(
    gene_symbol_with_phosphosite = ifelse(
      !is.na(gene_symbol),
      paste0(gene_symbol, "-", phosphosites),  # Combine gene symbol and cleaned phosphosites
      feature_id  # Fallback to feature_id if gene_symbol is NA
    )
  )
precawg_phos_da_sig_ee <- precawg_phos_da %>%
  filter(contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise") %>%
  filter(adj_p_value<0.05)%>%
  pull(feature_id)

mm_ph_ee <- mm_ph %>%
  mutate(da = ifelse(feature_id %in% precawg_phos_da_sig_ee, "yes", "no"))

mm_ph_summary_ee <- mm_ph_ee %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Calculate total number of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0)) # Ensure all modules are included

#RE
precawg_phos_da_sig_re <- precawg_phos_da %>%
  filter(contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise") %>%
  filter(adj_p_value<0.05)%>%
  pull(feature_id)
mm_ph_re <- mm_ph %>%
  mutate(da = ifelse(feature_id %in% precawg_phos_da_sig_re, "yes", "no"))

mm_ph_summary_re <- mm_ph_re %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Calculate total number of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0)) # Ensure all modules are included

# integrated plot
mm_ph_combined <- bind_rows(
  mm_ph_summary_ee %>% mutate(Modality = "EE"),
  mm_ph_summary_re %>% mutate(Modality = "RE")
)
mm_ph_matrix <- mm_ph_combined %>% #This is from wgcna_downstream.R
  filter(module != "Ph0") %>%  # Exclude "T0", keep T1, T2, ...
  group_by(module, Modality) %>%  # Group by module & Modality
  summarise(yes_count = sum(yes_count), .groups = "drop") %>%  # Ensure unique pairs
  pivot_wider(names_from = Modality, values_from = yes_count, values_fill = 0) %>%
  column_to_rownames("module") %>%
  as.matrix()

# Define color scale
col_fun_ph_da <- colorRamp2(
  c(0, max(mm_ph_matrix, na.rm = TRUE)),
  c("white", "#F3A02B")  # White for low, Blue for high
)
colanno_ex <- columnAnnotation(
  Modality = colnames(mm_ph_matrix),  # Assign EE and RE labels
  col = list(Modality = c("EE" = "#d95f02", "RE" = "#1b9e77")),  # Correct color mapping
  annotation_legend_param = list(title = "Modality")
)
# Create heatmap
heat_ph_da <- Heatmap(
  mm_ph_matrix,
  name = "Yes Count",
  col = col_fun_ph_da,
  cluster_rows = FALSE,  # Keep module order fixed
  cluster_columns = FALSE,  # Keep EE and RE order
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "DA count"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),  # Grid lines
  top_annotation = colanno_ex,  # Add column annotation
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))  # Subtle gridlines
    grid.text(mm_ph_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))  # Display Yes Count
  }
)

# stitching
wgcna_figure_ph <- heat_ora_ph + heat_trait_mod_ph + heat_ph_da
draw(wgcna_figure_ph)


##Metabolomics
# 1. ORA
mm3 = module_membership3 %>%
  filter(module != "M0" ) %>%
  droplevels()
signatures3 = split(mm3$ID, mm3$module)
bckg_m <- as.character(colnames(datExpr3))

ora_wgcna_m <- lapply(signatures3, function(input_i) {
  run_ORA(input = as.character(input_i),
          background = bckg_m)
}) %>%
  bind_rows(.id = "Module") %>%
  mutate(log10p = -log10(p_value),
         Module = factor(Module, levels = unique(Module)))

top_pathways_m <- ora_wgcna_m %>%
  group_by(Module) %>%
  slice_min(p_value, n = 2) %>%  # Select most significant pathway per module
  ungroup() %>%
  pull(set_short)  # Extract pathway names

# Filter the original dataset to keep all rows associated with selected pathways
ora_top_per_module_m <- ora_wgcna_m %>%
  filter(set_short %in% top_pathways_m)

heatmap_mat_ora_m <- ora_top_per_module_m %>%
  select(Module, set_short, log10p) %>%
  pivot_wider(names_from = set_short, values_from = log10p, values_fill = NA) %>%
  column_to_rownames("Module") %>%
  as.matrix()
col_fun_ora_m <- colorRamp2(c(0, max(heatmap_mat_ora_m, na.rm = TRUE)), c("white", "#6D4B08"))

adjp_threshold <- 0.05

heat_ora_m <- Heatmap(
  heatmap_mat_ora_m,
  name = "-log10(P-value)",
  col = col_fun_ora_m,
  cluster_rows = TRUE,  # Cluster modules
  cluster_columns = TRUE,  # Cluster pathways
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,  # Ensure text rotates properly
  column_names_gp = gpar(fontsize = 9),  # Keep font size, but remove rotation from here
  heatmap_legend_param = list(title = "-log10(P-value)"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))
    if (heatmap_mat_ora_m[i, j] >= -log10(adjp_threshold)) {
      grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
    }
  }
)
draw(heat_ora_m)

## module (M)-trait correlation
colnames(MEs3) <- gsub("ME", "M", colnames(MEs3))  # Replace "ME" with "M"
MEs3 <- MEs3[, colnames(MEs3) != "M0"]             # Remove column T0

# Reorder columns from T1 to T13
desired_order <- paste0("M", 1:7)              # Generate order from T1 to T13
MEs3 <- MEs3[, desired_order] #proteomics ME ready
common_samples <- intersect(rownames(cli_full_z2), rownames(MEs3))

# Subset both data frames to keep only the common samples and ensure they are in the same order
cli_full_z2_m <- cli_full_z2[common_samples, ]
MEs3 <- MEs3[common_samples, ]
stopifnot(identical(rownames(cli_full_z2_m), rownames(MEs3)))

# Correlation
cor3 = bicorAndPvalue(cli_full_z2_m, MEs3) # alternative biweight midcor
df3 = melt(cor3$bicor) %>% dplyr::rename(bicor = value)
df3$pval = melt(cor3$p)$value
df3$obs = melt(cor3$nObs)$value
df3 <- df3 %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

traitXmod_m = dcast(df3, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod_m) = traitXmod_m$Var1
traitXmod_m$Var1 = NULL
traitXmod_m <- traitXmod_m %>%
  select(M1, M2, M3, M4, M5, M6, M7)
# set up signif overlay
sig_table_m = dcast(df3, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_m) = sig_table_m$Var1
sig_table_m$Var1 = NULL
sig_table_m <- sig_table_m %>%
  select(M1, M2, M3, M4, M5, M6, M7)

traitXmod_mat_m <- t(traitXmod_m)
sig_table_mat_m <- t(sig_table_m)

col_fun_m <- colorRamp2(seq(min(traitXmod_mat_m), max(traitXmod_mat_m), length.out = 10),
                        rev(brewer.pal(10, "RdBu")))

M_count <- module_membership3%>%
  filter(module %in% rownames(traitXmod_mat_m)) %>%  # Ensure matching modules
  group_by(module) %>%
  summarise(ngenes = n(), .groups = "drop") %>%
  mutate(module = factor(module, levels = rownames(traitXmod_mat_m))) %>%  # Order modules correctly
  arrange(module)

# Step 3: Create Right Annotation (Barplot)
right_annot_m <- rowAnnotation(
  "Feature Count" = anno_barplot(
    M_count$ngenes,  # Gene counts per module
    gp = gpar(fill = "#6D4B08", col = "black"),  # Bar color and outline
    width = unit(3, "cm"),
    axis_param = list(
      at = seq(0, 120, by = 40),      # Ticks at 0, 40, 80, 120
      labels = seq(0, 120, by = 40),
      labels_rot = 0                  # Horizontal labels
    )
  ),
  annotation_name_gp = gpar(fontsize = 10)
)
heat_trait_mod_m <- Heatmap(
  traitXmod_mat_m,
  name = "Correlation",
  col = col_fun_m,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Correlation"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),
  top_annotation = column_annotation,
  right_annotation = right_annot_m,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))

    if (!is.na(sig_table_mat_m[i, j])) {
      fill_rgb <- col2rgb(fill)
      luminance <- (0.299 * fill_rgb[1] + 0.587 * fill_rgb[2] + 0.114 * fill_rgb[3]) / 255
      text_color <- if (luminance < 0.5) "white" else "black"

      grid.text(sig_table_mat_m[i, j], x, y, gp = gpar(fontsize = 9, col = text_color))
    }
  }
)

## DA metabolomics
precawg_metab_da <- motrpac_da$adipose %>%
  keep(~ grepl("metab", .x$assay[1])) %>% # Select elements containing "metab" in assay column
  bind_rows()
mm_m <- tibble(
  feature_id = names(net3$colors),  # Extract feature names (metabolite IDs)
  module = as.character(net3$colors)  # Extract module assignments as characters
) %>%
  mutate(
    module = factor(paste0("M", module), levels = paste0("M", 0:7))  # Format module labels
  ) %>%
  arrange(module)

#EE
precawg_met_da_sig_ee <- precawg_metab_da %>%
  filter(contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                         "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                         "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise")) %>%
  filter(adj_p_value<0.05)%>%
  pull(feature_id)
mm_m_ee <- mm_m %>%
  mutate(da = ifelse(feature_id %in% precawg_met_da_sig_ee, "yes", "no"))

mm_m_summary_ee <- mm_m_ee %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Calculate total number of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0)) # Ensure all modules are included

#RE
precawg_met_da_sig_re <- precawg_metab_da %>%
  filter(contrast %in% c("group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                         "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                         "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise")) %>%
  filter(adj_p_value<0.05)%>%
  pull(feature_id)
mm_m_re <- mm_m %>%
  mutate(da = ifelse(feature_id %in% precawg_met_da_sig_re, "yes", "no"))

mm_m_summary_re <- mm_m_re %>%
  group_by(module) %>%
  summarise(
    yes_count = sum(da == "yes"),
    total_count = n() # Calculate total number of IDs in each module
  ) %>%
  mutate(percentage = (yes_count / total_count) * 100) %>% # Calculate percentage
  complete(module, fill = list(yes_count = 0, total_count = 0, percentage = 0)) # Ensure all modules are included

# integrated plot
mm_m_combined <- bind_rows(
  mm_m_summary_ee %>% mutate(Modality = "EE"),
  mm_m_summary_re %>% mutate(Modality = "RE")
)

mm_m_matrix <- mm_m_combined %>% #This is from wgcna_downstream.R
  filter(module != "M0") %>%  # Exclude "T0", keep T1, T2, ...
  group_by(module, Modality) %>%  # Group by module & Modality
  summarise(yes_count = sum(yes_count), .groups = "drop") %>%  # Ensure unique pairs
  pivot_wider(names_from = Modality, values_from = yes_count, values_fill = 0) %>%
  column_to_rownames("module") %>%
  as.matrix()

# Define color scale
col_fun_m_da <- colorRamp2(
  c(0, max(mm_m_matrix, na.rm = TRUE)),
  c("white", "#6D4B08")  # White for low, Blue for high
)
colanno_ex <- columnAnnotation(
  Modality = colnames(mm_m_matrix),  # Assign EE and RE labels
  col = list(Modality = c("EE" = "#d95f02", "RE" = "#1b9e77")),  # Correct color mapping
  annotation_legend_param = list(title = "Modality")
)
# Create heatmap
heat_m_da <- Heatmap(
  mm_m_matrix,
  name = "Yes Count",
  col = col_fun_m_da,
  cluster_rows = FALSE,  # Keep module order fixed
  cluster_columns = FALSE,  # Keep EE and RE order
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "DA count"),
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 0.5),  # Grid lines
  top_annotation = colanno_ex,  # Add column annotation
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.3))  # Subtle gridlines
    grid.text(mm_m_matrix[i, j], x, y, gp = gpar(fontsize = 10, col = "black"))  # Display Yes Count
  }
)

# Draw the heatmap
draw(heat_m_da)

# stitching
wgcna_figure_m <- heat_ora_m + heat_trait_mod_m + heat_m_da
draw(wgcna_figure_m)

# Figure 5B
draw(wgcna_figure_t)
draw(wgcna_figure_pr)
draw(wgcna_figure_ph)
draw(wgcna_figure_m)



### Figures 5C-D
kME_ph = signedKME(datExpr2, MEs2) #Extracting module connectivity (kME) for each feature in phospho module.
kME_ph <- rownames_to_column(kME_ph, var = "feature_id")
kME_ph <- kME_ph %>%
  # Join with HUMAN_FEATURE_TO_GENE to map feature_id to gene_symbol
  left_join(
    HUMAN_FEATURE_TO_GENE %>% dplyr::select(feature_id, gene_symbol),
    by = "feature_id"
  ) %>%
  # Create gene_symbol_with_phosphosite
  mutate(
    gene_symbol_with_phosphosite = ifelse(
      !is.na(gene_symbol),
      paste0(
        gene_symbol,
        "-",
        gsub("^[^_]*_", "", feature_id) %>%          # Extract portion after the last underscore
          gsub("([0-9]+)[a-zA-Z]", "\\1;", .) %>%   # Insert semicolon after each number-letter pair
          gsub("[a-zA-Z]+$", "", .) %>%            # Remove trailing letters
          gsub(";$", "", .)                        # Remove trailing semicolon
      ),
      feature_id  # Fallback to feature_id if gene_symbol is NA
    )
  )

# Label phosphosite whether it is a DA or not.
mm_ph <- mm_ph %>%
  mutate(
    da = case_when(
      feature_id %in% precawg_phos_da_sig_ee &
        feature_id %in% precawg_phos_da_sig_re ~ "EE & RE",
      feature_id %in% precawg_phos_da_sig_ee ~ "EE",
      feature_id %in% precawg_phos_da_sig_re ~ "RE",
      TRUE ~ NA_character_  # Default value for unmatched cases
    )
  )
# extract ph10 information
ph10_da <- mm_ph %>%
  # Filter for module "Ph10" and non-NA `da`
  filter(module == "Ph10", !is.na(da)) %>%
  # Join with kME_ph to include kME10 values
  left_join(kME_ph %>% select(feature_id, kME10), by = "feature_id") %>%
  # Separate into positive and negative lists
  mutate(kME_category = ifelse(kME10 > 0, "positive", "negative"))
ph10_features <- mm_ph %>%
  filter(module == "Ph10") %>%
  pull(gene_symbol_with_phosphosite)
ph10_genes <- kME_ph %>%
  filter(gene_symbol_with_phosphosite %in% ph10_features)
ph10_genes <- ph10_genes %>%
  mutate(DA_status = ifelse(gene_symbol_with_phosphosite %in% ph10_da$gene_symbol_with_phosphosite, "DA", "non-DA"))


# Step 2: Safely create unique levels
gene_levels <- ph10_genes %>%
  group_by(gene_symbol_with_phosphosite) %>%
  slice_min(order_by = kME10, n = 1, with_ties = FALSE) %>%  # keep 1 row per gene
  ungroup() %>%
  arrange(kME10) %>%
  pull(gene_symbol_with_phosphosite)

# Step 3: Apply factor levels
ph10_genes$gene_symbol_with_phosphosite <- factor(ph10_genes$gene_symbol_with_phosphosite, levels = gene_levels)
ph10_genes <- ph10_genes %>%
  mutate(DA_status = ifelse(gene_symbol_with_phosphosite %in% ph10_da$gene_symbol_with_phosphosite, "DA", "non-DA")) %>%
  arrange(DA_status == "DA")

# Figure 5C (left)
#pdf("Ph10_scatter.pdf", width = 5, height = 4)
ggplot(ph10_genes, aes(x = gene_symbol_with_phosphosite, y = kME10)) +
  geom_point(aes(color = DA_status), size = 3, alpha = 0.8) +
  geom_label_repel(
    data = ph10_genes %>% filter(DA_status == "DA"),
    aes(label = gene_symbol_with_phosphosite),
    size = 4,
    color = "black",
    fill = "white",
    point.padding = 0.6,        # Pull label away from point
    box.padding = 0.6,          # Padding around label box
    label.padding = unit(0.25, "lines"),  # Text-to-box padding
    nudge_y = 0.05,             # Slight vertical offset
    segment.curvature = -0.2,
    segment.ncp = 5,
    segment.angle = 25,
    force = 2,
    max.overlaps = Inf
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_color_manual(values = c("DA" = "red", "non-DA" = "grey"), name = "Feature Status") +
  labs(
    title = "kME for Ph10 phosphosites",
    x = NULL,
    y = "kME"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#dev.off()

# Convert your palette into a color function
color_palette <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50)
col_fun <- circlize::colorRamp2(seq(-1, 1, length.out = 50), color_palette)

### multiple correlation : Ph10 DA with HOMA IR, ADIPO IR, and HbA1c
common_ids <- intersect(colnames(phos_pre), rownames(cli_full))

# Subset data for common IDs
phos_pre_feat <- phos_pre[ph10_da$feature_id, common_ids, drop = FALSE] #P15104_S343 is GLUL-S343. Q9UHB6_S343s is LIMA1-S343
clinic_feat <- cli_full[common_ids, c("HOMA_IR", "Adipo_IR", "HbA1c")]
clinic_feat <- scale(clinic_feat)

phos_pre_feat <- as.matrix(phos_pre_feat)  # Convert to matrix
clinic_feat <- as.matrix(clinic_feat)  # Convert to matrix

# Compute spearman correlation and p-values
cor_results <- corAndPvalue(t(phos_pre_feat), clinic_feat, use = "pairwise.complete.obs", method = "spearman")

cor_matrix <- cor_results$cor
pval_matrix <- cor_results$p

ordered_indices <- order(cor_matrix[, "HOMA_IR"], decreasing = TRUE)
cor_matrix <- cor_matrix[ordered_indices, ]

# Apply the same order to pval_matrix
pval_matrix <- pval_matrix[ordered_indices, ]

# Convert p-values to significance markers
sig_markers <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
sig_markers[pval_matrix < 0.05] <- "*"
sig_markers[pval_matrix < 0.01] <- "**"
sig_markers[pval_matrix < 0.001] <- "***"

rownames(cor_matrix) <- ph10_da$gene_symbol_with_phosphosite[match(rownames(cor_matrix), ph10_da$feature_id)]

# Correlational heatmap of ph10 DA genes with IR.
ph10_ir <- Heatmap(
  cor_matrix,
  col = col_fun,
  name = "Spearman",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_names_side = "left",
  cluster_columns = FALSE,
  show_heatmap_legend = TRUE,
  row_dend_side = "left",
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(sig_markers[i, j]) && sig_markers[i, j] != "") {
      # Convert fill to RGB
      rgb_vals <- col2rgb(fill) / 255
      luminance <- sum(rgb_vals * c(0.299, 0.587, 0.114))  # perceptual brightness

      text_color <- ifelse(luminance < 0.5, "white", "black")

      grid.text(
        sig_markers[i, j],
        x = x, y = y,
        gp = gpar(fontsize = 11, col = text_color)
      )
    }
  },
  column_names_gp = gpar(fontsize = 10)
)

# Custom function for generating 4hr only feature logFC heatmap.
create_4hr_heatmap <- function(data, contrasts, features = NULL, gene_col, ordered_rows = NULL) {
  # Step 1: Filter and annotate
  heatmap_data <- data %>%
    filter(contrast %in% contrasts) %>%
    mutate(
      group = ifelse(grepl("^group_timepointADUEndur", contrast), "EE", "RE"),
      timeline = "4hrPost",
      contrast = factor(contrast, levels = contrasts)
    ) %>%
    arrange(contrast)

  # Optional feature filtering
  if (!is.null(features)) {
    heatmap_data <- heatmap_data %>%
      filter(!!sym(gene_col) %in% features)
  }

  # Deduplicate based on *lowest* adjusted p-value
  heatmap_data <- heatmap_data %>%
    group_by(!!sym(gene_col), contrast) %>%
    slice_min(order_by = adj_p_value, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Step 2: Pivot to matrix
  logFC_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, logFC) %>%
    pivot_wider(
      names_from = contrast,
      values_from = logFC,
      values_fill = 0
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()

  significance_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, adj_p_value) %>%
    pivot_wider(
      names_from = contrast,
      values_from = adj_p_value,
      values_fill = Inf
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()

  significance_matrix <- ifelse(significance_matrix < 0.05, "*", "")

  # Optional row ordering
  if (!is.null(ordered_rows)) {
    logFC_matrix <- logFC_matrix[ordered_rows, , drop = FALSE]
    significance_matrix <- significance_matrix[ordered_rows, , drop = FALSE]
  }

  # Step 3: Annotations
  col_annotations <- data.frame(
    contrast = contrasts,
    group = ifelse(grepl("^group_timepointADUEndur", contrasts), "EE", "RE"),
    timeline = "4hrPost"
  )

  col_annotations$timeline <- factor(col_annotations$timeline, levels = c("4hrPost"))

  col_annotation <- HeatmapAnnotation(
    Exercise = col_annotations$group,
    Timeline = col_annotations$timeline,
    col = list(
      Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77"),
      Timeline = c("4hrPost" = "#882E72")
    ),
    annotation_legend_param = list(
      Exercise = list(title = "Modality"),
      Timeline = list(title = "Timepoint")
    ),
    show_annotation_name = FALSE
  )

  # Step 4: Heatmap
  heatmap <- Heatmap(
    logFC_matrix,
    name = "logFC",
    col = colorRamp2(
      c(-max(abs(logFC_matrix)), 0, max(abs(logFC_matrix))),
      c("blue", "white", "red")
    ),
    cluster_rows = FALSE, # Turn to FALSE for wgcna main figures
    cluster_columns = FALSE,
    show_row_names = FALSE, # Turn to FALSE for wgcna main figures
    show_column_names = FALSE,
    border = TRUE,
    column_split = col_annotations$group,
    column_title = NULL,
    top_annotation = col_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (significance_matrix[i, j] == "*") {
        grid.text("*", x, y, gp = gpar(fontsize = 11))
      }
    }
  )

  return(heatmap)
}
ordered_rownames <- rownames(draw(ph10_ir)@ht_list[[1]]@matrix)

contrasts_to_include2 <- c(
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise"
)
# Step 2: Recreate the second heatmap using that order
ph10_heatmap_all_ordered <- create_4hr_heatmap(
  data = precawg_phos_da,
  contrasts = contrasts_to_include2,
  features = ordered_rownames,
  gene_col = "gene_symbol_with_phosphosite",
  ordered_rows = ordered_rownames  # <- row order override
)

#pdf("ph10_110725.pdf", width = 5, height = 4)
draw(ph10_ir + ph10_heatmap_all_ordered) # Figure 5C (right)
#dev.off()

## Pr5 with BCAA
precawg_prot_da_sig_ee_v2 <- precawg_prot_da %>% # feature ids of EE protein DA
  filter(contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise") %>%
  filter(adj_p_value<0.05)%>%
  pull(feature_id)
precawg_prot_da_sig_re_v2 <- precawg_prot_da %>% # feature ids of RE protein DA
  filter(contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise") %>%
  filter(adj_p_value<0.05)%>%
  pull(feature_id)

# module connectivity for protein wgcna
kME_pr = signedKME(datExpr1, MEs1)
kME_pr <- rownames_to_column(kME_pr, var = "feature_id")

mm_pr <- mm_pr %>%
  mutate(
    da = case_when(
      feature_id %in% precawg_prot_da_sig_ee_v2 &
        feature_id %in% precawg_prot_da_sig_re_v2 ~ "EE & RE",
      feature_id %in% precawg_prot_da_sig_ee_v2 ~ "EE",
      feature_id %in% precawg_prot_da_sig_re_v2 ~ "RE",
      TRUE ~ NA_character_  # Default value for unmatched cases
    )
  )
pr5_da <- mm_pr %>%
  filter(module == "Pr5", !is.na(da)) %>%
  # Join with kME_ph to include kME10 values
  left_join(kME_pr %>% select(feature_id, kME5), by = "feature_id") %>%
  # Separate into positive and negative lists
  mutate(kME_category = ifelse(kME5 > 0, "positive", "negative"))

pr5_features <- mm_pr %>%
  filter(module == "Pr5") %>%
  pull(feature_id)
pr5_genes <- kME_pr %>%
  filter(feature_id %in% pr5_features) %>%
  left_join(HUMAN_FEATURE_TO_GENE %>% select(feature_id, gene_symbol), by = "feature_id") %>%
  mutate(DA_status = ifelse(feature_id %in% pr5_da$feature_id, "DA", "non-DA")) %>%
  arrange(DA_status == "DA")
pr5_genes$feature_id <- factor(pr5_genes$feature_id, levels = pr5_genes$feature_id[order(pr5_genes$kME5)])

# Figure 5D (left)
#pdf("Pr5_scatter.pdf", width = 5, height = 4)
ggplot(pr5_genes, aes(x = feature_id, y = kME5)) +
  geom_point(aes(color = DA_status), size = 3, alpha = 0.7) +
  geom_label_repel(
    data = pr5_genes %>% filter(DA_status == "DA"),
    aes(label = gene_symbol),
    size = 4,
    nudge_x = 0.4,
    color = "black",
    fill = "white",
    box.padding = 0.2,
    point.padding = 0.3,
    min.segment.length = 0.1,
    max.overlaps = 15,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  labs(
    title = "kME for Pr5 proteins",
    x = NULL,
    y = "kME",
    color = "Feature Status"
  ) +
  scale_color_manual(values = c("DA" = "red", "non-DA" = "grey")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#dev.off()

# correlation between Pr5 DA and BCAA
common_ids <- intersect(colnames(prot_pre), rownames(cli_full))

# Subset data for common IDs
prot_pre_feat <- prot_pre[pr5_da$feature_id, common_ids, drop = FALSE] #P15104_S343 is GLUL-S343. Q9UHB6_S343s is LIMA1-S343
clinic_feat <- cli_full[common_ids, c("Isoleucine", "Leucine", "Valine")]
phos_pre_feat <- as.matrix(prot_pre_feat)  # Convert to matrix
clinic_feat <- as.matrix(clinic_feat)  # Convert to matrix

# Compute bicorrelation and p-values
cor_results <- corAndPvalue(t(prot_pre_feat), clinic_feat, use = "pairwise.complete.obs", method = "spearman")
cor_matrix <- cor_results$cor
pval_matrix <- cor_results$p

ordered_indices <- order(cor_matrix[, "Isoleucine"], decreasing = FALSE)
cor_matrix <- cor_matrix[ordered_indices, ]

# Apply the same order to pval_matrix
pval_matrix <- pval_matrix[ordered_indices, ]

# Convert p-values to significance markers
sig_markers <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
sig_markers[pval_matrix < 0.05] <- "*"
sig_markers[pval_matrix < 0.01] <- "**"
sig_markers[pval_matrix < 0.001] <- "***"

rownames(cor_matrix) <- pr5_da$gene_symbol[match(rownames(cor_matrix), pr5_da$feature_id)]

range_vals <- range(cor_matrix, na.rm = TRUE)
col_fun <- circlize::colorRamp2(seq(range_vals[1], range_vals[2], length.out = 50), color_palette)

pr5_bcaa <- Heatmap(
  cor_matrix,
  col = col_fun,
  name = "Spearman",
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_names_side = "left",
  cluster_columns = FALSE,
  show_heatmap_legend = TRUE,
  row_dend_side = "left",
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(sig_markers[i, j]) && sig_markers[i, j] != "") {
      # Convert fill to RGB
      rgb_vals <- col2rgb(fill) / 255
      luminance <- sum(rgb_vals * c(0.299, 0.587, 0.114))  # perceptual brightness

      text_color <- ifelse(luminance < 0.5, "white", "black")

      grid.text(
        sig_markers[i, j],
        x = x, y = y,
        gp = gpar(fontsize = 11, col = text_color)
      )
    }
  },
  column_names_gp = gpar(fontsize = 10)
)

ordered_rownames_pr5 <- rownames(draw(pr5_bcaa)@ht_list[[1]]@matrix)

pr5_heatmap_all_ordered <- create_4hr_heatmap(
  data = precawg_prot_da,
  contrasts = contrasts_to_include2,
  features = rownames(cor_matrix),
  gene_col = "gene_symbol",
  ordered_rows = ordered_rownames_pr5
)
draw(pr5_bcaa + pr5_heatmap_all_ordered) # Figure 5D (right)

#### Figure S5: Supplementary to WGCNA
# Figure S5A: Intercorrelation of module eigengene at baseline
common_subjects <- Reduce(intersect, list(rownames(MEs), rownames(MEs1), rownames(MEs2), rownames(MEs3)))
alltog_t2 <- MEs[common_subjects, ]
alltog_pr2 <- MEs1[common_subjects, ]
alltog_ph2 <- MEs2[common_subjects, ]
alltog_m2 <- MEs3[common_subjects, ]


# Step 3: Verify that all datasets now have the same subjects
identical(rownames(alltog_t2), rownames(alltog_pr2)) &&
  identical(rownames(alltog_t2), rownames(alltog_ph2)) &&
  identical(rownames(alltog_t2), rownames(alltog_m2))
#TRUE

# Correlation
cor.mtest <- function(mat, method = "pearson") {
  n <- ncol(mat)
  p_mat <- matrix(NA, n, n)
  colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = method)
      p_mat[i, j] <- p_mat[j, i] <- tmp$p.value
    }
  }
  diag(p_mat) <- 0
  return(list(p = p_mat))
}
# Combine into one big cross-omic ME matrix
all_MEs <- cbind(alltog_t2, alltog_pr2, alltog_ph2, alltog_m2)

# Correlation matrix
ME_corr <- cor(all_MEs, method = "pearson", use = "pairwise.complete.obs")

# Optional: compute p-values as well
ME_pval <- cor.mtest(all_MEs, method = "pearson")$p

module_prefix <- gsub("[0-9]+$", "", colnames(ME_corr))  # Strip numeric suffix
module_colors <- c(
  "T" = "#377EB8",     # Transcript
  "Pr" = "#228833",    # Protein
  "Ph" = "#F3A02B",    # Phospho
  "M" = "#6D4B08"      # Metabolite
)
get_prefix <- function(x) {
  sub("[0-9]+$", "", x)  # e.g., "Ph13" → "Ph"
}

# Assign color to each module
module_prefixes <- sapply(rownames(ME_corr), get_prefix)
row_module_colors <- module_colors[module_prefixes]

module_prefixes_col <- sapply(colnames(ME_corr), get_prefix)
col_module_colors <- module_colors[module_prefixes_col]

col_fun <- colorRamp2(
  c(-1, 0, 1),
  rev(brewer.pal(11, "RdBu")[c(1, 6, 11)])  # RdBu is reversed to have red = high, blue = low
)
# Row annotation
row_ha <- rowAnnotation(
  Omic = anno_simple(
    sapply(rownames(ME_corr), get_prefix),
    col = module_colors
  ),
  annotation_legend_param = list(title = "Omics Type")
)

# Column annotation
omic_levels <- c("T", "Pr", "Ph", "M")
col_omics <- factor(gsub("[0-9]+", "", colnames(ME_corr)), levels = omic_levels)
col_ha <- HeatmapAnnotation(
  Omic = col_omics,
  col = list(Omic = module_colors),
  annotation_legend_param = list(title = "Assay"),
  show_annotation_name = FALSE
)
ht <- Heatmap(
  ME_corr,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  top_annotation = col_ha,
  left_annotation = row_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(ME_pval[i, j]) && ME_pval[i, j] < 0.05) {
      grid.text("*", x = x, y = y, gp = gpar(fontsize = 12, col = "black"))
    }
  }
)

#pdf("cross_module_heatmap.pdf", width = 10, height = 8)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right") # Figure S5A
#dev.off()

## FIgure S5B: Correlation between deconvolution outcome and module eigengene at baseline.
# 1. Trans.
cor_decon_t = bicorAndPvalue(precawg_decon, MEs) # alternative biweight midcor2
df_decon_t = melt(cor_decon_t$bicor) %>% dplyr::rename(bicor = value)
df_decon_t$pval = melt(cor_decon_t$p)$value
df_decon_t$obs = melt(cor_decon_t$nObs)$value
df_decon_t <- df_decon_t %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

deconXmod0 = dcast(df_decon_t, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(deconXmod0) = deconXmod0$Var1
deconXmod0$Var1 = NULL

# set up signif overlay
sig_table_decon_t = dcast(df_decon_t, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_decon_t) = sig_table_decon_t$Var1
sig_table_decon_t$Var1 = NULL

# 2. Prot.
common_samples <- intersect(rownames(precawg_decon), rownames(MEs1))

# Subset both data frames to keep only the common samples and ensure they are in the same order
decon_pre_prot <- precawg_decon[common_samples, ]
cor_decon_pr = bicorAndPvalue(decon_pre_prot, MEs1) # alternative biweight midcor
df_decon_pr = melt(cor_decon_pr$bicor) %>% dplyr::rename(bicor = value)
df_decon_pr$pval = melt(cor_decon_pr$p)$value
df_decon_pr$obs = melt(cor_decon_pr$nObs)$value
df_decon_pr <- df_decon_pr %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

# set up correlation data frame for heatmap
deconXmod1 = dcast(df_decon_pr, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(deconXmod1) = deconXmod1$Var1
deconXmod1$Var1 = NULL

# set up signif overlay
sig_table_decon_pr = dcast(df_decon_pr, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_decon_pr) = sig_table_decon_pr$Var1
sig_table_decon_pr$Var1 = NULL

# 3. Phosphoproteomics.
common_samples <- intersect(rownames(precawg_decon), rownames(MEs2))

# Subset both data frames to keep only the common samples and ensure they are in the same order
decon_pre_ph <- precawg_decon[common_samples, ]

cor_decon_ph = bicorAndPvalue(decon_pre_ph, MEs2) # alternative biweight midcor
df_decon_ph = melt(cor_decon_ph$bicor) %>% dplyr::rename(bicor = value)
df_decon_ph$pval = melt(cor_decon_ph$p)$value
df_decon_ph$obs = melt(cor_decon_ph$nObs)$value
df_decon_ph <- df_decon_ph %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

# set up correlation data frame for heatmap
deconXmod2 = dcast(df_decon_ph, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(deconXmod2) = deconXmod2$Var1
deconXmod2$Var1 = NULL

# set up signif overlay
sig_table_decon_ph = dcast(df_decon_ph, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_decon_ph) = sig_table_decon_ph$Var1
sig_table_decon_ph$Var1 = NULL


# 4. Metabolomics
common_samples <- intersect(rownames(precawg_decon), rownames(MEs3))

# Subset both data frames to keep only the common samples and ensure they are in the same order
decon_pre_m <- precawg_decon[common_samples, ]
cor_decon_m = bicorAndPvalue(decon_pre_m, MEs3) # alternative biweight midcor
df_decon_m = melt(cor_decon_m$bicor) %>% dplyr::rename(bicor = value)
df_decon_m$pval = melt(cor_decon_m$p)$value
df_decon_m$obs = melt(cor_decon_m$nObs)$value
df_decon_m <- df_decon_m %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))

# set up correlation data frame for heatmap
deconXmod3 = dcast(df_decon_m, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(deconXmod3) = deconXmod3$Var1
deconXmod3$Var1 = NULL

# set up signif overlay
sig_table_decon_m = dcast(df_decon_m, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table_decon_m) = sig_table_decon_m$Var1
sig_table_decon_m$Var1 = NULL

## Combine decon x ME correlations
deconXmod_combined <- cbind(deconXmod0, deconXmod1, deconXmod2, deconXmod3)
deccon_sig_table_combined <- cbind(sig_table_decon_t, sig_table_decon_pr, sig_table_decon_ph, sig_table_decon_m)
#heatmap
deconXmod_combined <- t(deconXmod_combined)
deccon_sig_table_combined <- t(deccon_sig_table_combined)

deconXmod_combined <- as.data.frame(deconXmod_combined) %>%
  rownames_to_column(var = "Rowname") %>%            # Move rownames to a column
  mutate(Assay = case_when(
    grepl("^T", Rowname)  ~ "Transcriptomics",
    grepl("^Pr", Rowname) ~ "Proteomics",
    grepl("^Ph", Rowname) ~ "Phosphoproteomics",
    grepl("^M", Rowname)  ~ "Metabolomics"
  )) %>%
  column_to_rownames(var = "Rowname")

correct_order <- c(
  paste0("T", 1:12),   # T1 - T12
  paste0("Pr", 1:14),  # Pr1 - Pr14
  paste0("Ph", 1:18),  # Ph1 - Ph18
  paste0("M", 1:7)     # M1 - M7
)

# Convert row names to factors with correct order
deconXmod_combined <- deconXmod_combined[order(factor(rownames(deconXmod_combined), levels = correct_order)), ]
deccon_sig_table_combined <- deccon_sig_table_combined[order(factor(rownames(deccon_sig_table_combined), levels = correct_order)), ]

# Reset row names to maintain ordering
rownames(deconXmod_combined) <- correct_order[correct_order %in% rownames(deconXmod_combined)]
rownames(deccon_sig_table_combined) <- correct_order[correct_order %in% rownames(deccon_sig_table_combined)]

assay_colors <- setNames(
  c("#4c81b7", "#6ca767", "#d4884a", "#b34f44"),  # Colors for each assay
  c("Transcriptomics", "Proteomics", "Phosphoproteomics", "Metabolomics")  # Assay names
)
row_annotation <- rowAnnotation(
  Assay = deconXmod_combined$Assay,
  col = list(
    Assay = assay_colors)
)

#pdf(file = 'wgcnaXdecon_complexheat_110725.pdf', width = 4, height = 12)
Heatmap(
  as.matrix(deconXmod_combined[, 1:(ncol(deconXmod_combined) - 1)]),
  name = "Correlation",
  col = col_fun,
  right_annotation = row_annotation,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, col = "black"),
  column_names_side = "top",
  show_heatmap_legend = TRUE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    sig_marker <- deccon_sig_table_combined[i, j]

    if (!is.na(sig_marker) && sig_marker != "") {
      rgb_col <- col2rgb(fill) / 255
      luminance <- 0.299 * rgb_col[1, ] + 0.587 * rgb_col[2, ] + 0.114 * rgb_col[3, ]
      text_color <- ifelse(luminance < 0.5, "white", "black")  # Threshold for brightness

      grid.text(sig_marker, x, y, gp = gpar(fontsize = 10, col = text_color))
    }
  }
)
#dev.off()

# Figure S5C: Covariate (age, sex, WC)-adjusted correlation between clinical variables and module eigengene
ME_residuals <- lapply(as.data.frame(MEs), function(y) {
  lm(y ~ cli_full_z2_t$Sex + cli_full_z2_t$Age + cli_full_z2_t$WC)$residuals
})
traits_residuals <- lapply(as.data.frame(cli_full_z2_t), function(y) {
  lm(y ~ cli_full_z2_t$Sex + cli_full_z2_t$Age + cli_full_z2_t$WC)$residuals
})
ME_residuals_matrix <- do.call(cbind, ME_residuals)
traits_residuals_matrix <- do.call(cbind, traits_residuals)

cor_adj = bicorAndPvalue(traits_residuals_matrix, ME_residuals_matrix)
df = melt(cor_adj$bicor) %>% dplyr::rename(bicor = value)
df$pval = melt(cor_adj$p)$value
df$obs = melt(cor_adj$nObs)$value
df <- df %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
# set up correlation data frame for heatmap
traitXmod = dcast(df, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod) = traitXmod$Var1
traitXmod$Var1 = NULL

# set up signif overlay
sig_table = dcast(df, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table) = sig_table$Var1
sig_table$Var1 = NULL

traitXmod_trans_adj <- traitXmod
sig_table_trans_adj <- sig_table

###Proteomics
ME1_residuals <- lapply(as.data.frame(MEs1), function(y) {
  lm(y ~ cli_full_z2_pr$Sex + cli_full_z2_pr$Age + cli_full_z2_pr$WC)$residuals
})
traits_residuals <- lapply(as.data.frame(cli_full_z2_pr), function(y) {
  lm(y ~ cli_full_z2_pr$Sex + cli_full_z2_pr$Age + cli_full_z2_pr$WC)$residuals
})
ME1_residuals_matrix <- do.call(cbind, ME1_residuals)
traits_residuals_matrix <- do.call(cbind, traits_residuals)

cor_adj = bicorAndPvalue(traits_residuals_matrix, ME1_residuals_matrix) # alternative biweight midcor
df = melt(cor_adj$bicor) %>% dplyr::rename(bicor = value)
df$pval = melt(cor_adj$p)$value
df$obs = melt(cor_adj$nObs)$value
df <- df %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
# set up correlation data frame for heatmap
traitXmod = dcast(df, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod) = traitXmod$Var1
traitXmod$Var1 = NULL

# set up signif overlay
sig_table = dcast(df, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table) = sig_table$Var1
sig_table$Var1 = NULL

traitXmod_prot_adj <- traitXmod
sig_table_prot_adj <- sig_table

# Phosphoproteomics
ME2_residuals <- lapply(as.data.frame(MEs2), function(y) {
  lm(y ~ cli_full_z2_ph$Sex + cli_full_z2_ph$Age + cli_full_z2_ph$WC)$residuals
})
traits_residuals <- lapply(as.data.frame(cli_full_z2_ph), function(y) {
  lm(y ~ cli_full_z2_ph$Sex + cli_full_z2_ph$Age + cli_full_z2_ph$WC)$residuals
})
ME2_residuals_matrix <- do.call(cbind, ME2_residuals)
traits_residuals_matrix <- do.call(cbind, traits_residuals)

cor_adj = bicorAndPvalue(traits_residuals_matrix, ME2_residuals_matrix) # alternative biweight midcor
df = melt(cor_adj$bicor) %>% dplyr::rename(bicor = value)
df$pval = melt(cor_adj$p)$value
df$obs = melt(cor_adj$nObs)$value
df <- df %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
# set up correlation data frame for heatmap
traitXmod = dcast(df, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod) = traitXmod$Var1
traitXmod$Var1 = NULL

# set up signif overlay
sig_table = dcast(df, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table) = sig_table$Var1
sig_table$Var1 = NULL

traitXmod_phos_adj <- traitXmod
sig_table_phos_adj <- sig_table

# Metabolomics
ME3_residuals <- lapply(as.data.frame(MEs3), function(y) {
  lm(y ~ cli_full_z2_m$Sex + cli_full_z2_m$Age + cli_full_z2_m$WC)$residuals
})
traits_residuals <- lapply(as.data.frame(cli_full_z2_m), function(y) {
  lm(y ~ cli_full_z2_m$Sex + cli_full_z2_m$Age + cli_full_z2_m$WC)$residuals
})
ME3_residuals_matrix <- do.call(cbind, ME3_residuals)
traits_residuals_matrix <- do.call(cbind, traits_residuals)

cor_adj = bicorAndPvalue(traits_residuals_matrix, ME3_residuals_matrix) # alternative biweight midcor
df = melt(cor_adj$bicor) %>% dplyr::rename(bicor = value)
df$pval = melt(cor_adj$p)$value
df$obs = melt(cor_adj$nObs)$value
df <- df %>%
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
# set up correlation data frame for heatmap
traitXmod = dcast(df, Var1 ~ Var2, value.var = 'bicor', fun.aggregate = mean)
rownames(traitXmod) = traitXmod$Var1
traitXmod$Var1 = NULL

# set up signif overlay
sig_table = dcast(df, Var1 ~ Var2, value.var = 'sig')
rownames(sig_table) = sig_table$Var1
sig_table$Var1 = NULL

traitXmod_metab_adj <- traitXmod
sig_table_metab_adj <- sig_table

##### Combine and visualize
#Use these for adjusted correlation (sex, age, WC)
traitXmod_adj <- cbind(traitXmod_trans_adj, traitXmod_prot_adj, traitXmod_phos_adj, traitXmod_metab_adj)
traitXmod_adj <- as.data.frame(t(traitXmod_adj))
sig_table_adj <- cbind(sig_table_trans_adj, sig_table_prot_adj, sig_table_phos_adj, sig_table_metab_adj)
sig_table_adj <- as.data.frame(t(sig_table_adj))

new_order <- c(
   "BMI", "SBP", "DBP", "HR",
  "Steps", "Movement", "Energy_exp",
  "VO2max_rel", "VO2max", "O2_Pulse",
  "Grip_strength", "Knee_torque",
  "HOMA_IR", "Adipo_IR",
  "HbA1c", "Insulin", "Glucose", "Lactate", "Glycerol", "NEFA", "Trig", "Cholesterol", "HDL", "LDL", "KET", "Glucagon", "Cortisol",
  "Isoleucine", "Leucine", "Valine"
)
traitXmod_adj <- traitXmod_adj[,new_order]
sig_table_adj <- sig_table_adj[,new_order]

traitXmod_adj <- traitXmod_adj %>%
  rownames_to_column(var = "Rowname") %>%            # Move rownames to a column
  mutate(Assay = case_when(
    grepl("^T", Rowname)  ~ "Transcriptomics",
    grepl("^Pr", Rowname) ~ "Proteomics",
    grepl("^Ph", Rowname) ~ "Phosphoproteomics",
    grepl("^M", Rowname)  ~ "Metabolomics"
  )) %>%
  column_to_rownames(var = "Rowname")                # Return rownames

color_palette <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))
col_fun <- color_palette(10)
assay_colors <- setNames(
  c("#377EB8", "#228833", "#F3A02B", "#6D4B08"),  # Colors for each assay
  c("Transcriptomics", "Proteomics", "Phosphoproteomics", "Metabolomics")  # Assay names
)

x_axis_colors <- c(
  "BMI" = "#FCC737",
  "SBP" = "#FCC737",
  "DBP" = "#FCC737",
  "HR" = "#FCC737",
  "Knee_torque" = "#0A5EB0",
  "Grip_strength" = "#0A5EB0",
  "VO2max" = "#FF8000",
  "VO2max_rel" = "#FF8000",
  "O2_Pulse" = "#FF8000",
  "Steps" = "#BDE8CA",
  "Movement" = "#BDE8CA",
  "Energy_exp" = "#BDE8CA",
  "HOMA_IR" = "#41B3A2",
  "Adipo_IR" = "#41B3A2",
  "HbA1c" = "#F95454",
  "Glucose" = "#F95454",
  "Lactate" = "#F95454",
  "Glycerol" = "#F95454",
  "NEFA" = "#F95454",
  "Trig" = "#F95454",
  "KET" = "#F95454",
  "Insulin" = "#F95454",
  "Glucagon" = "#F95454",
  "Cortisol" = "#F95454",
  "Cholesterol" = "#F95454",
  "HDL" = "#F95454",
  "LDL" = "#F95454",
  "Isoleucine" = "#D7C3F1",
  "Leucine" = "#D7C3F1",
  "Valine" = "#D7C3F1"
)
# Create row annotations
row_annotation <- rowAnnotation(
  Assay = traitXmod_adj$Assay,
  col = list(
    Assay = assay_colors)
)

#column annotations
column_groups <- c(
  rep("Anthropometrics", 4),
  rep("Muscle strength", 2),
  rep("Cardiorespiration", 3),
  rep("Accelerometry", 3),
  rep("IR_index", 2),
  rep("Blood", 13),
  rep("Plasma BCAA", 3)
)
column_annotation <- HeatmapAnnotation(
  Clinical_traits = column_groups,
  col = list(Clinical_traits = structure(
    c("#FCC737", "#0A5EB0", "#FF8000", "#BDE8CA", "#41B3A2", "#F95454", "#D7C3F1"),
    names = c("Anthropometrics", "Muscle strength", "Cardiorespiration",
              "Accelerometry", "IR_index", "Blood", "Plasma BCAA")
  )),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)
# Figure S5C
#pdf(file = 'wgcnaXtrait_complexheat_adj_110725.pdf', width = 12, height = 14)
Heatmap(
  as.matrix(traitXmod_adj[, 1:(ncol(traitXmod_adj) - 1)]),  # Exclude annotation column
  name = "Correlation",
  col = col_fun,
  top_annotation = column_annotation,
  right_annotation = row_annotation,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10, col = "black"),
  column_names_side = "top",
  show_heatmap_legend = TRUE,

  # Dynamically set asterisk color based on fill brightness
  cell_fun = function(j, i, x, y, width, height, fill) {
    sig_marker <- sig_table_adj[i, j]
    if (!is.na(sig_marker) && sig_marker != "") {
      # Convert fill to RGB and compute brightness
      rgb_val <- col2rgb(fill) / 255
      brightness <- 0.299 * rgb_val[1, ] + 0.587 * rgb_val[2, ] + 0.114 * rgb_val[3, ]
      text_col <- ifelse(brightness < 0.5, "white", "black")  # Dark background → white text
      grid.text(sig_marker, x, y, gp = gpar(fontsize = 10, col = text_col))
    }
  }
)
#dev.off()
