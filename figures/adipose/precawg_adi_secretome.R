library(MotrpacHumanPreSuspension)
library(TMSig)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(purrr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(tibble)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(stats)
library(WGCNA)
library(reshape2)
allowWGCNAThreads()
library(foreach)
library(doParallel)
library(readr)
library(forcats)

motrpac_qc <- load_qc()
motrpac_da <- load_differential_analysis()

blood_pr <- unique(motrpac_qc$blood$`prot-ol`$feature_metadata$assay)
blood_metabolites <- motrpac_da$blood %>%
  keep(~ grepl("metab", .x$assay[1])) %>% # Select elements containing "metab" in assay column
  bind_rows() %>%
  pull(unique(feature_id))


hum_sec <- read.table("sa_location_Secreted.tsv", header = TRUE, sep = "\t") # from HPA
hum_sec_geneid <- hum_sec$Gene
hum_sec_int <- union(blood_pr, hum_sec_geneid)

#Transcriptomics
t_vol2 <- t_vol %>% # Obtain all_vol from precawg_adi_da.R
  filter(facet_label %in% c("EE 45min vs Control",
                            "EE 4hr vs Control",
                            "EE 24hr vs Control",
                            "RE 45min vs Control",
                            "RE 4hr vs Control",
                            "EE 24hr vs Control")) %>%
  filter(adj_p_value<0.05)
#prot
pr_vol2 <- pr_vol %>%
  filter(facet_label %in% c("EE 4hr vs Control",
                            "RE 4hr vs Control")) %>%
  filter(adj_p_value<0.05)

combined_dat <- bind_rows(
  t_vol2 %>% mutate(uniprot = NA),  # Add `uniprot` column with NA to match `all_volcano_prot2`
  pr_vol2
)

plot_data <- combined_dat %>%
  mutate(facet_label = gsub("vs Control", "", facet_label),
         facet_label = trimws(facet_label),  # Clean up extra spaces
         is_hum_sec = ifelse(gene_symbol %in% hum_sec_int, "Blood protein", "Other")) %>%
  group_by(facet_label, is_hum_sec) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(facet_label = c("EE 45min", "EE 4hr", "EE 24hr", "RE 45min", "RE 4hr", "RE 24hr"),
           is_hum_sec = c("Blood protein", "Other"),
           fill = list(count = 0))  # Ensure all combinations exist

# Set the desired order of facet labels
plot_data <- plot_data %>%
  mutate(facet_label = factor(facet_label, levels = c("EE 45min", "EE 4hr", "EE 24hr", 
                                                      "RE 45min", "RE 4hr", "RE 24hr")),
         group = ifelse(grepl("^EE", facet_label), "EE", "RE"))  # Add EE/RE group

blood_protein_genes <- combined_dat %>%
  mutate(facet_label = gsub("vs Control", "", facet_label),
         facet_label = trimws(facet_label),  # Clean up extra spaces
         is_hum_sec = ifelse(gene_symbol %in% hum_sec_int, "Blood protein", "Other")) %>%
  filter(is_hum_sec == "Blood protein") %>%
  dplyr::select(assay, facet_label, gene_symbol) %>%
  distinct() 
length(unique(blood_protein_genes$gene_symbol)) # 60 secretome candidates

# Figure 6B
ggplot(plot_data, aes(x = facet_label, y = count, fill = is_hum_sec)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(
    values = c("Blood protein" = "purple", "Other" = "gray"),
    labels = c("Blood protein" = "Yes", "Other" = "No")  # Update legend labels
  ) +
  labs(
    y = "DA transcripts + proteins",
    fill = "Detected in blood"  # Update legend title
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(1, 1, 2, 1, "cm")
  ) +
  geom_text(
    aes(label = ifelse(is_hum_sec == "Blood protein" & count > 0, count, "")), 
    position = position_stack(vjust = 0.5), size = 5, color = "white"
  ) +
  # Add EE/RE color boxes below the x-axis
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -5, ymax = -30, alpha = 1, fill = "#d95f02") +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -5, ymax = -30, alpha = 1, fill = "#1b9e77") +
  annotate("text", x = 2, y = -8.5, label = "", color = "white", size = 4, fontface = "bold") +
  coord_cartesian(clip = "off") +  # Allow elements outside the plot area
  # Separate legend for EE and RE
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(title = "Group", override.aes = list(fill = c("#d95f02", "#1b9e77")))
  ) +
  scale_color_manual(values = c(EE = "#d95f02", RE = "#1b9e77"))


m_vol2 <- m_vol %>%
  filter(facet_label %in% c("EE 45min vs Control",
                            "EE 4hr vs Control",
                            "EE 24hr vs Control",
                            "RE 45min vs Control",
                            "RE 4hr vs Control",
                            "EE 24hr vs Control")) %>%
  filter(adj_p_value<0.05)

plot_data_metab <- m_vol2 %>%
  mutate(
    facet_label = gsub("vs Control", "", facet_label),
    facet_label = trimws(facet_label),  # Clean up extra spaces
    is_blood_sec = ifelse(feature_id %in% blood_metabolites, "Human Secretome", "Other")  # Check against blood_metabolites
  ) %>%
  group_by(facet_label, is_blood_sec) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(
    facet_label = c("EE 45min", "EE 4hr", "EE 24hr", "RE 45min", "RE 4hr", "RE 24hr"),
    is_blood_sec = c("Human Secretome", "Other"),
    fill = list(count = 0)
  )

# Set the desired order of facet labels
plot_data_metab <- plot_data_metab %>%
  mutate(facet_label = factor(facet_label, levels = c("EE 45min", "EE 4hr", "EE 24hr", 
                                                      "RE 45min", "RE 4hr", "RE 24hr")),
         group = ifelse(grepl("^EE", facet_label), "EE", "RE"))  # Add EE/RE group
# Figure 6C
ggplot(plot_data_metab, aes(x = facet_label, y = count, fill = is_blood_sec)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(
    values = c("Human Secretome" = "purple", "Other" = "gray"),
    labels = c("Human Secretome" = "Yes", "Other" = "No")  # Update legend labels
  ) +
  labs(
    y = "DA metabolites",
    fill = "Detected in blood"  # Update legend title
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(1, 1, 2, 1, "cm")
  ) +
  geom_text(
    aes(label = ifelse(is_blood_sec == "Human Secretome" & count > 0, count, "")), 
    position = position_stack(vjust = 0.5), size = 5, color = "white"
  ) +
  # Add EE/RE color boxes below the x-axis
  annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -3, ymax = -11, alpha = 1, fill = "#d95f02") +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -3, ymax = -11, alpha = 1, fill = "#1b9e77") +
  annotate("text", x = 2, y = -8.5, label = "", color = "white", size = 4, fontface = "bold") +
  coord_cartesian(clip = "off") +  # Allow elements outside the plot area
  # Separate legend for EE and RE
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(title = "Group", override.aes = list(fill = c("#d95f02", "#1b9e77")))
  ) +
  scale_color_manual(values = c(EE = "#d95f02", RE = "#1b9e77"))

### Figure S6A: plot secretome DAs. 
# Custom functions for DA heatmap
create_feature_heatmap <- function(data, contrasts, order, features = NULL, gene_col) {
  # Step 1: Filter data for the specified contrasts
  heatmap_data <- data %>%
    filter(contrast %in% contrasts) %>%
    mutate(
      group = ifelse(grepl("^group_timepointADUEndur", contrast), "EE", "RE"),
      timeline = case_when(
        grepl("post_15_30_45_min", contrast) ~ "Post 45 min",
        grepl("post_3.5_4_hr", contrast) ~ "Post 4 hr",
        grepl("post_24_hr", contrast) ~ "Post 24 hr"
      ),
      contrast = factor(contrast, levels = order) # Ensure correct contrast order
    ) %>%
    arrange(contrast)
  
  # If features are provided, filter data to include only those features
  if (!is.null(features)) {
    heatmap_data <- heatmap_data %>%
      filter(!!sym(gene_col) %in% features)
  }
  
  # Step 2: Pivot data to wide format
  logFC_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, logFC) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = logFC,
      values_fill = 0
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  # Create a significance mask for asterisks
  significance_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, adj_p_value) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = adj_p_value,
      values_fill = Inf
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  significance_matrix <- ifelse(significance_matrix < 0.05, "*", "")
  
  # Step 3: Column annotations
  col_annotations <- data.frame(
    contrast = order,
    group = ifelse(grepl("^group_timepointADUEndur", order), "EE", "RE"),
    timeline = rep(c("Post 45 min", "Post 4 hr", "Post 24 hr"), times = 2) # Repeat for EE and RE
  )
  
  # Set order for timeline legend
  col_annotations$timeline <- factor(col_annotations$timeline, levels = c("Post 45 min", "Post 4 hr", "Post 24 hr"))
  
  col_annotation <- HeatmapAnnotation(
    Exercise = col_annotations$group,
    Timeline = col_annotations$timeline,
    col = list(
      Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77"),
      Timeline = c("Post 45 min" = "#AE76A3", "Post 4 hr" = "#882E72", "Post 24 hr" = "#61194F")
    ),
    annotation_legend_param = list(
      Exercise = list(title = "Modality"),
      Timeline = list(title = "Timepoint", at = c("Post 45 min", "Post 4 hr", "Post 24 hr"))
    ),
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm"),
    height=unit(0.5, "mm")
  )
  
  # Step 4: Heatmap creation
  heatmap <- Heatmap(
    logFC_matrix,
    name = "logFC",
    col = colorRamp2(c(-max(abs(logFC_matrix)), 0, max(abs(logFC_matrix))), c("blue", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_dend_side = "right", # Move row dendrogram to the right
    row_names_side = "left", # Move row labels (gene names) to the left
    column_split = col_annotations$group,
    column_title = NULL,
    top_annotation = col_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (significance_matrix[i, j] == "*") {
        # Convert fill (RGB) to luminance
        rgb_vals <- col2rgb(fill)
        luminance <- (0.299 * rgb_vals[1, ] + 0.587 * rgb_vals[2, ] + 0.114 * rgb_vals[3, ]) / 255
        
        # Decide text color based on brightness
        text_color <- ifelse(luminance < 0.5, "white", "black")
        
        grid.text(
          significance_matrix[i, j],
          x,
          y,
          gp = gpar(fontsize = 20, col = text_color)
        )
      }
    }
  )
  
  return(heatmap)
}
# For protein data
create_4hr_heatmap <- function(data, contrasts, features = NULL, gene_col) {
  
  # Step 1: Filter data for the specified contrasts
  heatmap_data <- data %>%
    filter(contrast %in% contrasts) %>%
    mutate(
      group = ifelse(grepl("^group_timepointADUEndur", contrast), "EE", "RE"),
      timeline = "4hrPost",
      contrast = factor(contrast, levels = contrasts)
    ) %>%
    arrange(contrast)
  
  # If features are provided, filter to those features
  if (!is.null(features)) {
    heatmap_data <- heatmap_data %>%
      filter(!!sym(gene_col) %in% features)
  }
  
  # STEP to solve duplicate rows: keep the one with smallest adj_p_value
  heatmap_data <- heatmap_data %>%
    group_by(!!sym(gene_col), contrast) %>%
    slice_min(adj_p_value, with_ties = FALSE) %>%
    ungroup()
  
  # Step 2: Pivot to wide format (logFC matrix)
  logFC_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, logFC) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = logFC,
      values_fill = 0
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  # Step 3: Create significance matrix for asterisks
  significance_matrix <- heatmap_data %>%
    select(all_of(gene_col), contrast, adj_p_value) %>%
    tidyr::pivot_wider(
      names_from = contrast,
      values_from = adj_p_value,
      values_fill = Inf
    ) %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  significance_matrix <- apply(significance_matrix, c(1, 2), function(p) {
    if (p < 0.05) {
      "*"
    } else if (p < 0.1) {
      "\u2020"  # Unicode for dagger (†)
    } else {
      ""
    }
  })
  
  # Step 4: Column annotations
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
    show_annotation_name = FALSE,
    simple_anno_size = unit(2, "mm"),
    height=unit(0.5, "mm")
  )
  
  # Step 5: Heatmap
  heatmap <- Heatmap(
    logFC_matrix,
    name = "logFC",
    col = colorRamp2(
      c(-max(abs(logFC_matrix)), 0, max(abs(logFC_matrix))),
      c("blue", "white", "red")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_dend_side = "right",
    row_names_side = "left",
    column_split = col_annotations$group,
    column_title = NULL,
    top_annotation = col_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      symbol <- significance_matrix[i, j]
      if (symbol != "") {
        font_color <- if (symbol == "*") "white" else "black"
        font_size  <- if (symbol == "*") 20 else 18
        grid.text(
          label = symbol,
          x = x,
          y = y,
          gp = gpar(fontsize = font_size, col = font_color)
        )
      }
    }
  )
  
  return(heatmap)
}
contrasts_to_include <- c(
  "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
)
contrast_order <- c(
  "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
)
contrasts_to_include2 <- c(
  "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
  "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise"
)
sec_trans_da <- t_vol2 %>%
  filter(gene_symbol %in% blood_protein_genes$gene_symbol) %>%
  pull(gene_symbol)

# Figure S6A (top)
sec_trans_da_plot <- create_feature_heatmap(
  data = precawg_trans_da,
  contrasts = contrasts_to_include,
  order = contrast_order,
  features = sec_trans_da,
  gene_col = "gene_symbol" # or "feature_id"
)

sec_prot_da <- pr_vol2 %>%
  filter(gene_symbol %in% blood_protein_genes$gene_symbol) %>%
  pull(gene_symbol)

# Figure S6A (bottom)
sec_prot_da_plot <- create_4hr_heatmap(
  data = precawg_prot_da,
  contrasts = contrasts_to_include2,
  features = sec_prot_da,
  gene_col = "gene_symbol" # or "feature_id"
)

## Correlation of feature expression between adipose (60 candidates from transcript and protein) and blood protein at baseline
#plasma prot qc
blood_prot <- motrpac_qc$blood$`prot-ol`$qc_norm
blood_prot_meta <- motrpac_qc$blood$`prot-ol`$sample_metadata
blood_prot_feat <- motrpac_qc$blood$`prot-ol`$feature_metadata
# Ensure vialLabel is character
vial_to_pid <- motrpac_qc$blood$`prot-ol`$sample_metadata %>%
  filter(Timepoint == "pre_exercise", visitcode == "ADU_BAS") %>%
  dplyr::select(vialLabel, pid) %>%
  distinct() %>%
  mutate(vialLabel = as.character(vialLabel))  # Convert vialLabel to character

# Filter blood_prot for pre_exercise vialLabels
blood_prot_pre <- blood_prot %>%
  as.data.frame() %>%
  dplyr::select(all_of(vial_to_pid$vialLabel))  # Now both are character
colnames(blood_prot_pre) <- vial_to_pid$pid[match(colnames(blood_prot_pre), vial_to_pid$vialLabel)]

# Handle duplicates by averaging before setting rownames
blood_prot_pre <- blood_prot_pre %>%
  rownames_to_column(var = "feature_id") %>%  # Move rownames to a column
  left_join(blood_prot_feat %>% select(feature_id, assay), by = "feature_id") %>%  # Map feature_id to assay
  group_by(assay) %>%  # Group by assay to handle duplicates
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%  # Average only numeric columns
  ungroup() %>%  # Remove grouping
  column_to_rownames(var = "assay")  # Set assay as the new rownames


blood_prot_pre <- t(blood_prot_pre) #blood protein baseline ready

#Trans qc baseline
adi_t <- motrpac_qc$adipose$`transcript-rna-seq`$qc_norm
vial_to_pid <- motrpac_qc$adipose$`transcript-rna-seq`$sample_metadata %>%
  filter(Timepoint == "pre_exercise", visitcode == "ADU_BAS") %>%
  dplyr::select(vialLabel, pid) %>%
  distinct() %>%
  mutate(vialLabel = as.character(vialLabel))
adi_t_pre <- adi_t %>%
  as.data.frame() %>%
  dplyr::select(all_of(vial_to_pid$vialLabel))  # Now both are character
colnames(adi_t_pre) <- vial_to_pid$pid[match(colnames(adi_t_pre), vial_to_pid$vialLabel)]
adi_t_pre_with_symbol <- adi_t_pre %>%
  as.data.frame() %>%  # Ensure it's a data frame
  rownames_to_column(var = "feature_id") %>%  # Convert rownames to a column
  left_join(HUMAN_FEATURE_TO_GENE %>% select(feature_id, gene_symbol), 
            by = "feature_id") %>%  # Map ENSG IDs to gene symbols
  filter(!is.na(gene_symbol)) %>%  # Remove rows without matching gene symbols
  group_by(gene_symbol) %>%  # Group to handle duplicate gene symbols
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%  # Take max for duplicates
  column_to_rownames(var = "gene_symbol") # adipose transcriptome baseline ready. 
adi_t_pre_with_symbol <- t(adi_t_pre_with_symbol)


adi_subset <- adi_t_pre_with_symbol[, colnames(adi_t_pre_with_symbol) %in% colnames(blood_prot_pre)]

common_rows <- intersect(rownames(blood_prot_pre), rownames(adi_subset))

blood_prot_pre <- blood_prot_pre[common_rows, , drop = FALSE]
adi_subset <- adi_subset[common_rows, , drop = FALSE]

trans_sec_cor <- bicorAndPvalue(blood_prot_pre, adi_subset, use = "pairwise.complete.obs")
df = melt(trans_sec_cor$bicor) %>% dplyr::rename(bicor = value)
df$pval = melt(trans_sec_cor$p)$value
df$obs = melt(trans_sec_cor$nObs)$value
df <- df %>% 
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
df_filtered <- df %>% 
  filter(as.character(Var1) == as.character(Var2)) %>%
  filter(Var1 %in% combined_dat$gene_symbol)

# Figure S6B
adi_trans_blood_prot_plot <- ggplot(df_filtered, aes(x = reorder(Var1, -bicor), y = bicor)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "purple") +  # Purple bars with black borders
  geom_text(aes(label = sig), vjust = -0.5, size = 5, color = "black") +  # Significance on top
  labs(y = "Bicor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        legend.position = "none")  # Remove legend

## Metab. 
metab_sec <- m_vol2 %>%
  mutate(
    facet_label = gsub("vs Control", "", facet_label),
    facet_label = trimws(facet_label),  # Clean up extra spaces
    is_blood_sec = ifelse(feature_id %in% blood_metabolites, "Detected in blood", "Other")  # Check against blood_metabolites
  )
metab_sec_da <- metab_sec %>%
  filter(is_blood_sec == "Detected in blood")
metab_sec_da_list <- metab_sec_da$feature_id

blood_metab_qc <- list()

# Loop through all objects starting with "metab"
for (metab_name in names(motrpac_qc$blood)) {
  if (startsWith(metab_name, "metab")) {
    # Extract qc_norm and sample_metadata for the current metab object
    qc_norm <- motrpac_qc$blood[[metab_name]]$qc_norm
    sample_metadata <- motrpac_qc$blood[[metab_name]]$sample_metadata
    
    # Filter sample_metadata to keep only pre-exercise samples with visitcode "ADU_BAS"
    pre_exercise_samples <- sample_metadata$Timepoint == "pre_exercise" & sample_metadata$visitcode == "ADU_BAS"
    filtered_metadata <- sample_metadata[pre_exercise_samples, ]
    
    # Subset qc_norm to include only pre-exercise vialLabels
    pre_exercise_qc_norm <- qc_norm[, colnames(qc_norm) %in% filtered_metadata$vialLabel, drop = FALSE]
    
    # Rename vialLabels (columns) to their corresponding pid from sample_metadata
    colnames(pre_exercise_qc_norm) <- filtered_metadata$pid[match(colnames(pre_exercise_qc_norm), filtered_metadata$vialLabel)]
    
    # Add the filtered qc_norm to the list
    blood_metab_qc[[metab_name]] <- pre_exercise_qc_norm
  }
}

# Combine all dataframes by rows, allowing for differing columns
all_columns <- unique(unlist(lapply(blood_metab_qc, colnames)))
blood_metab_qc <- lapply(blood_metab_qc, function(df) {
  missing_columns <- setdiff(all_columns, colnames(df))
  df[, missing_columns] <- NA  # Add missing columns filled with NA
  df <- df[, all_columns]  # Reorder columns to match the full set
  return(df)
})

# Combine all dataframes by rows
blood_metab_qc_pre <- do.call(rbind, blood_metab_qc)

# View the final dataframe
cleaned_rownames <- sub("^[^.]+\\.", "", rownames(blood_metab_qc_pre))
blood_metab_qc_pre$feature_id <- cleaned_rownames

blood_metab_qc_pre <- blood_metab_qc_pre %>%
  filter(feature_id %in% metab_sec_da_list)
blood_metab_qc_pre <- blood_metab_qc_pre %>%
  group_by(feature_id) %>%
  summarise(across(everything(), ~ ifelse(all(is.na(.x)), NA, max(.x, na.rm = TRUE)))) %>%
  ungroup()  # Remove grouping
blood_metab_qc_pre <- as.data.frame(blood_metab_qc_pre)
rownames(blood_metab_qc_pre) <- blood_metab_qc_pre$feature_id
blood_metab_qc_pre <- blood_metab_qc_pre[, -which(names(blood_metab_qc_pre) == "feature_id")] 

# Adipose metab qc data
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

filtered_metab_norm <- final_qc_norm_aggregated[rownames(final_qc_norm_aggregated) %in% precawg_metab_da$feature_id, ] 

# Align rownames between blood metab and adipose metab
common_rownames <- intersect(rownames(blood_metab_qc_pre), rownames(filtered_metab_norm))
blood_metab_qc_pre <- blood_metab_qc_pre[common_rownames, , drop = FALSE]
filtered_metab_norm <- filtered_metab_norm[common_rownames, , drop = FALSE]

# Align colnames
common_colnames <- intersect(colnames(blood_metab_qc_pre), colnames(filtered_metab_norm))
blood_metab_qc_pre <- blood_metab_qc_pre[, common_colnames, drop = FALSE]
filtered_metab_norm <- filtered_metab_norm[, common_colnames, drop = FALSE]
blood_metab_qc_pre <- t(blood_metab_qc_pre)
filtered_metab_norm <- t(filtered_metab_norm)
# Ensure they are now aligned
stopifnot(identical(rownames(blood_metab_qc_pre), rownames(filtered_metab_norm)))
stopifnot(identical(colnames(blood_metab_qc_pre), colnames(filtered_metab_norm)))

metab_sec_cor <- bicorAndPvalue(blood_metab_qc_pre, filtered_metab_norm, use = "pairwise.complete.obs")
df_m = melt(metab_sec_cor$bicor) %>% dplyr::rename(bicor = value)
df_m$pval = melt(metab_sec_cor$p)$value
df_m$obs = melt(metab_sec_cor$nObs)$value
df_m <- df_m %>% 
  mutate(sig = ifelse(pval < 0.001, '***',
                      ifelse(pval < 0.01, '**',
                             ifelse(pval < 0.05, '*', ''))))
df_m_filtered <- df_m %>% filter(as.character(Var1) == as.character(Var2))

# Figure S6C
adi_met_blood_met_plot <- ggplot(df_m_filtered, 
                                 aes(x = fct_reorder(Var1, desc(bicor)), y = bicor)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "purple") +  # Purple bars with black borders
  geom_text(aes(label = sig), vjust = -0.5, size = 3, color = "black") +  # Significance on top
  labs(y = "Bicor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_blank(),
        legend.position = "none")  # Remove legend


### Precawg adipose secretomes in GTEx
load('GTEx NA included env.RData')

GTEx_subfiltered$gene_tissue=NULL
GTEx_subfiltered = as.data.frame(t(GTEx_subfiltered))
tissue_names <- gsub("^[^_]+_", "", colnames(GTEx_subfiltered))

# Define the 90 genes (list from adipose)
gene_list <- unique(blood_protein_genes$gene_symbol)  # Add all 60 genes
gtex_gene_names <- unique(gsub("_.*", "", colnames(GTEx_subfiltered)))

# Check overlap between gene_list and GTEx gene names
matched_genes <- intersect(gene_list, gtex_gene_names)
print(matched_genes)

# Define the target tissues to compare against
target_tissues <- c("Adipose - Visceral (Omentum)", "Brain - Hypothalamus", 
                    "Brain - Hippocampus", "Small Intestine - Terminal Ileum",
                    "Stomach", "Thyroid", "Pancreas", "Spleen", "Muscle - Skeletal",
                    "Pituitary", "Artery - Coronary", "Liver", "Kidney - Cortex",
                    "Heart - Left Ventricle", "Colon - Transverse", "Colon - Sigmoid",
                    "Adrenal Gland", "Artery - Aorta")


target_tissue_cols <- colnames(GTEx_subfiltered)[grepl(paste(target_tissues, collapse = "|"), colnames(GTEx_subfiltered))]

# Function to filter genes
filter_low_variance_genes <- function(data) {
  zero_threshold <- 0.9  # Remove genes (columns) where >90% of values are zero
  mad_threshold <- 0.01  # Remove genes (columns) with MAD < 0.01
  
  # Calculate proportion of zeros per gene (column-wise)
  gene_zero_proportion <- colMeans(data == 0, na.rm = TRUE)
  
  # Compute MAD for each gene (column-wise)
  gene_mad <- apply(data, 2, mad, na.rm = TRUE)  # Change `1` to `2` for columns
  
  # Keep genes (columns) with <80% zeros and MAD > 0.01
  valid_genes <- (gene_zero_proportion < zero_threshold) & (gene_mad > mad_threshold)
  
  return(data[, valid_genes, drop = FALSE])  # Keep only valid columns (genes)
}


# Apply filtering to GTEx data
GTEx_filtered <- filter_low_variance_genes(GTEx_subfiltered)
cat("Number of genes retained after filtering:", ncol(GTEx_filtered), "\n") #Number of genes retained after filtering: 1110969 

compute_bicor_adipose <- function(gene, data, target_tissues) {
  # Step 1: Extract Adipose - Subcutaneous expression for the gene
  gene_colname <- paste0(gene, "_Adipose - Subcutaneous")
  
  if (!gene_colname %in% colnames(data)) {
    cat(paste0("Gene not found in dataset: ", gene, "\n"))
    return(NULL)
  }
  
  tissue1 <- data[, gene_colname, drop = FALSE]
  
  # Step 2: Escape special characters (like parentheses)
  target_tissues_fixed <- gsub("([()])", "\\\\\\1", target_tissues)  # Escape ( and )
  
  # Step 3: Extract expression for target tissues
  tissue2_cols <- colnames(data)[grepl(paste(target_tissues_fixed, collapse = "|"), colnames(data))]
  
  # Debugging: Print selected tissues
  cat("Processing", gene, "→ Selected Tissues:\n")
  print(tissue2_cols)
  
  # Debugging: Check if Adipose - Visceral (Omentum) is missing
  if (!any(grepl("Adipose - Visceral", tissue2_cols))) {
    cat("WARNING: Adipose - Visceral (Omentum) is STILL MISSING after regex fix!\n")
  }
  
  if (length(tissue2_cols) == 0) {
    cat(paste0("No matching tissues found for gene: ", gene, "\n"))
    return(NULL)
  }
  
  tissue2 <- data[, tissue2_cols, drop = FALSE]
  
  # Step 4: Remove zero-variance columns
  tissue_var <- apply(tissue2, 2, var, na.rm = TRUE)
  cat("Variance of tissues before filtering:\n")
  print(tissue_var)
  
  tissue2 <- tissue2[, tissue_var > 0, drop = FALSE]
  
  if (ncol(tissue2) == 0) {
    cat(paste0("All target tissues for gene ", gene, " have zero variance. Skipping...\n"))
    return(NULL)
  }
  
  # Step 5: Compute biweight midcorrelation
  full_cors <- bicorAndPvalue(tissue1, tissue2, use = "pairwise.complete.obs")
  
  if (is.null(full_cors$bicor) || is.null(full_cors$p)) {
    cat(paste0("Correlation failed for gene: ", gene, "\n"))
    return(NULL)
  }
  
  # Step 6: Reshape results
  cor_table <- reshape2::melt(full_cors$bicor)
  new_p <- reshape2::melt(full_cors$p)
  
  if (nrow(new_p) == 0) {
    cat(paste0("No valid p-values computed for gene: ", gene, "\n"))
    return(NULL)
  }
  
  cor_table$pvalue <- new_p$value
  colnames(cor_table) <- c('gene_tissue_1', 'gene_tissue_2', 'bicor', 'pvalue')
  
  # Step 7: Apply Benjamini-Hochberg (BH) correction per gene
  cor_table$qvalue <- p.adjust(cor_table$pvalue, method = "BH")
  
  return(cor_table)
}
correlation_results <- do.call(rbind, lapply(matched_genes, function(g) {
  compute_bicor_adipose(g, GTEx_filtered, target_tissues)
})) # This consumes time

# Extract exercise group (EE/RE) and omics type from blood_protein_genes
gene_metadata <- blood_protein_genes %>% 
  mutate(
    exercise_group = ifelse(grepl("EE", facet_label), "EE", "RE"),
    omics_type = ifelse(assay == "transcript-rna-seq", "Transcriptomics", "Proteomics")
  )
gene_metadata <- gene_metadata %>%
  group_by(gene_symbol) %>%
  summarize(
    exercise_group = ifelse(n() > 1, "Both", unique(exercise_group)),
    omics_type = unique(omics_type),
    .groups = "drop"
  )
gene_metadata <- as.data.frame(gene_metadata)

rownames(gene_metadata) <- gene_metadata$gene_symbol

correlation_results <- correlation_results %>%
  mutate(
    tissue = gsub(".*_", "", gene_tissue_2)  # Extract only the tissue name
  )

significant_results <- correlation_results %>% filter(qvalue < 0.05)

# Create heatmap matrix with count of significant genes per tissue
heatmap_matrix <- significant_results %>%
  group_by(gene_tissue_1, tissue) %>%
  summarise(num_genes = n(), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = num_genes, values_fill = 0) %>%
  rename(gene_symbol = gene_tissue_1) %>%
  column_to_rownames("gene_symbol")

heatmap_matrix <- as.matrix(heatmap_matrix)
heatmap_gene_names <- gsub("_.*", "", rownames(heatmap_matrix))
rownames(heatmap_matrix) <- heatmap_gene_names

gene_metadata <- gene_metadata[rownames(heatmap_matrix), ]

### Call in COMPARTMENTS data. 
compartments <-read_tsv(
  "https://download.jensenlab.org/human_compartment_integrated_full.tsv",
  col_names = c("feature_id", "gene_symbol", "GO", "location", "score"),
  show_col_types = FALSE
)

target_locations <- c("Extracellular region", "Extracellular space", "Extracellular exosome", "Extracellular vesicle") # These are the 4 locations that has 'extracellular' regions. 

# Filter the dataset for the selected genes and locations
compartments_filtered <- compartments %>%
  filter(gene_symbol %in% rownames(heatmap_matrix) & location %in% target_locations)

# Pick the highest extracellular score among 4 locations. 
compartments_filtered_max <- compartments_filtered %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() 
gene_metadata <- gene_metadata %>%
  left_join(compartments_filtered_max %>% select(gene_symbol, score), by = "gene_symbol")

# Define the white-to-green color scale for score
score_colors <- colorRamp2(
  c(min(gene_metadata$score, na.rm = TRUE), max(gene_metadata$score, na.rm = TRUE)),
  c("white", "green")
)


anno_colors <- list(
  Omics = c("Transcriptomics" = "#377EB8", "Proteomics" = "#228833"),
  Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77", "Both" = "black")
)
col_fun <- colorRamp2(
  seq(0, max(heatmap_matrix, na.rm = TRUE), length.out = 9),  # Auto-breaks across 9 levels
  brewer.pal(9, "YlGnBu") # Reversed YlGnBu color scale (Dark Blue → Light Yellow)
)

row_anno <- rowAnnotation(
  Omics = gene_metadata$omics_type,
  Exercise = gene_metadata$exercise_group,
  Score = anno_barplot(
    gene_metadata$score, 
    gp = gpar(fill = score_colors(gene_metadata$score)), 
    border = FALSE
  ),
  col = list(
    Omics = c("Transcriptomics" = "#377EB8", "Proteomics" = "#228833"),
    Exercise = c("EE" = "#d95f02", "RE" = "#1b9e77", "Both" = "black")
  )
)
# Figure 6D
Heatmap(
  heatmap_matrix,
  name = "Correlating genes",
  col = col_fun,
  rect_gp = gpar(col = "black"),  # Add grid lines
  left_annotation = row_anno,    # Correct annotation
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE
)
dev.off()

## Enriched pathways with CEACAM8 and BST2 from GTEx. 
correlation_results2 <- correlation_results %>%
  mutate(
    gene_symbol = gsub("_.*", "", gene_tissue_2)  # Remove everything after the first "_"
  )

#target_genes <- c("JUN", "CX3CL1", "COL6A1", "COL6A2", "ST3GAL1")
target_genes <- c("CEACAM8","BST2")

# Focusing on top tissues that have most correlating genes with subcutaneous adipose. 
target_tissues <- c("Muscle - Skeletal", "Adipose - Visceral (Omentum)", "Heart - Left Ventricle", "Thyroid")

# Filter, reshape, and set proper column names
target_genes_cor <- correlation_results2 %>%
  filter(
    gene_tissue_1 %in% paste0(target_genes, "_Adipose - Subcutaneous"),  # Filter genes
    tissue %in% target_tissues  # Filter target tissues
  ) %>%
  select(bicor, gene_symbol, tissue, gene_tissue_1) %>%  # Keep necessary columns
  mutate(gene_tissue = paste0(gene_tissue_1, "_", tissue)) %>%  # Create new column names
  select(bicor, gene_symbol, gene_tissue) %>%  # Keep relevant columns
  pivot_wider(names_from = gene_tissue, values_from = bicor) %>%  # Reshape to wide format
  column_to_rownames("gene_symbol")  

# prepare genesets
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
target_genes_cor_z <- as.data.frame(apply(target_genes_cor, 2, function(x) {
  z <- 0.5 * log((1 + x) / (1 - x))
  z[!is.finite(z)] <- NA  # remove +/- Inf
  return(z)
}))

gtex_camera_res <- cameraPR.matrix(
  statistic = as.matrix(target_genes_cor_z),
  index = genesets_list,
  use.ranks = FALSE,        # Use numeric values (not ranks)
  inter.gene.cor = 0.01,    # Default correlation adjustment
  alternative = "two.sided",
  min.size = 10             # Ignore very small sets
)

gtex_top_res <- gtex_camera_res %>%
  mutate(GeneSet_clean = gsub("^GOBP_", "", GeneSet)) %>%
  group_by(Contrast)

top_terms <- gtex_top_res %>%
  arrange(desc(ZScore)) %>% slice_head(n=2) %>%
  bind_rows(gtex_top_res %>% arrange(ZScore) %>% slice_head(n=1)) %>%
  pull(GeneSet_clean) %>% unique()

# 2. Column reordering info
column_info <- gtex_top_res %>%
  filter(GeneSet_clean %in% top_terms) %>%
  distinct(Contrast) %>%
  mutate(
    Gene   = sub("_Adipose - Subcutaneous.*", "", Contrast),
    Tissue = sub(".*_Adipose - Subcutaneous_", "", Contrast)
  ) %>%
  arrange(
    factor(Gene,   levels = c("CEACAM8", "BST2")),
    factor(Tissue, levels = c("Adipose - Visceral (Omentum)",
                              "Muscle - Skeletal",
                              "Heart - Left Ventricle",
                              "Thyroid"))
  )

heatmap_columns <- column_info$Contrast
column_order_numeric <- match(heatmap_columns, unique(gtex_top_res$Contrast))

# Annotation vectors aligned with heatmap_columns
column_genes   <- setNames(column_info$Gene,   column_info$Contrast)
column_tissues <- setNames(column_info$Tissue, column_info$Contrast)

# Colors
gene_colors <- c(
  "CEACAM8" = "#66C2A5",
  "BST2"    = "#FC8D62"
)

tissue_colors <- c(
  "Adipose - Visceral (Omentum)" = "yellow",
  "Muscle - Skeletal"            = "#377eb8",
  "Heart - Left Ventricle"       = "#e41a1c",
  "Thyroid"                      = "purple"
)

# Column Annotation
col_anno <- columnAnnotation(
  Gene   = column_genes[heatmap_columns],
  Tissue = column_tissues[heatmap_columns],
  col = list(Gene = gene_colors, Tissue = tissue_colors)
)

# Split by gene
column_split_vec <- column_genes[heatmap_columns]

# Figure 6F
gtex_top_res %>% 
  filter(GeneSet_clean %in% top_terms) %>%
  enrichmap(
    n_top = Inf,
    set_column = "GeneSet_clean",
    statistic_column = "ZScore",
    contrast_column  = "Contrast",
    padj_column = "FDR",
    padj_legend_title = "BH Adjusted\nP-Value",
    heatmap_args = list(
      heatmap_legend_param = list(title = "Z-Score"),
      top_annotation = col_anno,
      column_order = column_order_numeric,
      column_split = column_split_vec,
      show_column_names = TRUE
    ))



