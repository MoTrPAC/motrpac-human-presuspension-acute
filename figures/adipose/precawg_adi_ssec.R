library(MotrpacHumanPreSuspensionAnalysis)
library(ggtext)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(ggtext)
library(patchwork)
library(RColorBrewer)
library(stats)
library(tibble)
library(circlize)
library(ComplexHeatmap)

files_path = file.path(here(), "figures", "adipose", "Files")
wat_combined_all <- readRDS(file.path(files_path, "wat_ssec_sexadj.RDS")) # data from PMCID: PMC12340562
#this file size is too large to store in the github repository, please refer to the manuscript to get access to this data.

wat_combined_all <- wat_combined_all %>%
  mutate(Facet_Label = paste0("WATSC → ", Target_Tissue))
training_colors <- brewer.pal(8, "YlGnBu")  # Generates a 4-color gradient from the palette
condition_colors <- c(
  "CON"  = "#000000B3",     # Control
  "TR1W" = training_colors[2],
  "TR2W" = training_colors[4],
  "TR4W" = training_colors[6],
  "TR8W" = training_colors[8]
)

# Figure 7B
d_facet <- ggplot(
  wat_combined_all,
  aes(x = Ssec, fill = Condition, color = Condition)
) +
  geom_density(alpha = 0.4, size = 1) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  facet_wrap(~ Facet_Label, ncol = 4, scales = "free") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    title = NULL,
    x = "Ssec",
    y = "Density"
  )

# Perform wilcoxon test to compare CON and TR9W focusing on 60 secretome candidate genes.
# blood_protein_genes is obtained from precawg_adi_secretome.R.
#str(blood_protein_genes)
#tibble [63 × 3] (S3: tbl_df/tbl/data.frame)
#$ assay      : Factor w/ 2 levels "transcript-rna-seq",..: 1 1 1 1 1 1 1 1 1 1 ...
#$ facet_label: chr [1:63] "EE 45min" "EE 45min" "EE 45min" "EE 45min" ...
#$ gene_symbol: chr [1:63] "CCN1" "CCN2" "ADAMTS1" "AREG" ...

wat_filtered <- wat_combined_all %>%
  filter(Human_Gene %in% unique(blood_protein_genes$gene_symbol)) %>%
  rename(Tissue = Facet_Label)

wilcoxon_test_ssec <- function(df) {

  df_filtered <- df %>%
    filter(Condition %in% c("CON", "TR8W"))

  # If either group is missing → return NA
  if (!all(c("CON", "TR8W") %in% unique(df_filtered$Condition))) {
    return(data.frame(Tissue = unique(df$Tissue), P_Value = NA))
  }

  # Require paired genes (each gene appearing in both CON & TR8W)
  wide_df <- df_filtered %>%
    select(Human_Gene, Condition, Ssec) %>%
    pivot_wider(names_from = Condition, values_from = Ssec)

  # Remove incomplete pairs
  wide_df <- wide_df %>% drop_na()

  if (nrow(wide_df) < 3) {
    return(data.frame(Tissue = unique(df$Tissue), P_Value = NA))
  }

  test_result <- wilcox.test(wide_df$CON,
                             wide_df$TR8W,
                             paired = TRUE)

  data.frame(Tissue = unique(df$Tissue),
             P_Value = test_result$p.value)
}

wilcoxon_results <- wat_filtered %>%
  group_by(Tissue) %>%
  group_split() %>%
  map_df(wilcoxon_test_ssec) %>%
  mutate(
    adjp = p.adjust(P_Value, method = "BH"),
    LogP = -log10(adjp),
    Significance = case_when(
      adjp < 0.001 ~ "***",
      adjp < 0.01  ~ "**",
      adjp < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Figure 7C
ggplot(wilcoxon_results,
       aes(x = reorder(Tissue, LogP), y = LogP)) +
  geom_col(fill = "gray80", color = "black") +
  geom_text(aes(label = Significance),
            hjust = -0.2, size = 6, fontface = "bold") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Control 8-week vs. Trained 8-week",
       x = NULL,
       y = "-log10(adjusted p-value)") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Function for paired t-test
paired_t_test_ssec <- function(df) {

  df_filtered <- df %>%
    filter(Condition %in% c("CON", "TR8W"))

  # Ensure both groups exist
  if (!all(c("CON", "TR8W") %in% unique(df_filtered$Condition))) {
    return(data.frame(Tissue = unique(df$Tissue), P_Value = NA))
  }

  # Paired genes (must appear once in CON & once in TR8W)
  wide_df <- df_filtered %>%
    select(Human_Gene, Condition, Ssec) %>%
    pivot_wider(names_from = Condition, values_from = Ssec) %>%
    drop_na()

  if (nrow(wide_df) < 3) {
    return(data.frame(Tissue = unique(df$Tissue), P_Value = NA))
  }

  # Paired t-test
  test_result <- t.test(wide_df$CON, wide_df$TR8W, paired = TRUE)

  data.frame(Tissue = unique(df$Tissue),
             P_Value = test_result$p.value)
}
t_test_results <- wat_combined_all %>%
  filter(Human_Gene %in% unique(blood_protein_genes$gene_symbol)) %>%
  rename(Tissue = Facet_Label) %>%
  group_by(Tissue) %>%
  group_split() %>%
  map_df(paired_t_test_ssec) %>%
  mutate(adjp = p.adjust(P_Value, method = "BH"))

# Prepare for plotting
ssec_combined <- wat_combined_all %>%
  filter(Human_Gene %in% unique(blood_protein_genes$gene_symbol),
         Condition %in% c("CON", "TR8W")) %>%
  rename(Tissue = Facet_Label) %>%
  left_join(t_test_results, by = "Tissue") %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)))

adjp_labels <- t_test_results %>%
  mutate(label = paste0("adj p = ", signif(adjp, 2)))

# Figure S7A
ggplot(ssec_combined, aes(x = Condition, y = Ssec, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~Tissue, ncol = 5, scales = "fixed") +
  scale_fill_manual(
    values = c("CON" = "pink", "TR8W" = "#002612"),
    name = "Condition"
  ) +
  coord_cartesian(ylim = c(0.2, 1.5)) +
  labs(y = "Ssec", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    strip.text = element_text(size = 10, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) +
  geom_text(
    data = adjp_labels,
    aes(x = 1.5, y = 1.48, label = label),
    inherit.aes = FALSE,
    size = 3.5
  )


# Custom function to plot top 15 different features with largest Ssec difference between CON and TR8W.
plot_top_rank_change <- function(origin_tissue, target_tissue) {
  # Step 1: Subset and rank
  con_data <- wat_combined_all %>%
    filter(Condition == "CON", Origin_Tissue == origin_tissue, Target_Tissue == target_tissue) %>%
    mutate(Rank_control = rank(-Ssec, na.last = "keep"))

  tr8w_data <- wat_combined_all %>%
    filter(Condition == "TR8W", Origin_Tissue == origin_tissue, Target_Tissue == target_tissue) %>%
    mutate(Rank_tr8 = rank(-Ssec, na.last = "keep")) %>%
    select(RAT_SYMBOL, Rank_tr8, Ssec, score, Human_Gene)

  # Step 2: Merge
  ranked_genes <- merge(
    con_data,
    tr8w_data[, c("RAT_SYMBOL", "Rank_tr8", "Ssec", "score", "Human_Gene")],
    by = "RAT_SYMBOL", suffixes = c("_control", "_tr8")
  ) %>%
    mutate(
      Rank_change = abs(Rank_tr8 - Rank_control)
    )

  # Step 3: Top 15 genes
  top_10_genes <- ranked_genes %>%
    filter(Human_Gene_tr8 %in% unique(blood_protein_genes$gene_symbol)) %>%
    arrange(desc(Rank_change)) %>%
    slice(1:15) %>%
    mutate(
      display_label = factor(RAT_SYMBOL, levels = .$RAT_SYMBOL)
    )

  # Step 4: Main Rank Plot (with arrow, flipped Y)
  p1 <- ggplot(top_10_genes, aes(x = display_label)) +
    geom_segment(
      aes(
        y = Rank_control + 0.2 * (Rank_tr8 - Rank_control),
        yend = Rank_tr8 - 0.2 * (Rank_tr8 - Rank_control),
        x = display_label, xend = display_label
      ),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      color = "gray40", size = 1.2
    )+
    geom_point(aes(y = Rank_control), shape = 24, fill = "pink", size = 3, color = "black") +
    geom_point(aes(y = Rank_tr8), shape = 21, fill = "#002612", size = 4, color = "black") +
    scale_y_reverse() +
    labs(title = paste(origin_tissue, "→", target_tissue),
         y = "Ssec Rank",
         x = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    )

  # Step 5: Horizontal bar plot for score
  p2 <- ggplot(top_10_genes, aes(x = display_label, y = score_tr8, fill = score_tr8)) +
    geom_col(width = 0.6, color = "black") +
    scale_fill_gradient(low = "lightgreen", high = "darkgreen", na.value = "white") +
    labs(y = "Extracellular\nScore", x = NULL) +  # Line break
    theme_minimal() +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 7),
      legend.position = "none",
      axis.title.y = element_text(size = 9)  # Smaller font for Y-axis title
    )

  # Combine plots with patchwork
  p1 / p2 + plot_layout(heights = c(3, 1))
}

# Figure 7E
plot_top_rank_change("WAT-SC", "SKM-VL")
plot_top_rank_change("WAT-SC", "LIVER")
plot_top_rank_change("WAT-SC", "KIDNEY")


## Show Ssec changes for 60 secretome candidates in all tissue.
wat_combined_all2 <- wat_combined_all %>%
  filter(Human_Gene %in% blood_protein_genes$gene_symbol)

# Step 3: Pivot to wide format to get TR8W - CON per gene and tissue
ssec_diff_60 <- wat_combined_all2 %>%
  select(Human_Gene, Target_Tissue, Condition, Ssec) %>%
  pivot_wider(names_from = Condition, values_from = Ssec) %>%
  mutate(diff = TR8W - CON) %>%
  select(Human_Gene, Target_Tissue, diff)

# Step 4: Create matrix: rows = gene, cols = target tissue
ssec_diff_60_matrix <- ssec_diff_60 %>%
  pivot_wider(names_from = Target_Tissue, values_from = diff) %>%
  column_to_rownames("Human_Gene") %>%
  as.matrix()

# Step 5: Z-score across columns (per gene)
ssec_z_matrix <- t(scale(t(ssec_diff_60_matrix)))

# Step 6: Color function
col_fun <- colorRamp2(
  c(min(ssec_z_matrix, na.rm = TRUE), 0, max(ssec_z_matrix, na.rm = TRUE)),
  c("pink", "white", "#002612")
)

# Figure S7B
Heatmap(
  ssec_z_matrix,
  name = "Z-score\nSsec TR8W-CON",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(min(ssec_z_matrix, na.rm = TRUE), 0, max(ssec_z_matrix, na.rm = TRUE)),
    labels = round(c(min(ssec_z_matrix, na.rm = TRUE), 0, max(ssec_z_matrix, na.rm = TRUE)), 1)
  )
)
