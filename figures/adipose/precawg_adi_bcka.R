## Circulating BCKA to adipose protein synthesis
library(readxl)
library(MotrpacRatTraining6mo)
library(tidyr)
library(MotrpacHumanPreSuspensionAnalysis)
library(TMSig)
library(dplyr)
library(patchwork)
library(ggplot2)
library(purrr)

### Incorportate Walejko et al., Nat comm, 2021
adipose_files_path = file.path(here(), "figures", "adipose", "Files")
wal_phos_up <- read_excel(file.path(adipose_files_path, "wlejko_sup2.xlsx"), sheet = "Upregulated Phosphosites")
wal_phos_down <- read_excel(file.path(adipose_files_path, "wlejko_sup2.xlsx"), sheet = "Downregulated Phosphosites")

data("RAT_TO_HUMAN_GENE") #from MotrpacRatTraining6mo
uniprot_map_clean <- RAT_TO_HUMAN_GENE %>%
  separate_rows(RAT_SYMBOL, sep = ";") %>%
  distinct(RAT_SYMBOL, HUMAN_ORTHOLOG_SYMBOL)

# Step 2: Left join to your phospho data
wal_phos_up_annotated <- wal_phos_up %>%
  left_join(uniprot_map_clean, by = c("Gene ID (Uniprot)" = "RAT_SYMBOL"))
wal_phos_down_annotated <- wal_phos_down %>%
  left_join(uniprot_map_clean, by = c("Gene ID (Uniprot)" = "RAT_SYMBOL"))

wal_phos_up_prot <- wal_phos_up_annotated %>%
  pull(HUMAN_ORTHOLOG_SYMBOL) %>%
  unique() # 29 unique proteins

wal_phos_down_prot <- wal_phos_down_annotated %>%
  pull(HUMAN_ORTHOLOG_SYMBOL) %>%
  unique() # 79 unique proteins

precawg_phos_da_wal_up <- ph_vol %>% # ph_vol is from precawg_adi_da.R
  filter(gene_symbol %in% wal_phos_up_prot) %>%
  filter(adj_p_value<0.1)
setdiff(wal_phos_up_prot, unique(precawg_phos_da_wal_up$gene_symbol))
intersect(wal_phos_up_prot, unique(precawg_phos_da_wal_up$gene_symbol))
length(intersect(wal_phos_up_prot, unique(precawg_phos_da_wal_up$gene_symbol))) / length(wal_phos_up_prot)
# 15/28 DA phospho (51%)

precawg_phos_da_wal_down <- ph_vol %>%
  filter(gene_symbol %in% wal_phos_down_prot) %>%
  filter(adj_p_value<0.1)
setdiff(wal_phos_down_prot, unique(precawg_phos_da_wal_down$gene_symbol))
intersect(wal_phos_down_prot, unique(precawg_phos_da_wal_down$gene_symbol))
length(intersect(wal_phos_down_prot, unique(precawg_phos_da_wal_down$gene_symbol))) / length(wal_phos_down_prot)

# 34/78 DA phospho (43%)

# Merge up/downregulated features
wal_phos_da <- union(wal_phos_up_prot, wal_phos_down_prot) %>%
  discard(~ is.na(.x) || .x == "NA") # 102
precawg_phos_da_wal <- ph_vol %>%
  filter(gene_symbol %in% wal_phos_da) %>%
  filter(adj_p_value<0.1)
length(unique(precawg_phos_da_wal$gene_symbol)) #47. 47/102. 46$ being DA

#### ORA
wal_phos_signatures <- list(
  wal_phos_up = intersect(wal_phos_up_prot, unique(precawg_phos_da_wal_up$gene_symbol)),
  wal_phos_down = intersect(wal_phos_down_prot, unique(precawg_phos_da_wal_down$gene_symbol))
)
bckg_ph <- as.character(unique(precawg_phos_da$gene_symbol))

ora_precawg_wal <- lapply(wal_phos_signatures, function(input_i) {
  run_ORA(input = as.character(input_i),
          background = bckg_ph,
          overlap_cutoff = 0)
}) %>%
  bind_rows(.id = "Direction") %>%
  mutate(log10p = -log10(p_value),
         Direction = factor(Direction, levels = unique(Direction)))

top_terms <- ora_precawg_wal %>%
  group_by(Direction) %>%  # Group by Module to apply slice_min() within each group
  pull(set_short) %>%  # Extract the set_short column
  unique()

length(top_terms)
log10p_colors <- circlize::colorRamp2(c(0, max(ora_precawg_wal$log10p)), c("white", "#543483"))
ora_precawg_wal <- ora_precawg_wal %>%
  mutate(
    Direction = recode(Direction,
                       "wal_phos_down" = "Hypophosphorylation",
                       "wal_phos_up" = "Hyperphosphorylation"),
    Direction = factor(Direction, levels = c("Hypophosphorylation", "Hyperphosphorylation")),
    log10p = -log10(p_value)
  )
# Figure 5H
ora_precawg_wal %>%
  enrichmap(n_top = Inf,
            set_column = "set_short",
            statistic_column = "log10p",
            contrast_column = "Direction",
            padj_column = "p_value",
            padj_legend_title = "P-Value",
            padj_cutoff = 0.001,
            plot_sig_only = TRUE,
            colors = c("white", "#543483"),
            heatmap_args = list(heatmap_legend_param = list(
              title = latex2exp::TeX("$\\bf{$-log$_{10}($P-Value$)}$")
            )))



# Figure 5G: Transcriptional enrichment of protein translation
selected_contrasts <- c(
  "Endur.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)",
  "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
  "Endur.post_24_hr - Control.post_24_hr (delta-delta)",
  "Resist.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)",
  "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)",
  "Resist.post_24_hr - Control.post_24_hr (delta-delta)"
)
protsyn_t <- CAMERA_RESULTS %>%
  filter(grepl("ribosome|translation", set_short, ignore.case = TRUE)) %>%
  filter(tissue == 'adipose', assay == "transcript-rna-seq", database == "REACTOME") %>%
  filter(!grepl("SARS_COV", set_short, ignore.case = TRUE)) %>%
  filter(contrast_short %in% selected_contrasts) %>%
  mutate(
    contrast_short = factor(contrast_short, levels = selected_contrasts)
  )


# Define exercise & control colors
Exercise_colors <- c("EE" = "#d95f02", "RE" = "#1b9e77")
time_colors <- c("45minPost" = "#AE76A3", "4hrPost" = "#882E72", "24hrPost" = "#61194F")  # Timepoint colors

# Extract Exercise Type & Time for annotation
protsyn_t <- protsyn_t %>%
  mutate(
    Timepoint = case_when(
      grepl("24_hr", contrast) ~ "24hrPost",
      grepl("4_hr", contrast) ~ "4hrPost",
      grepl("45_min", contrast) ~ "45minPost",
      TRUE ~ NA_character_  # fallback
    ),
    Exercise = case_when(
      grepl("Endur", contrast_short) ~ "EE",
      grepl("Resist", contrast_short) ~ "RE"
    )
  )
col_anno <- columnAnnotation(
  Exercise = rep(c("EE", "RE"), each = 3),
  Time = rep(c("45minPost", "4hrPost", "24hrPost"), times = 2),
  col = list(
    Exercise = Exercise_colors,
    Time = time_colors
  ),
  annotation_legend_param = list(
    Exercise = list(title = "Modality", at = c("EE", "RE")),
    Time = list(title = "Timepoint", at = c("45minPost", "4hrPost", "24hrPost"))
  )
)
# Figure 5G
protsyn_t_plot <- protsyn_t %>%
  enrichmap(
    n_top = Inf,
    plot_sig_only = TRUE,
    set_column = "set_short",
    statistic_column = "z.std",
    contrast_column = "contrast",
    padj_column = "adj_p_value",
    padj_legend_title = "BH Adjusted\nP-Value",
    heatmap_args = list(
      cluster_rows = TRUE,
      top_annotation = col_anno,
      column_split = rep(c("EE", "RE"), each = 3),
      show_column_names = FALSE
    )
  )

