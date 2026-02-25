#This code uses a configuration file to ensure local file structure and local machine settings can be used properly.

#' Steps for using the MoTrPAC config file:
#' 1. Open the configuration template file in the root folder of this repo called `example_config.json` available in `precovid-analyses/`
#' 2. Save a copy of this template to some other location in your computer. Many choose to place it in their home directory: Example: `~/motrpac_config.json`
#' 3. Customize the fields for your machine and save changes. Remove any "comments" as JSON does not permit comments
#' 4. This script expects the following elements in the config file:
##' a. "gsutil_user" and "gsutil_path": this is the way your machine calls the gsutil, usually "/usr/bin/gsutil"
##' b. "local_path": this is the absolute file path for the `precovid-analyses` repo. Please do not include the final "/"
#' 5. Add the path to this config file below 
# config <- jsonlite::read_json(path = "//mnt/c/Users/dankatz/Documents/motrpac_config.json")
rm(list=ls())
setwd("C:/Users/gayat/Documents/MoTrPAC/precovid_analyses")
config <- rjson::fromJSON(file = 'C:/Users/gayat/Documents/MoTrPAC/motrpac_precovid_config_GI.json')

session_packages = c(
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "cowplot",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "ggrepel",
  "hrbrthemes",
  "magrittr",
  "MotrpacHumanPreSuspension",
  "MotrpacHumanPreSuspensionAnalysis",
  "patchwork",
  "purrr",
  "stringr",
  "tidyverse",
  "TMSig"
)
for (lib_name in session_packages){
  tryCatch({library(lib_name,character.only = TRUE)}, error = function(e) {
    print(paste("Cannot load",lib_name,", please install"))
  })
}

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("union","base")
conflict_prefer("intersect","base")
conflict_prefer("setdiff", "base")
conflict_prefer("HUMAN_FEATURE_TO_GENE","MotrpacHumanPreSuspensionAnalysis")
conflict_prefer("CONTRAST_CONVERTER","MotrpacHumanPreSuspensionAnalysis")

repo_local_dir <- config$precovid_repo_path
gsutil_cmd <- config$gsutil_user

################################################################################
sel_tissue <- "muscle"
sel_ome <- "all"
contrasts <- c("EE post 15 min",
               "EE post 3.5 hr",
               "EE post 24 hr",
               "RE post 15 min",
               "RE post 3.5 hr",
               "RE post 24 hr")

all_da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = sel_ome,
                                                                        selected_tissues = sel_tissue,
                                                                        epigen = TRUE,
                                                                        combine_with_featgene = T,
                                                                        single_matrix = T) %>%
  filter(contrast_type=="exercise_with_controls") %>%
  mutate(contrast_short=gsub(" - .*","",contrast_short))

all_sig <- all_da %>%
  filter(adj_p_value < 0.05)

all_sig2 <- all_sig %>%
  mutate(contrast_short=case_when(contrast_short=="Endur.post_15_30_45_min"~"EE post 15 min",
                                  contrast_short=="Endur.post_3.5_4_hr"~"EE post 3.5 hr",
                                  contrast_short=="Endur.post_24_hr"~"EE post 24 hr",
                                  contrast_short=="Resist.post_15_30_45_min"~"RE post 15 min",
                                  contrast_short=="Resist.post_3.5_4_hr"~"RE post 3.5 hr",
                                  contrast_short=="Resist.post_24_hr"~"RE post 24 hr")) %>%
  mutate(ome=case_when(assay=="prot-ol"~"Proteomics",
                       assay=="prot-pr"~"Proteomics",
                       assay=="prot-ph"~"Phosphoproteomics",
                       assay=="transcript-rna-seq"~"Transcriptomics",
                       assay=="epigen-atac-seq"~"Chromatin Accessibility (ATAC)",
                       assay=="epigen-methylcap-seq"~"Methylation",
                       .default = "Metabolomics")) %>%
  mutate(ex_modality=gsub("_.*","",contrast_short))

upset_data_all_features <- all_sig2 %>%
  #mutate(feature_contrast = paste0(feature_id,"_",contrast_short)) %>%
  #distinct(feature_contrast, .keep_all = T) %>%
  mutate(feature_id2=feature_id) %>%
  select(c(feature_id,feature_id2,contrast_short)) %>%
  pivot_wider(names_from = contrast_short,
              values_from = feature_id2) %>%
  mutate(across(-feature_id, ~ ifelse(is.na(.), 0, 1)))

upset_data_all_features2 <- inner_join(upset_data_all_features,
                                       all_sig2 %>%
                                         select(c(feature_id,ome)),
                                       by="feature_id") %>%
  distinct(feature_id, .keep_all = T) 


#################################################################################
# Scatter plots of shared features at each timepoints
filter_timepoint <- function(data, EE_values, RE_values) {
  data %>%
    filter(
      `EE post 15 min` == EE_values[1],
      `EE post 3.5 hr` == EE_values[2],
      `EE post 24 hr` == EE_values[3],
      `RE post 15 min` == RE_values[1],
      `RE post 3.5 hr` == RE_values[2],
      `RE post 24 hr` == RE_values[3]
    )
}

# Define patterns for each timepoint
pattern_15min <- list(EE = c(1, 0, 0), RE = c(1, 0, 0))
pattern_3.5hr <- list(EE = c(0, 1, 0), RE = c(0, 1, 0))
pattern_24hr <- list(EE = c(0, 0, 1), RE = c(0, 0, 1))

# Apply the function for each group
ee_re_15min_only <- filter_timepoint(
  upset_data_all_features2,
  pattern_15min$EE,
  pattern_15min$RE
)

ee_re_3.5hr_only <- filter_timepoint(
  upset_data_all_features2,
  pattern_3.5hr$EE,
  pattern_3.5hr$RE
)

ee_re_24hr_only <- filter_timepoint(
  upset_data_all_features2,
  pattern_24hr$EE,
  pattern_24hr$RE
)

###########
# get shared features at each timepoint
get_shared_features <- function(data, HUMAN_FEATURE_TO_GENE, all_sig2, timepoint, timepoint_label) {
  data %>%
    select(feature_id, ome) %>%
    left_join(
      HUMAN_FEATURE_TO_GENE %>% select(feature_id, gene_symbol), 
      by = "feature_id"
    ) %>%
    semi_join(
      all_sig2 %>% 
        filter(Timepoint == timepoint) %>% 
        select(feature_id),
      by = "feature_id"
    ) %>%
    left_join(
      all_sig2 %>%
        filter(Timepoint == timepoint) %>%
        select(feature_id, randomGroupCode, logFC),
      by = "feature_id"
    ) %>%
    select(feature_id, gene_symbol, ome, randomGroupCode, logFC) %>%
    mutate(
      gene_symbol = if_else(is.na(gene_symbol), feature_id, gene_symbol),
      randomGroupCode = recode(randomGroupCode,
                               "ADUEndur" = "EE",
                               "ADUResist" = "RE"
      )
    )
}


shared_15min <- get_shared_features(
  ee_re_15min_only,
  HUMAN_FEATURE_TO_GENE,
  all_sig2,
  "post_15_30_45_min",
  "15min"
)

shared_3.5hr <- get_shared_features(
  ee_re_3.5hr_only,
  HUMAN_FEATURE_TO_GENE,
  all_sig2,
  "post_3.5_4_hr",
  "3.5hr"
)

shared_24hr <- get_shared_features(
  ee_re_24hr_only,
  HUMAN_FEATURE_TO_GENE,
  all_sig2,
  "post_24_hr",
  "24hr"
)


###
# ome colors
ome_colors <- c(
  'Transcriptomics' = "#4477AA",
  'Metabolomics' = "#6D4B08",
  'Proteomics' = "#228833",
  'Phosphoproteomics' = "#F3A02B",
  'Chromatin Accessibility (ATAC)' = "#882255",
  'Methylation' = "#D687B5"
)


# function to create scatter plot
plot_logFC_scatter <- function(input_df, plot_title, top_n = 5) {

  scatter_data <- input_df %>%
    select(feature_id, gene_symbol, ome, randomGroupCode, logFC) %>%
    pivot_wider(
      names_from = randomGroupCode,
      values_from = logFC,
      names_prefix = "logFC_"
    )
  
  #Quadrants, distance
  scatter_data <- scatter_data %>%
    mutate(
      quadrant = case_when(
        logFC_EE > 0 & logFC_RE > 0 ~ "Q1",
        logFC_EE < 0 & logFC_RE > 0 ~ "Q2",
        logFC_EE < 0 & logFC_RE < 0 ~ "Q3",
        logFC_EE > 0 & logFC_RE < 0 ~ "Q4",
        TRUE ~ NA_character_ # Handles zero case
      ),
      distance = sqrt(logFC_EE^2 + logFC_RE^2)
    )
  
  # top N per quadrant
  top_features_by_quadrant <- scatter_data %>%
    filter(!is.na(quadrant)) %>%
    group_by(quadrant) %>%
    arrange(desc(distance)) %>%
    slice_head(n = top_n) %>%
    pull(gene_symbol)
  
  if (plot_title!="15 min") {
    top_features_by_quadrant <- c(top_features_by_quadrant,
                                  (input_df %>% 
                                     filter(ome!="Transcriptomics") %>% 
                                     pull(gene_symbol)))
  }
  
  # labeling
  scatter_data <- scatter_data %>%
    mutate(highlight = ifelse(gene_symbol %in% top_features_by_quadrant, "yes", "no"))
  
  # compute axis limits
  all_vals <- c(scatter_data$logFC_EE, scatter_data$logFC_RE)
  axis_max <- max(abs(all_vals), na.rm = TRUE)
  axis_limit <- c(-axis_max, axis_max)
  
  p <- ggplot(scatter_data, aes(
    x = logFC_EE,
    y = logFC_RE,
    color = ome
  )) +
    geom_point(size = 3, alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 1) +
    geom_text_repel(
      data = subset(scatter_data, highlight == "yes"),
      aes(label = gene_symbol),
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = "grey60",
      size = 3,
      show.legend = FALSE
    ) +
    labs(
      title = plot_title,
      x = "logFC (EE)",
      y = "logFC (RE)",
      color = "Ome"
    ) +
    scale_color_manual(values = ome_colors, drop = FALSE) +
    theme_cowplot(font_size = 14) +
    #coord_fixed(xlim = axis_limit, ylim = axis_limit)
    coord_fixed(xlim = c(-5,5), ylim = c(-5,5))
  
  return(p)
}

p1 <- plot_logFC_scatter(input_df = shared_15min, plot_title = "15 min")
p2 <- plot_logFC_scatter(input_df = shared_3.5hr, plot_title = "3.5 hr")
p3 <- plot_logFC_scatter(input_df = shared_24hr, plot_title = "24 hr")

combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 1)

ggsave("precovid_muscle_figure2C_scatterPlot.pdf", 
       combined_plot, height = 16, width = 8, dpi = 1200)
