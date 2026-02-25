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
  "circlize",
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "dplyr",
  "ggplot2",
  "ggh4x",
  "ggrepel",
  "ggpubr",
  "hrbrthemes",
  "magrittr",
  "MotrpacHumanPreSuspension",
  "MotrpacHumanPreSuspensionAnalysis",
  "purrr",
  "RColorBrewer",
  "RefMet",
  "stringr",
  "tidyverse"
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

###############################################################################
sel_tissue <- "muscle"
sel_ome <- "metab"
sel_epigen <- FALSE
contrasts <- c("EE post 15 min",
               "EE post 3.5 hr",
               "EE post 24 hr",
               "RE post 15 min",
               "RE post 3.5 hr",
               "RE post 24 hr")

metab_da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = sel_ome,
                                                                      selected_tissues = sel_tissue,
                                                                      epigen = sel_epigen,
                                                                      combine_with_featgene = T,
                                                                      single_matrix = T) %>%
  filter(contrast_type=="exercise_with_controls") %>%
  mutate(contrast_short=gsub(" - .*","",contrast_short))

metab_sig <- metab_da %>%
  filter(adj_p_value < 0.05)

metab_sig2 <- metab_sig %>%
  mutate(contrast_short=case_when(contrast_short=="Endur.post_15_30_45_min"~"EE post 15 min",
                                  contrast_short=="Endur.post_3.5_4_hr"~"EE post 3.5 hr",
                                  contrast_short=="Endur.post_24_hr"~"EE post 24 hr",
                                  contrast_short=="Resist.post_15_30_45_min"~"RE post 15 min",
                                  contrast_short=="Resist.post_3.5_4_hr"~"RE post 3.5 hr",
                                  contrast_short=="Resist.post_24_hr"~"RE post 24 hr")) %>%
  mutate(ome=case_when(assay=="prot-ol"~"Proteomics OL",
                       assay=="prot-pr"~"Proteomics PR",
                       assay=="prot-ph"~"Phosphoproteomics",
                       assay=="transcript-rna-seq"~"Transcriptomics",
                       assay=="epigen-atac-seq"~"ATAC-seq",
                       assay=="epigen-methylcap-seq"~"MethylCap-seq",
                       .default = "Metabolomics")) %>%
  mutate(ex_modality=gsub("_.*","",contrast_short)) %>%
  mutate(refmet_name = case_when(is.na(refmet_name) ~ feature_id,
                                 .default = refmet_name))

RefMet_mapped <- refmet_map_df(unique(metab_sig2$refmet_name)) %>%
  mutate(Sub.class = case_when(Sub.class=="Carnitines" ~ "Acyl carnitines",
                               .default = Sub.class)) %>%
  rename(refmet_name=Standardized.name)

metab_sig_class <- full_join(metab_sig2,
                              RefMet_mapped %>%
                                select(refmet_name,Sub.class),
                              by = "refmet_name") %>%
  distinct(feature_id,contrast, .keep_all = TRUE) %>%
  mutate(Timepoint = recode(Timepoint,
                            post_15_30_45_min = "post 15 min",
                            post_3.5_4_hr = "post 3.5 hr",
                            post_24_hr = "post 24 hr"))

################################################################################
# function to extract chain length and unsaturation (only works for acyl carnitines!)
extract_acylcarnitine_info <- function(feature_id) {
  # Function to parse each individual feature id
  parse_feature <- function(fid) {
    chain_length <- NA
    double_bonds <- NA
    
    # Free carnitine
    if (fid == "Carnitine") {
      chain_length <- 0
      double_bonds <- 0
    }
    # Special case: "3-Dehydroxycarnitine"
    else if (fid == "3-Dehydroxycarnitine") {
      chain_length <- 4
      double_bonds <- 0
    }
    # Branched chain: "CAR DC3:0;2Me"
    else if (grepl("CAR DC3:0", fid)) {
      chain_length <- 3
      double_bonds <- 0
    }
    # Common CAR pattern: "CAR X:Y", "CAR X:Y;OH", "CAR X:Y;3Me"
    else if (grepl("CAR ", fid)) {
      matches <- regmatches(fid, regexec("CAR ([0-9]+):([0-9]+)", fid))
      if (length(matches[[1]]) == 3) {
        chain_length <- as.numeric(matches[[1]][2])
        double_bonds <- as.numeric(matches[[1]][3])
      }
    }
    
    return(c(chain_length, double_bonds))
  }
  
  # Apply parsing to every feature_id
  parsed <- t(sapply(feature_id, parse_feature))
  
  # result data frame
  df <- data.frame(
    feature_id = feature_id,
    chain_length = parsed[,1],
    double_bonds = parsed[,2],
    stringsAsFactors = FALSE
  )
  return(df)
}

################################################################################
# Acylcarnitine plot
metab_class <- "Acyl carnitines"
features_to_plot <- metab_sig_class %>%
  filter(Sub.class==metab_class) %>%
  select(feature_id) %>%
  distinct() %>%
  pull()

acylcar_chain <- extract_acylcarnitine_info(features_to_plot)

sel_metabs_da <- metab_da %>%
  filter(feature_id %in% features_to_plot) %>%
  separate(contrast_short, into = c("Modality","Timepoint"),
           remove = FALSE, extra = "merge") %>%
  mutate(
    Timepoint = case_when(
      Timepoint == "post_15_30_45_min" ~ "post 15/30/45 min",
      Timepoint == "post_3.5_4_hr" ~ "post 3.5/4 hr", 
      TRUE ~ Timepoint
    ),
    Timepoint = gsub("_", " ", Timepoint)
  ) %>%
  filter(assay == "metab") %>%
  group_by(feature_id) %>%
  filter(any(adj_p_value < 0.05)) %>%
  ungroup()

sel_metab_small <- sel_metabs_da %>% 
  select(feature_id, 
         randomGroupCode, 
         Timepoint, 
         contrast_short, 
         z.std, 
         adj_p_value)

merged_data <- sel_metab_small %>%
  left_join(acylcar_chain, 
            by = "feature_id") %>%
  mutate(randomGroupCode = case_when(randomGroupCode == "ADUEndur" ~ "EE",
                                     randomGroupCode == "ADUResist" ~ "RE"),
         Timepoint = case_when(Timepoint == "post 15/30/45 min" ~ "post 15 min",
                               Timepoint == "post 3.5/4 hr" ~ "post 3.5 hr",
                               .default = Timepoint)) %>%
  mutate(sig_group = case_when(
    adj_p_value < 0.05 & randomGroupCode == "EE" ~ "EE_sig",
    adj_p_value < 0.05 & randomGroupCode == "RE" ~ "RE_sig",
    TRUE ~ "NS"
  )) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("post 15 min",
                                                  "post 3.5 hr",
                                                  "post 24 hr")))


# Plot
ggplot(merged_data, aes(x = chain_length, y = z.std)) +
  geom_point(
    aes(fill = sig_group),
    shape = 21, size = 3, color = "black", stroke = 0.5
  ) +
  geom_text_repel(
    data = merged_data %>% filter(abs(z.std) > 4),
    aes(label = feature_id),
    size = 4,
    fontface = "italic",
    color = "black",
    max.overlaps = 10,
    min.segment.length = 0,
    box.padding = 0.3,
    show.legend = FALSE
  ) +
  geom_smooth(method = "loess", se = TRUE, color = "grey20", linetype = "solid", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_grid(randomGroupCode ~ Timepoint, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "RE_sig" = "#1b9e77",    
      "EE_sig" = "#d95f02",    
      "NS"     = "white"       
    ),
    breaks = c("EE_sig", "RE_sig", "NS"),
    labels = c(
      "EE_sig"  = "EE sig",
      "RE_sig"  = "RE sig",
      "NS"      = "Not significant"
    ),
    name = "Significance"
  ) +
  labs(
    x = "Chain Length",
    y = "z.std"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10)
  )
# ggsave(paste0("muscle_acylcarnitines_plot_by_chain_length.jpeg"),
#        height = 8, width = 12, units = "in", dpi = 1200)
ggsave(paste0("precovid_muscle_figureS2B_acylCar_plot.pdf"),
       height = 8, width = 12, units = "in", dpi = 1200)