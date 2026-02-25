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
  "dplyr",
  "ggplot2",
  "ggpubr",
  "hrbrthemes",
  "magrittr",
  "MotrpacHumanPreSuspension",
  "MotrpacHumanPreSuspensionAnalysis",
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

##########
data <- read.csv("precovid_muscle_ORA_newDA.csv", 
                 header = T, stringsAsFactors = F, 
                 check.names = F)
curated_terms <- read.csv("precovid_muscle_ORA_newDA_curated.csv", 
                          header = T, stringsAsFactors = F, 
                          check.names = F)
# term_order <- read.csv("precovid_muscle_ORA_newDA_curated_order.csv", 
#                        header = T, stringsAsFactors = F, 
#                        check.names = F)

ORA_selected <- data %>% 
  filter(set %in% curated_terms$set & 
           contrast %in% curated_terms$contrast)

plot_df1 <- ORA_selected %>%
  select(contrast, set_short, statistic_column, adj_p_value) %>%
  pivot_wider(id_cols = set_short,
              names_from = contrast, 
              values_from = c(adj_p_value, statistic_column),
              names_sep = "@") %>% 
  tidyr::pivot_longer(cols = -set_short,
                      names_to = c(".value", "contrast"),
                      names_sep = "@") %>%
  mutate(contrast = case_when(contrast == "ee_re_15min_only" ~ "EE_RE_15min",
                              contrast == "ee_re_3.5hr_only" ~ "EE_RE_3.5hr",
                              contrast == "ee_re_24hr_only" ~ "EE_RE_24hr")) %>%
  mutate(contrast = factor(contrast, 
                           levels = c("EE_RE_15min",
                                      "EE_RE_3.5hr",
                                      "EE_RE_24hr"))) %>%
  arrange(contrast) 
plot_df1$DE_status <- ifelse(plot_df1$adj_p_value < 0.05, "Significant", "Not significant")
plot_df1$set_short <- factor(plot_df1$set_short,
                             levels = rev(curated_terms$set_short))


# Dot plot
# jpeg(file = "precovid_muscle_figure2B_ORA_bubblePlot.jpg",
#      width = 6.5*1.25,height=8*1.25,units = "in",res=1200)
pdf(file = "precovid_muscle_figure2B_ORA_bubblePlot.pdf",
    width = 6.5*1.25,height=8*1.25)
ggplot(plot_df1, aes(
  x = contrast,
  y = set_short,
  size = statistic_column,
  fill = -log10(adj_p_value),
  color = DE_status                
)) +
  geom_point(
    shape = 21,
    stroke = 0.8
  ) +
  # scale_fill_viridis_c(
  #   option = "cividis",
  #   name = expression(-log[10]*"(adj p-value)"),
  #   direction = -1,
  #   guide = guide_colourbar(order = 1)
  # ) +
  scale_fill_gradient(
    name = expression(-log[10]*"(adj p-value)"),
    low = "white",
    high = "#57347d",
    guide = guide_colourbar(order = 1)
  ) +
  scale_color_manual(
    name = "Significance",
    values = c("Significant" = "black", 
               "Not significant" = "transparent")
  ) +
  scale_size_continuous(
    name = "Enrichment statistic",
    guide = guide_legend(override.aes = list(fill = "grey60"))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    x = "Contrast",
    y = "Pathway / Gene Set"
    #title = "Over-representation analysis (ORA) results"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, angle = 90,
                               hjust = 1, vjust = 0.5),
    legend.position = "right",
    axis.title = element_blank()
  ) +
  guides(color = "none") 

dev.off()
