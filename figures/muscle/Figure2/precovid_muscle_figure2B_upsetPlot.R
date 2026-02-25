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


# jpeg(file = "precovid_muscle_figure2B_upsetPlot.jpg",
#      width = 9*1.25,height=4.05*1.25,units = "in",res=1200)
pdf(file = "precovid_muscle_figure2B_upsetPlot.pdf",
    width = 9.5*1.25,height=4*1.25)
ComplexUpset::upset(upset_data_all_features2 %>%
                      select(-c(feature_id)),
                    contrasts,
                    encode_sets=F,
                    width_ratio = 0.2,
                    height_ratio = 1,
                    set_sizes = F,
                    min_size=70,
                    #min_degree=2,
                    base_annotations=list(
                      'Number of \ndifferential features'=intersection_size(
                        counts=TRUE,
                        mapping=aes(fill=ome),
                        text=list(size=4,color="black")) +
                        scale_fill_manual(values=c(
                          'Metabolomics'=as.character("#6D4B08"), 
                          'Proteomics'=as.character("#228833"),
                          'Transcriptomics'=as.character("#4477AA"),
                          'Chromatin Accessibility (ATAC)'=as.character("#882255"),
                          'Phosphoproteomics'=as.character("#F3A02B"),
                          'Methylation'=as.character("#D687B5")))),
                    #fill="black")),
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(axis.title.x=element_blank())
                      )
                    ),
                    sort_sets=FALSE,
                    wrap = TRUE,
                    stripes = "grey95")
dev.off()
