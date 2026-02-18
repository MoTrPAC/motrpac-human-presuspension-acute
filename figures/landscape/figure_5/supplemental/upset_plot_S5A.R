#This code uses a configuration file to ensure local file structure and local machine settings can be used properly.

#' Steps for using the MoTrPAC config file:
#' 1. Open the configuration template file in the root folder of this repo called `example_config.json` available in `precovid-analyses/`
#' 2. Save a copy of this template to some other location in your computer. Many choose to place it in their home directory: Example: `~/motrpac_config.json`
#' 3. Customize the fields for your machine and save changes. Remove any "comments" as JSON does not permit comments
#' 4. This script expects the following elements in the config file:
##' a. "gsutil_user" and "gsutil_path": this is the way your machine calls the gsutil, usually "/usr/bin/gsutil"
##' b. "local_path": this is the absolute file path for the `precovid-analyses` repo. Please do not include the final "/"
#' 5. Add the path to this config file below 
# config <- config <- rjson::fromJSON(file = 'C:/Users/gayat/Documents/MoTrPAC/motrpac_precovid_config_GI.json')
rm(list=ls())
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
  "MotrpacHumanPreSuspensionData",
  "MotrpacHumanPreSuspension",
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
conflict_prefer("HUMAN_FEATURE_TO_GENE","MotrpacHumanPreSuspensionData")
conflict_prefer("CONTRAST_CONVERTER","MotrpacHumanPreSuspensionData")

repo_local_dir <- config$precovid_repo_path
gsutil_cmd <- config$gsutil_user

tissues<-c("Blood", "Muscle", "Adipose")
tissue_modality <- c("Muscle_EE","Muscle_RE",
                     "Blood_EE","Blood_RE",
                     "Adipose_EE","Adipose_RE")

all_da <- MotrpacHumanPreSuspensionData::load_differential_analysis(repo_local_dir = repo_local_dir, 
                                                                    selected_omes = "all",
                                                                    selected_tissues = "all",
                                                                    epigen = TRUE,
                                                                    combine_with_featgene = T,
                                                                    single_matrix = T) %>%
  filter(contrast_type=="exercise_with_controls") %>%
  mutate(contrast_short=gsub(" - .*","",contrast_short))

all_sig <- all_da %>%
  filter(adj_p_value < 0.05)

all_sig2 <- all_sig %>%
  mutate(contrast_short=case_when(contrast_short=="Endur.during_20_min"~"EE_during_20min",
                                  contrast_short=="Endur.during_40_min"~"EE_during_40min",
                                  contrast_short=="Endur.post_10_min"~"EE_post_10min",
                                  contrast_short=="Endur.post_15_30_45_min"~"EE_post_15_30_45_min",
                                  contrast_short=="Endur.post_3.5_4_hr"~"EE_post_3.5_4_hr",
                                  contrast_short=="Endur.post_24_hr"~"EE_post_24hr",
                                  contrast_short=="Resist.post_10_min"~"RE_post_10min",
                                  contrast_short=="Resist.post_15_30_45_min"~"RE_post_15_30_45_min",
                                  contrast_short=="Resist.post_3.5_4_hr"~"RE_post_3.5_4_hr",
                                  contrast_short=="Resist.post_24_hr"~"RE_post_24hr")) %>%
  mutate(Ome=case_when(assay=="prot-ol"~"Proteomics",
                       assay=="prot-pr"~"Proteomics",
                       assay=="prot-ph"~"Phosphoproteomics",
                       assay=="transcript-rna-seq"~"Transcriptomics",
                       assay=="epigen-atac-seq"~"Chromatin Accessibility (ATAC)",
                       assay=="epigen-methylcap-seq"~"Methylation",
                       .default = "Metabolomics")) %>%
  mutate(tissue=case_when(tissue=="adipose"~"Adipose",
                          tissue=="blood"~"Blood",
                          tissue=="muscle"~"Muscle")) %>%
  mutate(ex_modality=gsub("_.*","",contrast_short))

all_sig3 <- all_sig2 %>%
  mutate(feature_id = case_when(Ome == "Proteomics" ~ uniprot,
                                .default = feature_id))

all_sig4 <- all_sig3 %>%
  mutate(tissue_modality = paste0(tissue,"_",ex_modality))

################################################################################
# upset plot of DA features across tissues and ex modalities #
upset_data_features <- all_sig4 %>%
  mutate(feature_tissue = paste0(feature_id,"_",tissue_modality)) %>%
  distinct(feature_tissue, .keep_all = T) %>%
  mutate(feature_id2=feature_tissue) %>%
  select(c(feature_id,feature_id2,tissue_modality)) %>%
  pivot_wider(names_from = tissue_modality,
              values_from = feature_id2) %>%
  mutate(across(-feature_id, ~ ifelse(is.na(.), 0, 1))) 

upset_data_features2 <- inner_join(upset_data_features,
                                   all_sig4 %>%
                                     select(c(feature_id,Ome)),
                                   by="feature_id") %>%
  distinct(feature_id, .keep_all = T)

## upset plot ##
# jpeg(file = "landscape_figureS3A_upsetPlot_by_tissue_modality_features_newDA.jpg",
#      width = 10*1.25,height=5*1.25,units = "in",res=1200)
pdf(file = "landscape_figureS3A_upsetPlot_by_tissue_modality_features_newDA.pdf",
     width = 10*1.25,height=5*1.25)
ComplexUpset::upset(upset_data_features2 %>%
                      select(-c(feature_id)),
                    tissue_modality,
                    encode_sets=F,
                    width_ratio = 0.2,
                    height_ratio = 1,
                    set_sizes = F,
                    min_size=7,
                    #min_degree=2,
                    base_annotations=list(
                      'Number of \ndifferential features'=intersection_size(
                        counts=TRUE,
                        mapping=aes(fill=Ome),
                        text=list(size=3)) +
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
                    stripes=c(as.character(HUMAN_TISSUE_COLORS["muscle"]),
                              as.character(HUMAN_TISSUE_COLORS["muscle"]),
                              as.character(HUMAN_TISSUE_COLORS["blood"]),
                              as.character(HUMAN_TISSUE_COLORS["blood"]),
                              as.character(HUMAN_TISSUE_COLORS["adipose"]),
                              as.character(HUMAN_TISSUE_COLORS["adipose"])))
dev.off()