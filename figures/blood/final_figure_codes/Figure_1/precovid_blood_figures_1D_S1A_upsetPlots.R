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
sel_tissue <- "blood"
sel_ome <- "all"
contrasts <- c("EE during 20 min",
               "EE during 40 min",
               "EE post 10 min",
               "EE post 30 min",
               "EE post 3.5 hr",
               "EE post 24 hr",
               "RE post 10 min",
               "RE post 30 min",
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
  mutate(contrast_short=case_when(contrast_short=="Endur.during_20_min" ~ "EE during 20 min",
                                  contrast_short=="Endur.during_40_min" ~ "EE during 40 min",
                                  contrast_short=="Endur.post_10_min" ~ "EE post 10 min",
                                  contrast_short=="Endur.post_15_30_45_min"~"EE post 30 min",
                                  contrast_short=="Endur.post_3.5_4_hr"~"EE post 3.5 hr",
                                  contrast_short=="Endur.post_24_hr"~"EE post 24 hr",
                                  contrast_short=="Resist.post_10_min" ~ "RE post 10 min",
                                  contrast_short=="Resist.post_15_30_45_min"~"RE post 30 min",
                                  contrast_short=="Resist.post_3.5_4_hr"~"RE post 3.5 hr",
                                  contrast_short=="Resist.post_24_hr"~"RE post 24 hr")) %>%
  mutate(ome=case_when(assay=="prot-ol"~"Proteomics",
                       assay=="prot-pr"~"Proteomics",
                       assay=="prot-ph"~"Phosphoproteomics",
                       assay=="transcript-rna-seq"~"Transcriptomics",
                       assay=="epigen-atac-seq"~"Chromatin Accessibility (ATAC)",
                       assay=="epigen-methylcap-seq"~"Methylation",
                       .default = "Metabolomics")) %>%
  mutate(randomGroupCode=case_when(randomGroupCode=="ADUEndur"~"EE",
                                   randomGroupCode=="ADUResist"~"RE",
                                   .default = randomGroupCode)) %>%
  mutate(ex_modality=gsub("_.*","",contrast_short),
         direction=ifelse(z.std>0,"Up","Down"),
         mode_direction=paste0(randomGroupCode,"_",direction)) %>%
  filter(ome!="Methylation")


# Option 1: Combine timepoints
# upset_data_all_features <- all_sig2 %>%
#        mutate(feature_modality = paste0(feature_id,"_",mode_direction)) %>%
#        distinct(feature_modality, .keep_all = T) %>%
#        mutate(feature_id2=feature_modality) %>%
#        select(c(feature_id,feature_id2,mode_direction)) %>%
#        pivot_wider(names_from = mode_direction,
#                                      values_from = feature_id2) %>%
#        mutate(across(-feature_id, ~ ifelse(is.na(.), 0, 1)))
upset_data_all_features <- all_sig2 %>%
       mutate(feature_modality = paste0(feature_id,"_",randomGroupCode)) %>%
       distinct(feature_modality, .keep_all = T) %>%
       mutate(feature_id2=feature_modality) %>%
       select(c(feature_id,feature_id2,randomGroupCode,direction)) %>%
       pivot_wider(names_from = randomGroupCode,
                                     values_from = feature_id2) %>%
       mutate(across(-c(feature_id,direction), ~ ifelse(is.na(.), 0, 1)))

upset_data_all_features2 <- inner_join(upset_data_all_features,
                                       all_sig2 %>%
                                         select(c(feature_id,ome)),
                                       by="feature_id") %>%
  distinct(feature_id, .keep_all = T) 


# methyl.upset <- upset_data_all_features2 %>%
#   filter(ome=="Methylation")
rna.upset <- upset_data_all_features2 %>%
  filter(ome=="Transcriptomics")
prot.ol.upset <- upset_data_all_features2 %>%
  filter(ome=="Proteomics")
metab.upset <- upset_data_all_features2 %>%
  filter(ome=="Metabolomics") 
 

pl1 <- ComplexUpset::upset(rna.upset %>%
                      select(-c(feature_id,ome)),
                    c("EE","RE"),
                    encode_sets=F,
                    width_ratio = 0.2,
                    height_ratio = 1,
                    set_sizes = F,
                    #min_size=10,
                    #min_degree=2,
                    base_annotations=list(
                      'Number of \ndifferential transcripts'=intersection_size(
                        counts=TRUE,
                        #fill="#4477AA"),
                        text=list(size=3),
                        mapping = aes(fill = direction)) +
                        scale_fill_manual(values = c(Down="#404688", 
                                                     Up="#999933"))+
                        labs(fill="Direction")),
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(axis.title.x=element_blank())
                      )
                    ),
                    sort_sets=FALSE,
                    wrap = TRUE,
                    stripes=c("#d95f02","#1b9e77")) +
  ggtitle("DA transcripts")

pl2 <- ComplexUpset::upset(prot.ol.upset %>%
                             select(-c(feature_id,ome)),
                           c("EE","RE"),
                           encode_sets=F,
                           width_ratio = 0.2,
                           height_ratio = 1,
                           set_sizes = F,
                           #min_size=10,
                           #min_degree=2,
                           base_annotations=list(
                             'Number of \ndifferential proteins'=intersection_size(
                               counts=TRUE,
                               #fill="#4477AA"),
                               text=list(size=3),
                               mapping = aes(fill = direction)) +
                               scale_fill_manual(values = c(Down="#404688", 
                                                            Up="#999933"))+
                               labs(fill="Direction")),
                           themes=upset_modify_themes(
                             list(
                               'intersections_matrix'=theme(axis.title.x=element_blank())
                             )
                           ),
                           sort_sets=FALSE,
                           wrap = TRUE,
                           stripes=c("#d95f02","#1b9e77")) +
  ggtitle("DA proteins")

pl3 <- ComplexUpset::upset(metab.upset %>%
                             select(-c(feature_id,ome)),
                           c("EE","RE"),
                           encode_sets=F,
                           width_ratio = 0.2,
                           height_ratio = 1,
                           set_sizes = F,
                           #min_size=10,
                           #min_degree=2,
                           base_annotations=list(
                             'Number of \ndifferential metabolites'=intersection_size(
                               counts=TRUE,
                               #fill="#4477AA"),
                               text=list(size=3),
                               mapping = aes(fill = direction)) +
                               scale_fill_manual(values = c(Down="#404688", 
                                                            Up="#999933"))+
                               labs(fill="Direction")),
                           themes=upset_modify_themes(
                             list(
                               'intersections_matrix'=theme(axis.title.x=element_blank())
                             )
                           ),
                           sort_sets=FALSE,
                           wrap = TRUE,
                           stripes=c("#d95f02","#1b9e77")) +
  ggtitle("DA metabolites")

# pl4 <- ComplexUpset::upset(methyl.upset %>%
#                              select(-c(feature_id,ome)),
#                            c("ADUEndur", "ADUResist"),
#                            encode_sets=F,
#                            width_ratio = 0.2,
#                            height_ratio = 1,
#                            set_sizes = F,
#                            #min_size=10,
#                            #min_degree=2,
#                            base_annotations=list(
#                              'Number of \ndifferential sites'=intersection_size(
#                                counts=TRUE,
#                                text=list(size=3),
#                                fill="#D687B5")),
#                            themes=upset_modify_themes(
#                              list(
#                                'intersections_matrix'=theme(axis.title.x=element_blank())
#                              )
#                            ),
#                            sort_sets=FALSE,
#                            wrap = TRUE,
#                            stripes=c("#d95f02","#1b9e77")) +
#   ggtitle("DA methylated regions")

pdf("precovid_blood_figure1D_upsetPlot.pdf",
    height = 4, width = 12)
#jpeg("upset_plot_blood_v3.jpeg", height = 4, width = 12, units = "in", res = 1200)
(pl1 | pl2 | pl3) +
  plot_layout(guides = "collect") + 
  theme(legend.position = "left")
dev.off()



################################################################################
# Option 2: Separated timepoints
upset_data_all_features_v2 <- all_sig2 %>%
  mutate(feature_id2=feature_id) %>%
  select(c(feature_id,feature_id2,contrast_short)) %>%
  pivot_wider(names_from = contrast_short,
              values_from = feature_id2) %>%
  mutate(across(-feature_id, ~ ifelse(is.na(.), 0, 1)))

upset_data_all_features2_v2 <- inner_join(upset_data_all_features_v2,
                                       all_sig2 %>%
                                         select(c(feature_id,ome)),
                                       by="feature_id") %>%
  distinct(feature_id, .keep_all = T) 

# methyl.upset2 <- upset_data_all_features2_v2 %>%
#   filter(ome=="Methylation")
rna.upset2 <- upset_data_all_features2_v2 %>%
  filter(ome=="Transcriptomics")
prot.ol.upset2 <- upset_data_all_features2_v2 %>%
  filter(ome=="Proteomics")
metab.upset2 <- upset_data_all_features2_v2 %>%
  filter(ome=="Metabolomics") 


pl5 <- ComplexUpset::upset(metab.upset2 %>%
                      select(-c(feature_id,ome)),
                    contrasts,
                    encode_sets=F,
                    width_ratio = 0.2,
                    height_ratio = 1,
                    set_sizes = F,
                    min_size=5,
                    #min_degree=2,
                    base_annotations=list(
                      'Number of \ndifferential metabolites'=intersection_size(
                        counts=TRUE,
                        text=list(size=3),
                        fill="#6D4B08")),
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(axis.title.x=element_blank())
                      )
                    ),
                    sort_sets=FALSE,
                    wrap = TRUE,
                    stripes=c("#FDE725", "#BAD071", "#D1BBD7", 
                              "#AE76A3", "#D1BBD7", "#AE76A3",
                              "#882E72", "#61194F")) +
  ggtitle("DA metabolites")

pl6 <- ComplexUpset::upset(rna.upset2 %>%
                             select(-c(feature_id,ome)),
                           contrasts,
                           encode_sets=F,
                           width_ratio = 0.2,
                           height_ratio = 1,
                           set_sizes = F,
                           min_size=50,
                           #min_degree=2,
                           base_annotations=list(
                             'Number of \ndifferential transcripts'=intersection_size(
                               counts=TRUE,
                               text=list(size=3),
                               fill="#4477AA")),
                           themes=upset_modify_themes(
                             list(
                               'intersections_matrix'=theme(axis.title.x=element_blank())
                             )
                           ),
                           sort_sets=FALSE,
                           wrap = TRUE,
                           stripes=c("#FDE725", "#BAD071", "#D1BBD7", 
                                     "#AE76A3", "#882E72", "#D1BBD7",
                                     "#AE76A3", "#882E72")) +
  ggtitle("DA transcripts")

pl7 <- ComplexUpset::upset(prot.ol.upset2 %>%
                             select(-c(feature_id,ome)),
                           contrasts,
                           encode_sets=F,
                           width_ratio = 0.2,
                           height_ratio = 1,
                           set_sizes = F,
                           min_size=3,
                           #min_degree=2,
                           base_annotations=list(
                             'Number of \ndifferential proteins'=intersection_size(
                               counts=TRUE,
                               text=list(size=3),
                               fill="#44AA99")),
                           themes=upset_modify_themes(
                             list(
                               'intersections_matrix'=theme(axis.title.x=element_blank())
                             )
                           ),
                           sort_sets=FALSE,
                           wrap = TRUE,
                           stripes=c("#FDE725", "#BAD071", "#D1BBD7", 
                                     "#882E72", "#D1BBD7", "#AE76A3",
                                     "#882E72")) +
  ggtitle("DA proteins")

# pl8 <- ComplexUpset::upset(methyl.upset2 %>%
#                              select(-c(feature_id,ome)),
#                            contrasts,
#                            encode_sets=F,
#                            width_ratio = 0.2,
#                            height_ratio = 1,
#                            set_sizes = F,
#                            min_size=5,
#                            #min_degree=2,
#                            base_annotations=list(
#                              'Number of \ndifferential sites'=intersection_size(
#                                counts=TRUE,
#                                text=list(size=3),
#                                fill="#D687B5")),
#                            themes=upset_modify_themes(
#                              list(
#                                'intersections_matrix'=theme(axis.title.x=element_blank())
#                              )
#                            ),
#                            sort_sets=FALSE,
#                            wrap = TRUE,
#                            stripes=c("#FDE725", "#882E72", "#61194F",
#                                      "#AE76A3", "#882E72", "#61194F")) +
#   ggtitle("DA methylation regions")


pdf("precovid_blood_figureS1A_upsetPlot.pdf",
    width = 7*1.25,height=12.05*1.25)
# jpeg(file = "precovid_blood_figure_upsetPlot.jpg",
#      width = 7*1.25,height=12.05*1.25,units = "in",res=1200)
pl6 / pl7 / pl5
dev.off()


