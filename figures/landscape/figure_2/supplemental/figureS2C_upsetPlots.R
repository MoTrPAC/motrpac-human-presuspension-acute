rm(list=ls())
session_packages = c(
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "hrbrthemes",
  "magrittr",
  "MotrpacHumanPreSuspensionAnalysis",
  "patchwork",
  "purrr",
  "RefMet",
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


tissues<-c("Blood", "Muscle", "Adipose")
tissue_modality <- c("Muscle_EE","Muscle_RE",
                     "Blood_EE","Blood_RE",
                     "Adipose_EE","Adipose_RE")

all_da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis( selected_omes = "all",
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
# upset plot of DA features across tissues #
upset_data_all_features <- all_sig4 %>%
  mutate(feature_tissue = paste0(feature_id,"_",tissue)) %>%
  distinct(feature_tissue, .keep_all = T) %>%
  mutate(feature_id2=feature_id) %>%
  select(c(feature_id,feature_id2,tissue)) %>%
  pivot_wider(names_from = tissue,
              values_from = feature_id2) %>%
  mutate(across(-feature_id, ~ ifelse(is.na(.), 0, 1)))

upset_data_all_features2 <- inner_join(upset_data_all_features,
                                       all_sig4 %>%
                                         select(c(feature_id,Ome)),
                                       by="feature_id") %>%
  distinct(feature_id, .keep_all = T)

################################################################################
## ome-wise upset plot

# metab
metab.upset <- upset_data_all_features2 %>%
  filter(Ome=="Metabolomics")
RefMet_mapped <- refmet_map_df(metab.upset$feature_id)
metab.upset2 <- inner_join(metab.upset,
                           RefMet_mapped %>%
                             select(Input.name,
                                    Standardized.name,
                                    Super.class,
                                    Main.class,
                                    Sub.class) %>%
                             rename(feature_id=Input.name),
                           by="feature_id") %>%
  mutate(Super.class=case_when(Super.class=="-" ~ "Sphingolipids",
                               .default = Super.class),
         Main.class=case_when(Main.class=="-" ~ "Ceramides",
                              .default = Main.class),
         Sub.class=case_when(Sub.class=="-" ~ "Cer",
                             .default = Sub.class)) %>%
  mutate(Super.class2=case_when(Super.class %in% c("Fatty Acyls",
                                                   "Sphingolipids",
                                                   "Glycerophospholipids",
                                                   "Prenol Lipids",
                                                   "Sterol Lipids",
                                                   "Glycerolipids") ~ "Lipids",
                                Super.class %in% c("Organic nitrogen compounds",
                                                   "Alkaloids",
                                                   "Nucleic acids") ~ "Organic Nitrogen Compounds",
                                Super.class %in% c("Benzenoids",
                                                   "Organoheterocyclic compounds") ~ "Aromatic Compounds",
                                .default = Super.class))
refmet_class_cols <- c("#E69F00","#56B4E9","#808000",
                       "#F0E442","#CC79A7")
names(refmet_class_cols) <- unique(metab.upset2$Super.class2)


ComplexUpset::upset(metab.upset2 %>%
                      select(-c(feature_id,
                                Main.class,
                                Sub.class,
                                Super.class)),
                    tissues,
                    encode_sets=F,
                    width_ratio = 0.2,
                    height_ratio = 1,
                    set_sizes = F,
                    #min_size=10,
                    #min_degree=2,
                    base_annotations=list(
                      'Number of \ndifferential metabolites'=intersection_size(
                        counts=TRUE,
                        mapping=aes(fill=Super.class2),
                        text=list(size=3)) +
                        scale_fill_manual(values=refmet_class_cols)+
                        labs(fill="RefMet class")),
                    #fill="black")),
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(axis.title.x=element_blank())
                      )
                    ),
                    sort_sets=FALSE,
                    wrap = TRUE,
                    stripes=c(as.character(HUMAN_TISSUE_COLORS["blood"]),
                              as.character(HUMAN_TISSUE_COLORS["muscle"]),
                              as.character(HUMAN_TISSUE_COLORS["adipose"]))) +
  ggtitle("DA metabolites")

###################

metab.upset <- ComplexUpset::upset(upset_data_all_features2 %>%
                                     filter(Ome == "Metabolomics") %>%
                                     select(-c(feature_id)),
                                   tissues,
                                   encode_sets=F,
                                   width_ratio = 0.2,
                                   height_ratio = 1,
                                   set_sizes = F,
                                   #min_size=10,
                                   #min_degree=2,
                                   base_annotations=list(
                                     'Number of \ndifferential metabolites'=intersection_size(
                                       counts=TRUE,
                                       mapping=aes(fill=Ome),
                                       text=list(size=3),
                                       fill="#6D4B08")),
                                   #fill="black")),
                                   themes=upset_modify_themes(
                                     list(
                                       'intersections_matrix'=theme(axis.title.x=element_blank())
                                     )
                                   ),
                                   sort_sets=FALSE,
                                   wrap = TRUE,
                                   stripes=c(as.character(HUMAN_TISSUE_COLORS["blood"]),
                                             as.character(HUMAN_TISSUE_COLORS["muscle"]),
                                             as.character(HUMAN_TISSUE_COLORS["adipose"]))) +
  ggtitle("DA metabolites")


prot.upset <- ComplexUpset::upset(upset_data_all_features2 %>%
                                    filter(Ome == "Proteomics") %>%
                                    select(-c(feature_id)),
                                  tissues,
                                  encode_sets=F,
                                  width_ratio = 0.2,
                                  height_ratio = 1,
                                  set_sizes = F,
                                  #min_degree=2,
                                  base_annotations=list(
                                    'Number of differential \nproteins'=intersection_size(
                                      counts=TRUE,
                                      mapping=aes(fill=Ome),
                                      text=list(size=3),
                                      fill="#228833")),
                                  themes=upset_modify_themes(
                                    list(
                                      'intersections_matrix'=theme(axis.title.x=element_blank())
                                    )
                                  ),
                                  sort_sets="ascending",
                                  wrap = TRUE,
                                  stripes=c(as.character(HUMAN_TISSUE_COLORS["adipose"]),
                                            as.character(HUMAN_TISSUE_COLORS["blood"]),
                                            as.character(HUMAN_TISSUE_COLORS["muscle"]))) +
  ggtitle("DA Proteins")

protph.upset <- ComplexUpset::upset(upset_data_all_features2 %>%
                                      filter(Ome=="Phosphoproteomics") %>%
                                      select(-c(feature_id)),
                                    tissues,
                                    encode_sets=F,
                                    width_ratio = 0.2,
                                    height_ratio = 1,
                                    set_sizes = F,
                                    #min_degree=2,
                                    base_annotations=list(
                                      'Number of differential \nphosphosites'=intersection_size(
                                        counts=TRUE,
                                        text=list(size=3),
                                        fill="#FF7F00")),
                                    themes=upset_modify_themes(
                                      list(
                                        'intersections_matrix'=theme(axis.title.x=element_blank())
                                      )
                                    ),
                                    sort_sets="ascending",
                                    wrap = TRUE,
                                    stripes=c(as.character(HUMAN_TISSUE_COLORS["muscle"]),
                                              as.character(HUMAN_TISSUE_COLORS["adipose"]),
                                              as.character(HUMAN_TISSUE_COLORS["blood"]))) +
  ggtitle("DA Phosphosites")

atac.upset <- ComplexUpset::upset(upset_data_all_features2 %>%
                                    filter(Ome=="Chromatin Accessibility (ATAC)") %>%
                                    select(-c(feature_id)),
                                  tissues,
                                  encode_sets=F,
                                  width_ratio = 0.2,
                                  height_ratio = 1,
                                  set_sizes = F,
                                  #min_degree=2,
                                  base_annotations=list(
                                    'Number of differential \npeaks'=intersection_size(
                                      counts=TRUE,
                                      text=list(size=3),
                                      fill="#882255")),
                                  themes=upset_modify_themes(
                                    list(
                                      'intersections_matrix'=theme(axis.title.x=element_blank())
                                    )
                                  ),
                                  sort_sets="ascending",
                                  wrap = TRUE,
                                  stripes=c(as.character(HUMAN_TISSUE_COLORS["muscle"]),
                                            as.character(HUMAN_TISSUE_COLORS["blood"]),
                                            as.character(HUMAN_TISSUE_COLORS["adipose"]))) +
  ggtitle("Differentially Accessible Peaks")

rna.upset <- ComplexUpset::upset(upset_data_all_features2 %>%
                                   filter(Ome=="Transcriptomics") %>%
                                   select(-c(feature_id)),
                                 tissues,
                                 encode_sets=F,
                                 width_ratio = 0.2,
                                 height_ratio = 1,
                                 set_sizes = F,
                                 #min_degree=2,
                                 base_annotations=list(
                                   'Number of differential \ntranscripts'=intersection_size(
                                     counts=TRUE,
                                     text=list(size=3),
                                     fill="#377EB8")),
                                 themes=upset_modify_themes(
                                   list(
                                     'intersections_matrix'=theme(axis.title.x=element_blank())
                                   )
                                 ),
                                 sort_sets="ascending",
                                 wrap = TRUE,
                                 stripes=c(as.character(HUMAN_TISSUE_COLORS["muscle"]),
                                           as.character(HUMAN_TISSUE_COLORS["blood"]),
                                           as.character(HUMAN_TISSUE_COLORS["adipose"]))) +
  ggtitle("DA transcripts")

methylcap.upset <- ComplexUpset::upset(upset_data_all_features2 %>%
                                         filter(Ome=="Methylation") %>%
                                         select(-c(feature_id)),
                                       tissues,
                                       encode_sets=F,
                                       width_ratio = 0.2,
                                       height_ratio = 1,
                                       set_sizes = F,
                                       #min_degree=2,
                                       base_annotations=list(
                                         'Number of differential \npeaks'=intersection_size(
                                           counts=TRUE,
                                           text=list(size=3),
                                           fill="#D687B5")),
                                       themes=upset_modify_themes(
                                         list(
                                           'intersections_matrix'=theme(axis.title.x=element_blank())
                                         )
                                       ),
                                       sort_sets="ascending",
                                       wrap = TRUE,
                                       stripes=c(as.character(HUMAN_TISSUE_COLORS["blood"]),
                                                 as.character(HUMAN_TISSUE_COLORS["adipose"]),
                                                 as.character(HUMAN_TISSUE_COLORS["muscle"]))) +
  ggtitle("DA methylated regions")

combined_single_ome_upset_plot <- (atac.upset | methylcap.upset | rna.upset) /
  (prot.upset | protph.upset | metab.upset)

pdf(file = "landscape_figureS2C_ome_upsetPlot_newDA.pdf",
    width = 7*1.25,height=5*1.25)
combined_single_ome_upset_plot
dev.off()
