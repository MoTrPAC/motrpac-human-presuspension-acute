rm(list=ls())
session_packages = c(
  "ComplexHeatmap",
  "ComplexUpset",
  "conflicted",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "magrittr",
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


# repo_local_dir <- config$precovid_repo_path
# gsutil_cmd <- config$gsutil_user

tissues<-c("Blood", "Muscle", "Adipose")
tissue_modality <- c("Muscle_EE","Muscle_RE",
                     "Blood_EE","Blood_RE",
                     "Adipose_EE","Adipose_RE")

all_da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = "all",
                                                                        selected_tissues = "all",
                                                                        epigen = TRUE,
                                                                        combine_with_featgene = TRUE,
                                                                        single_matrix = TRUE) %>%
  filter(contrast_type=="exercise_with_controls") %>%
  mutate(contrast_short=gsub(" - .*","",contrast_short))

all_sig <- all_da %>%
  filter(adj_p_value < 0.05)

all_sig2 <- all_sig %>%
  mutate(contrast_short = contrast_short %>%
           gsub("^Endur", "EE", .) %>%
           gsub("^Resist", "RE", .) %>%
           gsub("\\.", "_", .)) %>%
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


pdf(file = "landscape_figure2B_upsetPlot_features_newDA.pdf",
    width = 6*1.25,height=4.05*1.25)
#jpeg(file = "landscape_figure2B_upsetPlot_features_newDA.jpg",
#     width = 6*1.25,height=4.05*1.25,units = "in",res=1200)
## upset plot ##
p = ComplexUpset::upset(upset_data_all_features2 %>%
                          select(-c(feature_id)),
                        tissues,
                        encode_sets=F,
                        width_ratio = 0.2,
                        height_ratio = 1,
                        set_sizes = FALSE,
                        #min_degree=2,
                        base_annotations=list(
                          'Number of \ndifferential features'=intersection_size(
                            counts=TRUE,
                            mapping=aes(fill=Ome),
                            text=list(size=3)) +
                            scale_fill_manual(values=HUMAN_OME_COLORS)),
                        #fill="black")),
                        themes=upset_modify_themes(
                          list(
                            'intersections_matrix'= theme(axis.title.x=element_blank())
                          )
                        ),
                        sort_sets="ascending",
                        wrap = TRUE,
                        stripes=c(as.character(HUMAN_TISSUE_COLORS["muscle"]),
                                  as.character(HUMAN_TISSUE_COLORS["blood"]),
                                  as.character(HUMAN_TISSUE_COLORS["adipose"])))



dev.off()


################################################################################
## Run ORA ##
# get background genes for ORA #
blood_bg_genes <- all_da %>%
  filter(tissue == "blood",
         !is.na(gene_symbol)) %>%
  distinct(feature_id, .keep_all = T) %>%
  select(feature_id,gene_symbol,assay)

mus_bg_genes <- all_da %>%
  filter(tissue == "muscle",
         !is.na(gene_symbol)) %>%
  distinct(feature_id, .keep_all = T) %>%
  select(feature_id,gene_symbol,assay)

adi_bg_genes <- all_da %>%
  filter(tissue == "adipose",
         !is.na(gene_symbol)) %>%
  distinct(feature_id, .keep_all = T) %>%
  select(feature_id,gene_symbol,assay)

## ORA ##
# muscle only
mus_only <- upset_data_all_features2 %>%
  filter(Blood==0 &
           Muscle==1 &
           Adipose==0)
mus_only2 <- inner_join(mus_only %>%
                          select(feature_id,Ome),
                        HUMAN_FEATURE_TO_GENE %>%
                          select(feature_id,
                                 gene_symbol),
                        by= "feature_id")
mus_only_ora <- run_ORA(input = as.character(mus_only2$gene_symbol),
                        background = as.character(mus_bg_genes$gene_symbol),
                        database = names(MOLECULAR_SIGNATURES)) %>%
  mutate(statistic_column = -log(p_value))

# blood only
blo_only <- upset_data_all_features2 %>%
  filter(Blood==1 &
           Muscle==0 &
           Adipose==0)
blo_only2 <- inner_join(blo_only %>%
                          select(feature_id,Ome),
                        HUMAN_FEATURE_TO_GENE %>%
                          select(feature_id,
                                 gene_symbol),
                        by= "feature_id")
blo_only_ora <- run_ORA(input = as.character(blo_only2$gene_symbol),
                        background = as.character(blood_bg_genes$gene_symbol),
                        database = names(MOLECULAR_SIGNATURES)) %>%
  mutate(statistic_column = -log(p_value))


# adipose only
adi_only <- upset_data_all_features2 %>%
  filter(Blood==0 &
           Muscle==0 &
           Adipose==1)
adi_only2 <- inner_join(adi_only %>%
                          select(feature_id,Ome),
                        HUMAN_FEATURE_TO_GENE %>%
                          select(feature_id,
                                 gene_symbol),
                        by= "feature_id")
adi_only_ora <- run_ORA(input = as.character(adi_only2$gene_symbol),
                        background = as.character(adi_bg_genes$gene_symbol),
                        database = names(MOLECULAR_SIGNATURES)) %>%
  mutate(statistic_column = -log(p_value))

# blood-muscle overlap
blo_mus <- upset_data_all_features2 %>%
  filter(Blood==1 &
           Muscle==1 &
           Adipose==0)
blo_mus2 <- inner_join(blo_mus %>%
                         select(feature_id,Ome),
                       HUMAN_FEATURE_TO_GENE %>%
                         select(feature_id,
                                gene_symbol),
                       by= "feature_id")
blo_mus_ora <- run_ORA(input = as.character(blo_mus2$gene_symbol),
                       background = as.character(intersect(blood_bg_genes$gene_symbol,
                                                           mus_bg_genes$gene_symbol)),
                       database = names(MOLECULAR_SIGNATURES)) %>%
  mutate(statistic_column = -log(p_value))

# adipose-muscle overlap
adi_mus <- upset_data_all_features2 %>%
  filter(Blood==0 &
           Muscle==1 &
           Adipose==1)
adi_mus2 <- inner_join(adi_mus %>%
                         select(feature_id,Ome),
                       HUMAN_FEATURE_TO_GENE %>%
                         select(feature_id,
                                gene_symbol),
                       by= "feature_id")
adi_mus_ora <- run_ORA(input = as.character(adi_mus2$gene_symbol),
                       background = as.character(intersect(adi_bg_genes$gene_symbol,
                                                           mus_bg_genes$gene_symbol)),
                       database = names(MOLECULAR_SIGNATURES)) %>%
  mutate(statistic_column = -log(p_value))


# adipose-muscle-blood overlap
adi_blo_mus <- upset_data_all_features2 %>%
  filter(Blood==1 &
           Muscle==1 &
           Adipose==1)
adi_blo_mus2 <- inner_join(adi_blo_mus %>%
                             select(feature_id,Ome),
                           HUMAN_FEATURE_TO_GENE %>%
                             select(feature_id,
                                    gene_symbol),
                           by= "feature_id")
adi_blo_mus_ora <- run_ORA(input = as.character(adi_blo_mus2$gene_symbol),
                           background = as.character(Reduce(intersect,
                                                            list(adi_bg_genes$gene_symbol,
                                                                 mus_bg_genes$gene_symbol,
                                                                 blood_bg_genes$gene_symbol))),
                           database = names(MOLECULAR_SIGNATURES)) %>%
  mutate(statistic_column = -log(p_value))


all_ORA_res <- c()
all_ORA_res[["Muscle_only"]] <- mus_only_ora
all_ORA_res[["Blood_only"]] <- blo_only_ora
all_ORA_res[["Blood_and_Muscle_only"]] <- blo_mus_ora
all_ORA_res[["Adipose_only"]] <- adi_only_ora
all_ORA_res[["Adipose_and_Muscle_only"]] <- adi_mus_ora
all_ORA_res[["Adipose_Blood_Muscle"]] <- adi_blo_mus_ora

all_ORA_res2 <- dplyr::bind_rows(all_ORA_res,
                                 .id = "contrast")
all_ORA_res2$contrast <- factor(all_ORA_res2$contrast,
                                levels = c("Muscle_only",
                                           "Blood_only",
                                           "Blood_and_Muscle_only",
                                           "Adipose_only",
                                           "Adipose_and_Muscle_only",
                                           "Adipose_Blood_Muscle"))

####
# ORA bar plot from the 50 features DA in all three tissues
da_3tissues_pathways = readxl::read_excel(path = file.path(here(), "figures", "landscape", "curated_pathways.xlsx"),
                                          sheet = "2B") %>%
  mutate(contrast="Adipose_Blood_Muscle") %>%
  select(contrast,everything())

# jpeg(file = "landscape_figureS2_DA_features_3tissues_ORA_heatmap_newDA_newORA.jpg",
#      width = 5*1.25,height=5*1.25,units = "in",res=1200)
pdf(file = "landscape_figure2B_DA_features_3tissues_ORA_barPlot_newDA_newORA.pdf",
    width = 5*1.25,height=5*1.25)
da_3tissues_pathways %>%
  arrange(statistic_column) %>%
  mutate(set_short=sub("^[^_]*_", "",set_short)) %>%
  mutate(set_short=gsub("_"," ",set_short)) %>%
  mutate(set_short=toupper(set_short)) %>%
  mutate(set_short=factor(set_short,levels=set_short)) %>%
  ggplot(aes(x=set_short, y=statistic_column)) +
  geom_bar(stat="identity",
           width=0.9,
           fill="bisque3") +
  geom_text(
    aes(label = set_short),
    hjust = 1.1,
    color = "grey20",
    size = 3
  ) +
  coord_flip() +
  #scale_y_reverse() +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("") +
  ylab("-log(p)")
dev.off()
