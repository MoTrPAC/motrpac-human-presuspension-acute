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

all_da <- MotrpacHumanPreSuspensionData::load_differential_analysis(selected_omes = "all",
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
# Heatmap of shared phosphosites between Adipose and Muscle

adi_mus_phospho <- upset_data_all_features2 %>%
  filter(Adipose==1 & Muscle==1 &
           Ome=="Phosphoproteomics")
adi_mus_phospho <- inner_join(adi_mus_phospho,
                              HUMAN_FEATURE_TO_GENE %>%
                                select(feature_id,gene_symbol),
                              by="feature_id")

extract_phosphosite <- function(x) {
  sub(".*_(.+).{1}$", "\\1", x)
}

adi_mus_phospho <- adi_mus_phospho %>%
  mutate(feature_id2 = paste0(gene_symbol," ",
                              extract_phosphosite(feature_id)))

features_to_plot <- adi_mus_phospho$feature_id
sel_genes_da <- all_da %>%
  filter(feature_id %in% features_to_plot) %>%
  separate(contrast_short,into = c("Modality","Timepoint"),
           remove = FALSE, extra = "merge") %>%
  mutate(gene_symbol = case_when(is.na(gene_symbol) ~ feature_id,
                                 .default = gene_symbol)) %>%
  mutate(Timepoint = case_when (Timepoint == "post_15_30_45_min" ~ "post 15/30/45 min",
                                Timepoint == "post_3.5_4_hr" ~ "post 3.5/4 hr",
                                .default = Timepoint)) %>%
  mutate(Timepoint = gsub("_", " ", Timepoint)) %>%
  filter(tissue!="blood") %>%
  mutate(Modality=case_when(Modality=="Endur"~"EE",
                            Modality=="Resist"~"RE"),
         tissue=str_to_title(tissue))

sel_genes_da2 <- inner_join(sel_genes_da %>%
                              select(feature_id,contrast_short,tissue,z.std) %>%
                              pivot_wider(names_from = c(tissue,contrast_short),
                                          values_from = z.std),
                            HUMAN_FEATURE_TO_GENE %>% select(feature_id,gene_symbol,
                                                             assay),
                            by = "feature_id") %>%
  mutate(gene_symbol = case_when(is.na(gene_symbol) ~ feature_id,
                                 .default = gene_symbol),
         feature_id2 = paste0("ph-",gene_symbol," ",
                              extract_phosphosite(feature_id))) %>%
  select(-c(feature_id,gene_symbol)) %>%
  column_to_rownames(var="feature_id2") %>%
  replace(is.na(.),0)

adj_p_mat <- inner_join(sel_genes_da %>%
                          select(feature_id,contrast_short,tissue,adj_p_value) %>%
                          pivot_wider(names_from = c(tissue,contrast_short),
                                      values_from = adj_p_value),
                        HUMAN_FEATURE_TO_GENE %>% select(feature_id,gene_symbol,
                                                         assay),
                        by = "feature_id") %>%
  mutate(gene_symbol = case_when(is.na(gene_symbol) ~ feature_id,
                                 .default = gene_symbol),
         feature_id2 = paste0("ph-",gene_symbol," ",
                              extract_phosphosite(feature_id))) %>%
  select(-c(feature_id,gene_symbol)) %>%
  column_to_rownames(var="feature_id2") %>%
  replace(is.na(.),1)

col.anno <- sel_genes_da %>%
  distinct(tissue,contrast_short,
           .keep_all = T) %>%
  select(tissue,contrast_short,Modality,Timepoint) %>%
  mutate(Timepoint=factor(Timepoint, levels = c("during 20 min",
                                                "during 40 min",
                                                "post 10 min",
                                                "post 15/30/45 min",
                                                "post 3.5/4 hr",
                                                "post 24 hr"))) %>%
  group_by(tissue,Modality,Timepoint) %>%
  arrange(tissue,Modality,Timepoint, .by_group = T) %>%
  mutate(tissue_contrast = paste0(tissue,"_",contrast_short)) %>%
  column_to_rownames(var="tissue_contrast")
sel_genes_da2 <- sel_genes_da2[,rownames(col.anno)]
adj_p_mat <- adj_p_mat[,rownames(col.anno)]

column_ha <- columnAnnotation(Tissue = col.anno$tissue,
                              Modality = col.anno$Modality,
                              Timepoint = col.anno$Timepoint,
                              col = list(Modality=c("EE"="#d95f02",
                                                    "RE"="#1b9e77"),
                                         Timepoint=c("post 10 min"="#D1BBD7",
                                                     "post 15/30/45 min"="#AE76A3",
                                                     "post 3.5/4 hr"="#882E72",
                                                     "post 24 hr"="#61194F"),
                                         Tissue=c("Adipose"="#ffffbf",
                                                  "Muscle"="#abd9e9")))

# jpeg(file = "landscape_figure2D_adi_mus_shared_phospho_newDA.jpg",
#      width = 5*1.25,height=7*1.25,units = "in",res=1200)
pdf(file = "landscape_figure2D_adi_mus_shared_phospho_newDA.pdf",
     width = 5*1.25,height=7*1.25)
ComplexHeatmap::Heatmap(sel_genes_da2,
                        top_annotation = column_ha,
                        #right_annotation = row_ha,
                        border = T,
                        cluster_rows = T,
                        cluster_columns = F,
                        column_split = c(rep("Adipose",2),
                                         rep("Muscle",6)),
                        column_gap = unit(2, "mm"),
                        row_names_gp = gpar(fontsize=12),
                        column_names_gp = gpar(fontsize=10),
                        rect_gp = gpar(col = "#2F4F4F",
                                       lwd = 1),
                        show_column_names = F,
                        name = "z.std",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if (adj_p_mat[i, j] < 0.05) {
                            grid.points(
                              x, y,
                              pch = 16,
                              size = unit(5, "pt"),
                              gp = gpar(col = "grey20")
                            )
                          }
                        }
                      )
dev.off()
