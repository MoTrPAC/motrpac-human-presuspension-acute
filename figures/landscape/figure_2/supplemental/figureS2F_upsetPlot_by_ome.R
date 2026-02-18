rm(list=ls())

session_packages = c(
  "MotrpacHumanPreSuspensionAnalysis",
  "tidyverse",
  "grid",
  "ComplexHeatmap",
  "UpsetR",
  "conflicted",
  "ComplexUpset"
)

# ----complexupset is giving problems with ggplot 4.0.0----
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    packageVersion("ggplot2") != "3.5.2") {
  stop("ggplot2 3.5.2 is required. Please install with:
       remotes::install_version('ggplot2', version = '3.5.2')")
}


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


####
tissues<-c("Blood", "Muscle", "Adipose")
sel_group <- c("Endur", "Resist")
human_omic_cols <- HUMAN_OME_COLORS

adipose_col_seq <- adipose_col_seq_ee <- c(as.character(human_omic_cols["Methylation"]),
                                           as.character(human_omic_cols["Transcriptomics"]),
                                           as.character(human_omic_cols["Phosphoproteomics"]),
                                           as.character(human_omic_cols["Proteomics"]))

adipose_col_seq_re <- c(as.character(human_omic_cols["Proteomics"]),
                        as.character(human_omic_cols["Methylation"]),
                        as.character(human_omic_cols["Phosphoproteomics"]),
                        as.character(human_omic_cols["Transcriptomics"]))


muscle_col_seq <- c(as.character(human_omic_cols["Proteomics"]),
                    as.character(human_omic_cols["Methylation"]),
                    as.character(human_omic_cols["ATAC"]),
                    as.character(human_omic_cols["Phosphoproteomics"]),
                    as.character(human_omic_cols["Transcriptomics"]))

blood_col_seq_ee <- c(as.character(human_omic_cols["Methylation"]),
                      as.character(human_omic_cols["Proteomics (Olink)"]),
                      as.character(human_omic_cols["Transcriptomics"]))

blood_col_seq <- blood_col_seq_re <- c(as.character(human_omic_cols["Proteomics (Olink)"]),
                                       as.character(human_omic_cols["Methylation"]),
                                       as.character(human_omic_cols["Transcriptomics"]))


#####

data <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_tissues = "all",
                                                                      single_matrix = F,
                                                                      epigen = T,
                                                                      combine_with_featgene = T)
prot_sig <- data.frame()

for(i in c('muscle', 'adipose')){

  prot_pr <- data[[i]][['prot-pr']] %>%
    dplyr::filter(adj_p_value < 0.05 &
                    contrast_type=="exercise_with_controls") %>%
    select(feature_id, logFC, contrast_short,tissue)
  prot_sig <- rbind(prot_sig, prot_pr)

}

ph_sig <- data.frame()

for(i in c('muscle', 'adipose')){

  ph <- data[[i]][['prot-ph']] %>%
    dplyr::filter(adj_p_value < 0.05 &
                    contrast_type=="exercise_with_controls") %>%
    select(feature_id, logFC, contrast_short,tissue)
  ph_sig <- rbind(ph_sig, ph)

}

prot_ol_sig <- data.frame()

for(i in c('blood')){

  prot_ol <- data[[i]][['prot-ol']] %>%
    dplyr::filter(adj_p_value < 0.05 &
                    contrast_type=="exercise_with_controls") %>%
    select(feature_id, logFC, contrast_short,tissue)
  prot_ol_sig <- rbind(prot_ol_sig, prot_ol)

}

transcript_sig <- data.frame()

for(i in c('muscle', 'adipose', 'blood')){

  transcript <- data[[i]][['transcript-rna-seq']] %>%
    dplyr::filter(adj_p_value < 0.05 &
                    contrast_type=="exercise_with_controls") %>%
    select(feature_id, logFC, contrast_short,tissue)
  transcript_sig <- rbind(transcript_sig, transcript)

}

atac_sig <- data.frame()

for(i in c('muscle', 'blood')){

  atac <- data[[i]][['epigen-atac-seq']] %>%
    dplyr::filter(adj_p_value < 0.05 &
                    contrast_type=="exercise_with_controls") %>%
    select(feature_id, logFC, contrast_short,tissue)
  atac_sig <- rbind(atac_sig, atac)

}

methyl_sig <- data.frame()

for(i in c('muscle', 'blood','adipose')){

  methyl <- data[[i]][['epigen-methylcap-seq']] %>%
    dplyr::filter(adj_p_value < 0.05 &
                    contrast_type=="exercise_with_controls") %>%
    select(feature_id, t, contrast_short,tissue)
  methyl_sig <- rbind(methyl_sig, methyl)

}


# function make upset plot for each tissue
make_upset_plot <- function(
    sel_tissue,
    HUMAN_FEATURE_TO_GENE,
    prot_sig,
    ph_sig,
    transcript_sig,
    prot_ol_sig,
    atac_sig,
    methyl_sig,
    sel_group,
    col_seq,
    sel_tissue_col
) {
  pr_all <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(assay =='prot-pr') %>%
    dplyr::filter(feature_id %in% (prot_sig %>%
                                     filter(
                                       tissue==sel_tissue &
                                         grepl(paste(sel_group, collapse = "|"), contrast_short)
                                     ) %>%
                                     select(feature_id) %>%
                                     pull()))

  ph_all <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(assay =='prot-ph') %>%
    dplyr::filter(feature_id %in% (ph_sig %>%
                                     filter(
                                       tissue==sel_tissue &
                                         grepl(paste(sel_group, collapse = "|"), contrast_short)
                                     ) %>%
                                     select(feature_id) %>%
                                     pull()))

  tr_all <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(assay =='transcript-rna-seq') %>%
    dplyr::filter(feature_id %in% (transcript_sig %>%
                                     filter(
                                       tissue==sel_tissue &
                                         grepl(paste(sel_group, collapse = "|"), contrast_short)
                                     ) %>%
                                     select(feature_id) %>%
                                     pull()))

  ol_all <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(assay =='prot-ol') %>%
    dplyr::filter(feature_id %in% (prot_ol_sig %>%
                                     filter(
                                       tissue==sel_tissue &
                                         grepl(paste(sel_group, collapse = "|"), contrast_short)
                                     ) %>%
                                     select(feature_id) %>%
                                     pull()))

  atac_seq_all <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(assay =='epigen-atac-seq') %>%
    dplyr::filter(feature_id %in% (atac_sig %>%
                                     filter(
                                       tissue==sel_tissue &
                                         grepl(paste(sel_group, collapse = "|"), contrast_short)
                                     ) %>%
                                     select(feature_id) %>%
                                     pull()))

  methyl_seq_all <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(assay =='epigen-methylcap-seq') %>%
    dplyr::filter(feature_id %in% (methyl_sig %>%
                                     filter(
                                       tissue==sel_tissue &
                                         grepl(paste(sel_group, collapse = "|"), contrast_short)
                                     ) %>%
                                     select(feature_id) %>%
                                     pull()))

  upset_data <- HUMAN_FEATURE_TO_GENE %>%
    dplyr::select(feature_id, ensembl_gene) %>%
    distinct(ensembl_gene, .keep_all = TRUE)

  upset_data$ATAC <- ifelse(upset_data$ensembl_gene %in% atac_seq_all$ensembl_gene, 1, 0)
  upset_data$Transcriptomics <- ifelse(upset_data$ensembl_gene %in% tr_all$ensembl_gene, 1, 0)
  upset_data$Proteomics <- ifelse(upset_data$ensembl_gene %in% pr_all$ensembl_gene, 1, 0)
  upset_data$OLINK <- ifelse(upset_data$ensembl_gene %in% ol_all$ensembl_gene, 1, 0)
  upset_data$Phosphoproteomics <- ifelse(upset_data$ensembl_gene %in% ph_all$ensembl_gene, 1, 0)
  upset_data$Methylation <- ifelse(upset_data$ensembl_gene %in% methyl_seq_all$ensembl_gene, 1, 0)

  upset_data <- upset_data %>%
    dplyr::filter(rowSums(across(where(is.numeric))) != 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-feature_id) %>%
    dplyr::filter(!(is.na(ensembl_gene))) %>%
    tibble::column_to_rownames('ensembl_gene') %>%
    dplyr::select(where(~ sum(.) != 0))

  # Return plot
  plot <- ComplexUpset::upset(
    upset_data,
    colnames(upset_data),
    encode_sets=FALSE,
    width_ratio=0.2,
    height_ratio=1,
    set_sizes=FALSE,
    base_annotations=list(
      'Number of \ndifferential genes'=intersection_size(
        counts=TRUE,
        fill=sel_tissue_col,
        color="black",
        text=list(size=4, color="black")
      )
    ),
    themes=upset_modify_themes(
      list(
        'intersections_matrix'=theme(axis.title.x=element_blank())
      )
    ),
    sort_intersections="descending",
    wrap=TRUE,
    stripes=col_seq
  ) +
    ggtitle(toupper(sel_tissue)) +
    theme(plot.title=element_text(hjust=0.5, size=16))

  return(plot)
}



upsetPlot_adipose <- make_upset_plot(
  "adipose",
  HUMAN_FEATURE_TO_GENE, prot_sig, ph_sig, transcript_sig, prot_ol_sig, atac_sig, methyl_sig,
  sel_group, adipose_col_seq, HUMAN_TISSUE_COLORS["adipose"]
)
upsetPlot_blood <- make_upset_plot(
  "blood",
  HUMAN_FEATURE_TO_GENE, prot_sig, ph_sig, transcript_sig, prot_ol_sig, atac_sig, methyl_sig,
  sel_group, blood_col_seq, HUMAN_TISSUE_COLORS["blood"]
)
upsetPlot_muscle <- make_upset_plot(
  "muscle",
  HUMAN_FEATURE_TO_GENE, prot_sig, ph_sig, transcript_sig, prot_ol_sig, atac_sig, methyl_sig,
  sel_group, muscle_col_seq, HUMAN_TISSUE_COLORS["muscle"]
)

combined_upsetPlot <- upsetPlot_adipose / upsetPlot_blood / upsetPlot_muscle
ggsave("combined_upset_plot_by_ome_3tissues.pdf", width = 6*1.25,height=9*1.25,
       units = "in")
