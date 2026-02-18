library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(MotrpacBicQC)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(PLIER)
library(gridExtra)
library(circlize)
library(ComplexHeatmap)
library(here)

# library(devtools)
# install_github("wgmao/PLIER")

setwd(file.path(here(), "figures", "landscape", "figure_5"))

# First, we must generate the molecular signature pathway matrix to be used as prior knowledge for PLIER
# This uses the MOLECULAR SIGNATURES variable in the MotrpacHumanPreSuspension package.

#Chris - updated Nov 11 2025. These are all the pathways related to transcriptomics.
#basically just excludes refmet, psp, ptmsigdb, cellmarker.
relevant_genesets = c("BIOCARTA", "KEGG_MEDICUS", "PID", "REACTOME", "WP", "GOBP", "GOCC", "GOMF", "MITOCARTA")
ourdatabases = MOLECULAR_SIGNATURES[names(MOLECULAR_SIGNATURES) %in% relevant_genesets]

combinedgenelist = unique(unlist(ourdatabases))
unlist_by_one = unlist(ourdatabases, recursive = FALSE)
names(unlist_by_one) = gsub("\\.", "_", names(unlist_by_one)) #by default unlist uses "." to combine, just replace w "_"
finalpathwaymat = sapply(unlist_by_one, function(x) as.integer(combinedgenelist %in% x))
rownames(finalpathwaymat) = combinedgenelist

# All the annotation colors for the plots we generate
precawg_ann_cols <- list("Tissue" = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS,
                         "Sex" = MotrpacHumanPreSuspensionAnalysis::HUMAN_SEX_COLORS[c("Male","Female")],
                         "Timepoint" = MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS[c("Pre_exercise","during_20_min","during_40_min","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr")],
                         "Modality" = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS)
names(precawg_ann_cols$Tissue) <- c("Blood","Muscle","Adipose")
names(precawg_ann_cols$Timepoint) <- c("Pre","D20M","D40M","P10M","P15-45M","P3.5/4H","P24H")
names(precawg_ann_cols$Modality) <- c("RE","EE","CON")

qc_norm = MotrpacHumanPreSuspensionData::load_qc(selected_omes = "transcript-rna-seq")

shared_transcripts = lapply(names(qc_norm), function (tissue){
  matrix_transcripts = rownames(qc_norm[[tissue]][["transcript-rna-seq"]][["qc_norm"]])
  matrix_transcripts
}) %>% reduce(intersect)

combined_matrixes = lapply(names(qc_norm), function(tissue){
  matrix_single = qc_norm[[tissue]][["transcript-rna-seq"]][["qc_norm"]] %>%
    filter(rownames(.) %in% shared_transcripts)
  t(scale(t(matrix_single))) %>% as.data.frame()
}) %>% reduce(cbind)

combined_meta = lapply(names(qc_norm), function(tissue){
  meta_single = qc_norm[[tissue]][["transcript-rna-seq"]][["sample_metadata"]] %>%
    select(Tissue,Sex,Timepoint,randomGroupCode,calculatedAge,BMI,pid)
  meta_single
}) %>% reduce(rbind)


# Data Frame connecting ensembl ids with gene symbols built off of Feature to Gene Matrix
#Chris - updated Nov 11 2025
gene_mapping = HUMAN_FEATURE_TO_GENE %>%
  dplyr::filter(assay == "transcript-rna-seq") %>%
  dplyr::mutate(
    feature_id = as.character(feature_id),
    gene_symbol = ifelse(is.na(gene_symbol), feature_id, as.character(gene_symbol)),
    ensembl_gene = as.character(ensembl_gene)
  ) %>%
  dplyr::filter(feature_id %in% rownames(combined_matrixes))
rownames(gene_mapping) = gene_mapping$feature_id

crosstissueRNAcomboPaths = finalpathwaymat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_symbol") %>%
  right_join(., gene_mapping %>% select(feature_id, gene_symbol), by = "gene_symbol") %>%
  tibble::column_to_rownames("feature_id") %>%
  dplyr::select(-"gene_symbol")

# Specifying the number of LVs we desire from PLIER. We typically use as input double the calculated number of PCs from the num.pc
# function; however, the combined RNAseq data is way too big to run num.pc efficiently, so we use 50 (100 LVs total) as an acceptable
# answer. This could be varied in the hopes of gaining more interesting LVs, but after testing, this was a sufficient number.

#note from chris: I get errors whenever I try to run num.pc or PLIER, because I get an error with `svd` not being supported.
#however, running svd or rsvd separately outside of the function seems to work. I don't know why

svd_calc = svd(combined_matrixes)
# crosstissuernacombovalue <- PLIER::num.pc(svd_calc, seed = 1) #212...too large to functionally annotate.
crosstissuernacombovalue <- 50 #we cap at 50 for this, which generates a max of 100 LVs

# Running PLIER with default parameters
crosstissuernacombofin.plierResult.all.7=PLIER(as.matrix(combined_matrixes),
                                               as.matrix(crosstissueRNAcomboPaths),
                                               k=2*crosstissuernacombovalue,
                                               trace=F,
                                               frac=0.7,
                                               scale=T,
                                               seed=1)

# Saving the PLIER output and load it in if needed. Can be computationally intensive
# saveRDS(crosstissuernacombofin.plierResult.all.7,file = "crosstissuernacombofin_plierresult.rds")
# crosstissuernacombofin.plierResult.all.7 = readRDS("crosstissuernacombofin_plierresult_102925.rds")

crosstissuernametadata <- crosstissuernametadata %>%
  arrange(Tissue, randomGroupCode, Timepoint, Sex) %>%
  mutate(
    Tissue = dplyr::recode(Tissue,
                           "Human Adipose Powder" = "Adipose",
                           "PaxGene RNA" = "Blood"
    ),
    Timepoint = dplyr::recode(Timepoint,
                              "pre_exercise" = "Pre",
                              "during_20_min" = "D20M",
                              "during_40_min" = "D40M",
                              "post_10_min" = "P10M",
                              "post_15_30_45_min" = "P15-45M",
                              "post_3.5_4_hr" = "P3.5/4H",
                              "post_24_hr" = "P24H"
    ),
    randomGroupCode = dplyr::recode(randomGroupCode,
                                    "ADUResist" = "RE",
                                    "ADUEndur" = "EE",
                                    "ADUControl" = "CON"
    )
  ) %>%
  rename(Modality = randomGroupCode) %>%
  mutate(pid = factor(pid))
combined_matrixes <- combined_matrixes[, rownames(crosstissuernametadata)]

# Generating heatmaps of z-scored expression across all samples for top genes associated with each LV
# Heatmaps of select LVs of interest are part of Supplemental Figure S5
for(i in 1:dim(crosstissuernacombofin.plierResult.all.7$Z)[2]){

  top10rows <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,i]),][1:10,])
  datatoplot <- combined_matrixes[top10rows,]

  png(filename=paste("crosstissuerna_PLIER_LV",toString(i),"Plot102925.png",sep=""),
      width = 8,height = 3,units = "in",res = 600)

  pheatmap(datatoplot,breaks = seq(-5, 5, length.out = 100),
           color = colorpanel(101,"blue", "white", "red"),
           show_colnames=F,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           annotation_col = crosstissuernametadata[,c("Sex","Timepoint","Modality","Tissue")],
           fontsize=14,
           labels_row = enstosym[gsub("\\..*","",top10rows),"gene_symbol"],
           annotation_colors = precawg_ann_cols,
           annotation_legend = FALSE)

  dev.off()
}

subjects_with_baseline <- crosstissuernametadata %>%
  filter(Timepoint == "Pre") %>%
  distinct(Tissue, pid)

adiposernasubjects <- subjects_with_baseline %>% filter(Tissue == "Adipose") %>% pull(pid)
musclernasubjects  <- subjects_with_baseline %>% filter(Tissue == "Muscle")  %>% pull(pid)
bloodrnasubjects   <- subjects_with_baseline %>% filter(Tissue == "Blood")   %>% pull(pid)

#####
# Making some line plots detailing individual LV response patterns used in the Figure 5A diagram
####

# The PLIER B matrix reflects the pattern of expression for each LV across all the samples in the input matrix.
crosstissuernabmat <- crosstissuernacombofin.plierResult.all.7$B
rownames(crosstissuernabmat) <- paste("LV",c(1:100),sep = "")

for (i in c(2, 6, 14, 16, 24, 31, 35, 89)) {
  png(
    file = paste0("crosstissue_rna_lv", i, "_lineplot_102925.png"),
    width = 6, height = 5, units = "in", res = 600
  )
  plot(crosstissuernabmat[i, colnames(combined_matrixes)], xlab = "Sample", ylab = "Relative Expression")
  lines(crosstissuernabmat[i, colnames(combined_matrixes)])
  dev.off()
}



# We're going to make a couple pathway enrichment plots for the three LVs highlighted in the Figure 5A diagram (14,16,24)

i = 14
lv14pathnames <- c("Cytosolic tRNA Aminoacylation",
                   "Ribosomal Biogenesis",
                   "snoRNA Binding",
                   "tRNA Processing",
                   "Preribosome")
png(file = "LV14_TopPathways_Barplot_072125.png",width = 8,height = 3,units = "in",res = 600)
ggbarplot(data.frame("Pathway" = rev(lv14pathnames),
                     "Uauc" = rev(crosstissuernacombofin.plierResult.all.7$Uauc[order(-crosstissuernacombofin.plierResult.all.7$Uauc[,14]),14][1:5])),x = "Pathway",y = "Uauc",fill = "Pathway",orientation = "horiz",palette = "Reds",legend = "none") + theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


i = 16
lv16pathnames <- c("Lymphoid Immunoregulation",
                   "Copy Number 17Q12 Variation Syndrome",
                   "T Cell Receptor Signaling",
                   "Alpha Beta T Cell Activation",
                   "B Cell Activation")
png(file = "LV16_TopPathways_Barplot_072125.png",width = 8,height = 3,units = "in",res = 600)
ggbarplot(data.frame("Pathway" = rev(lv16pathnames),
                     "Uauc" = rev(crosstissuernacombofin.plierResult.all.7$Uauc[order(-crosstissuernacombofin.plierResult.all.7$Uauc[,16]),16][1:5])),x = "Pathway",y = "Uauc",fill = "Pathway",orientation = "horiz",palette = "Reds",legend = "none") + theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


i = 24
lv24pathnames <- c("NGF Stimulated Transcription",
                   "ATF2 Pathway",
                   "Alpha Beta T Cell Activation",
                   "RNA Splicing Regulation",
                   "STK Signaling Regulation")
png(file = "LV24_TopPathways_Barplot_072125.png",width = 8,height = 3,units = "in",res = 600)
ggbarplot(data.frame("Pathway" = rev(lv24pathnames),
                     "Uauc" = rev(crosstissuernacombofin.plierResult.all.7$Uauc[order(-crosstissuernacombofin.plierResult.all.7$Uauc[,24]),24][1:5])),x = "Pathway",y = "Uauc",fill = "Pathway",orientation = "horiz",palette = "Reds",legend = "none") + theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


#####
# Now to generate the other important PLIER descriptor figures
####

# The PLIER Z matrix reflects the affinity of a feature to a given LV - important for identifying what features are markers for each LV
crosstissuernacombofin_Zmatscore <- scale(crosstissuernacombofin.plierResult.all.7$Z)

# Selecting features as aligned with a given LV if their z-scored affinity matrix value is greater than 3 (3 standard deviations above the mean LV affinity for that feature)
crosstissuernalv2markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,2]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[2],])
crosstissuernalv6markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,6]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[6],])
crosstissuernalv14markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,14]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[14],])
crosstissuernalv16markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,16]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[16],])
crosstissuernalv24markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,24]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[24],])
crosstissuernalv31markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,31]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[31],])
crosstissuernalv35markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,35]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[35],])
crosstissuernalv89markers <- rownames(crosstissuernacombofin.plierResult.all.7$Z[order(-crosstissuernacombofin.plierResult.all.7$Z[,89]),][1:colSums(crosstissuernacombofin_Zmatscore > 3)[89],])

write.csv(gsub("\\..*","",crosstissuernalv2markers),file = "newcrosstissuernalv2markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv6markers),file = "newcrosstissuernalv6markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv14markers),file = "newcrosstissuernalv14markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv16markers),file = "newcrosstissuernalv16markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv24markers),file = "newcrosstissuernalv24markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv31markers),file = "newcrosstissuernalv31markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv35markers),file = "newcrosstissuernalv35markers_102925.csv",row.names = F)
write.csv(gsub("\\..*","",crosstissuernalv89markers),file = "newcrosstissuernalv89markers_102925.csv",row.names = F)

crosstissuernalvmarkerlist <- list("LV2" = crosstissuernalv2markers,
                                   "LV6" = crosstissuernalv6markers,
                                   "LV14" = crosstissuernalv14markers,
                                   "LV16" = crosstissuernalv16markers,
                                   "LV24" = crosstissuernalv24markers,
                                   "LV31" = crosstissuernalv31markers,
                                   "LV35" = crosstissuernalv35markers,
                                   "LV89" = crosstissuernalv89markers)

#####
# Now, we want to plot how similar/different the responses are in the three tissues
# We will get a score for relationship/correlation between each tissue
# First we need to divide the B matrix into three matrices for the tissues
# Then we need to rename the matrix columns to reflect pid, timepoint
####

crosstissuernabmat_muscle <- crosstissuernabmat[,grep("muscle",colnames(crosstissuernabmat))]
crosstissuernabmat_adipose <- crosstissuernabmat[,grep("adipose",colnames(crosstissuernabmat))]
crosstissuernabmat_blood <- crosstissuernabmat[,grep("blood",colnames(crosstissuernabmat))]

colnames(crosstissuernabmat_muscle) <- gsub("muscle_","",colnames(crosstissuernabmat_muscle))
musclernasamplemetadata$pid_timepoint <- paste(musclernasamplemetadata$pid,musclernasamplemetadata$Timepoint,sep = "_")
colnames(crosstissuernabmat_muscle) <- musclernasamplemetadata[colnames(crosstissuernabmat_muscle),"pid_timepoint"]

colnames(crosstissuernabmat_adipose) <- gsub("adipose_","",colnames(crosstissuernabmat_adipose))
adiposernasamplemetadata$pid_timepoint <- paste(adiposernasamplemetadata$pid,adiposernasamplemetadata$Timepoint,sep = "_")
colnames(crosstissuernabmat_adipose) <- adiposernasamplemetadata[colnames(crosstissuernabmat_adipose),"pid_timepoint"]

colnames(crosstissuernabmat_blood) <- gsub("blood_","",colnames(crosstissuernabmat_blood))
bloodrnasamplemetadata$pid_timepoint <- paste(bloodrnasamplemetadata$pid,bloodrnasamplemetadata$Timepoint,sep = "_")
colnames(crosstissuernabmat_blood) <- bloodrnasamplemetadata[colnames(crosstissuernabmat_blood),"pid_timepoint"]

# Calculating the correlation of LV pattern across each tissue pair
crosstissuernabmat_tissuecordf <- data.frame(row.names = rownames(crosstissuernabmat),
                                             "Adipose_v_Blood" = rep(0,length(rownames(crosstissuernabmat))),
                                             "Adipose_v_Muscle" = rep(0,length(rownames(crosstissuernabmat))),
                                             "Blood_v_Muscle" = rep(0,length(rownames(crosstissuernabmat))))
for(i in 1:dim(crosstissuernabmat_tissuecordf)[1]){
  crosstissuernabmat_tissuecordf[i,"Adipose_v_Blood"] <- cor(crosstissuernabmat_adipose[i,intersect(colnames(crosstissuernabmat_adipose),colnames(crosstissuernabmat_blood))],
                                                             crosstissuernabmat_blood[i,intersect(colnames(crosstissuernabmat_adipose),colnames(crosstissuernabmat_blood))])
  crosstissuernabmat_tissuecordf[i,"Adipose_v_Muscle"] <- cor(crosstissuernabmat_adipose[i,intersect(colnames(crosstissuernabmat_adipose),colnames(crosstissuernabmat_muscle))],
                                                              crosstissuernabmat_muscle[i,intersect(colnames(crosstissuernabmat_adipose),colnames(crosstissuernabmat_muscle))])
  crosstissuernabmat_tissuecordf[i,"Blood_v_Muscle"] <- cor(crosstissuernabmat_blood[i,intersect(colnames(crosstissuernabmat_blood),colnames(crosstissuernabmat_muscle))],
                                                            crosstissuernabmat_muscle[i,intersect(colnames(crosstissuernabmat_blood),colnames(crosstissuernabmat_muscle))])
}
crosstissuernabmat_tissuecordf$LV <- rownames(crosstissuernabmat_tissuecordf)
crosstissuernabmat_tissuecordf$SigLV <- 1 + 1*(crosstissuernabmat_tissuecordf$LV %in% names(crosstissuernalvmarkerlist))

crosstissuernabmat_tissuecordf_forggplot <- data.frame("Correlation" = c(crosstissuernabmat_tissuecordf$Adipose_v_Blood,
                                                                         crosstissuernabmat_tissuecordf$Adipose_v_Muscle,
                                                                         crosstissuernabmat_tissuecordf$Blood_v_Muscle),
                                                       "Comparison" = c(rep("Adipose vs Blood",length(crosstissuernabmat_tissuecordf$Adipose_v_Blood)),
                                                                        rep("Adipose vs Muscle",length(crosstissuernabmat_tissuecordf$Adipose_v_Muscle)),
                                                                        rep("Blood vs Muscle",length(crosstissuernabmat_tissuecordf$Blood_v_Muscle))),
                                                       "LV" = rep(rownames(crosstissuernabmat_tissuecordf),3))
crosstissuernabmat_tissuecordf_forggplot$LV <- factor(crosstissuernabmat_tissuecordf_forggplot$LV,levels = crosstissuernabmat_tissuecordf_forggplot$LV[order(crosstissuernabmat_tissuecordf$Blood_v_Muscle)])

# Data frame of the count of markers associated with each LV
crosstissuernalvmembershipbarplotdf <- data.frame("LV" = c(rep("LV2",length(crosstissuernalvmarkerlist$LV2)),
                                                           rep("LV6",length(crosstissuernalvmarkerlist$LV6)),
                                                           rep("LV14",length(crosstissuernalvmarkerlist$LV14)),
                                                           rep("LV16",length(crosstissuernalvmarkerlist$LV16)),
                                                           rep("LV24",length(crosstissuernalvmarkerlist$LV24)),
                                                           rep("LV31",length(crosstissuernalvmarkerlist$LV31)),
                                                           rep("LV35",length(crosstissuernalvmarkerlist$LV35)),
                                                           rep("LV89",length(crosstissuernalvmarkerlist$LV89))))

crosstissuernalvmembershipbarplotdf$LV <- factor(crosstissuernalvmembershipbarplotdf$LV,levels = crosstissuernabmat_tissuecordf_forggplot$LV[order(crosstissuernabmat_tissuecordf$Blood_v_Muscle)])

# Figure 5C - LV tissue correlation
png(file = "Figure5C_crosstissuerna_lvtissuecorrelation_102925.png",width = 4,height = 4,units = "in",res = 600)
ggplot(crosstissuernabmat_tissuecordf_forggplot[crosstissuernabmat_tissuecordf_forggplot$LV %in% names(crosstissuernalvmarkerlist),]) + geom_hline(yintercept = 0,size = 1,linetype = "dashed") + geom_point(aes(x=LV,y=Correlation,color=Comparison),size = 3) + theme_classic() + scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7")) + theme(axis.title.x=element_blank(),legend.position = "top",legend,text = element_text(size=10),legend.key.spacing = unit(0,"cm")) + guides(col = guide_legend(ncol = 1))
dev.off()

# Figure 5D - membership barplot
png(file = "Figure5D_crosstissuerna_lvmembershipbarplot_102925.png",width = 4,height = 3,units = "in",res = 600)
ggplot(crosstissuernalvmembershipbarplotdf,aes(LV)) + geom_bar() + theme_classic() + ylab("Membership") + theme(axis.title.x=element_blank(),text = element_text(size=10))
dev.off()

#####
# Now, we want to make a plot of the significant (up or down) responses by each LV
# So we will need to use compare_means  to generate a p-val for each LV for each comparison
# - then create a heatmap with a column for each Tissue/Time-Point/Endurance+Resistance Groups
####

crosstissuerna_lvsigtestdf <- matrix(0L,nrow = 100,ncol = 22)
rownames(crosstissuerna_lvsigtestdf) <- paste("LV",c(1:100),sep = "")
colnames(crosstissuerna_lvsigtestdf) <- c("Adipose_45min_Endurance","Adipose_4hr_Endurance","Adipose_24hr_Endurance",
                                          "Adipose_45min_Resistance","Adipose_4hr_Resistance","Adipose_24hr_Resistance",
                                          "Blood_20minDuring_Endurance","Blood_40minDuring_Endurance","Blood_10min_Endurance",
                                          "Blood_30min_Endurance","Blood_4hr_Endurance","Blood_24hr_Endurance",
                                          "Blood_10min_Resistance","Blood_30min_Resistance","Blood_4hr_Resistance","Blood_24hr_Resistance",
                                          "Muscle_15min_Endurance","Muscle_3_5hr_Endurance","Muscle_24hr_Endurance",
                                          "Muscle_15min_Resistance","Muscle_3_5hr_Resistance","Muscle_24hr_Resistance")

for(i in 1:dim(crosstissuerna_lvsigtestdf)[1]){

  metamerge <- trimcrosstissuernametadata
  metamerge$B <- crosstissuernacombofin.plierResult.all.7$B[i,rownames(trimcrosstissuernametadata)]
  metamerge$pid <- factor(metamerge$pid,levels = unique(metamerge$pid))
  metamerge$Diff <- 0

  metamergeadipose <- metamerge[metamerge$Tissue %in% "Adipose",]
  metamergemuscle <- metamerge[metamerge$Tissue %in% "Muscle",]
  metamergeblood <- metamerge[metamerge$Tissue %in% "Blood",]

  metamergeadipose <- metamergeadipose[metamergeadipose$pid %in% adiposernasubjects,]
  metamergemuscle <- metamergemuscle[metamergemuscle$pid %in% musclernasubjects,]
  metamergeblood <- metamergeblood[metamergeblood$pid %in% bloodrnasubjects,]

  for(j in 1:dim(metamergeadipose)[1]){
    ourpid <- metamergeadipose[j,"pid"]
    ourtissue <- metamergeadipose[j,"Tissue"]
    ourbaseline <- metamergeadipose[metamergeadipose$pid %in% ourpid & metamergeadipose$Timepoint %in% "Pre" & metamergeadipose$Tissue %in% ourtissue,"B"]
    metamergeadipose[j,"Diff"] <- metamergeadipose[j,"B"] - ourbaseline
  }

  for(j in 1:dim(metamergemuscle)[1]){
    ourpid <- metamergemuscle[j,"pid"]
    ourtissue <- metamergemuscle[j,"Tissue"]
    ourbaseline <- metamergemuscle[metamergemuscle$pid %in% ourpid & metamergemuscle$Timepoint %in% "Pre" & metamergemuscle$Tissue %in% ourtissue,"B"]
    metamergemuscle[j,"Diff"] <- metamergemuscle[j,"B"] - ourbaseline
  }

  for(j in 1:dim(metamergeblood)[1]){
    ourpid <- metamergeblood[j,"pid"]
    ourtissue <- metamergeblood[j,"Tissue"]
    ourbaseline <- metamergeblood[metamergeblood$pid %in% ourpid & metamergeblood$Timepoint %in% "Pre" & metamergeblood$Tissue %in% ourtissue,"B"]
    metamergeblood[j,"Diff"] <- metamergeblood[j,"B"] - ourbaseline
  }

  adiposetestout <- compare_means(Diff ~ Modality,data = metamergeadipose[!(metamergeadipose$Timepoint %in% c("Pre")),],ref.group = "CON",group.by = "Timepoint")
  muscletestout <- compare_means(Diff ~ Modality,data = metamergemuscle[!(metamergemuscle$Timepoint %in% c("Pre")),],ref.group = "CON",group.by = "Timepoint")
  bloodtestout <- compare_means(Diff ~ Modality,data = metamergeblood[!(metamergeblood$Timepoint %in% c("Pre")),],ref.group = "CON",group.by = "Timepoint")

  crosstissuerna_lvsigtestdf[i,"Adipose_45min_Endurance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "EE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[1,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Adipose_45min_Resistance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "RE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[2,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Adipose_4hr_Endurance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "EE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[3,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Adipose_4hr_Resistance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "RE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[4,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Adipose_24hr_Endurance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "EE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[5,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Adipose_24hr_Resistance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "RE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[6,"p.adj"][[1]]))

  crosstissuerna_lvsigtestdf[i,"Muscle_15min_Endurance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "EE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[1,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Muscle_15min_Resistance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "RE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[2,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Muscle_3_5hr_Endurance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "EE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[3,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Muscle_3_5hr_Resistance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "RE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[4,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Muscle_24hr_Endurance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "EE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[5,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Muscle_24hr_Resistance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "RE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[6,"p.adj"][[1]]))

  crosstissuerna_lvsigtestdf[i,"Blood_20minDuring_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[1,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_40minDuring_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[2,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_10min_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[3,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_10min_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[4,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_30min_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[5,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_30min_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[6,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_4hr_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[7,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_4hr_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[8,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_24hr_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[9,"p.adj"][[1]]))
  crosstissuerna_lvsigtestdf[i,"Blood_24hr_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[10,"p.adj"][[1]]))

}


crosstissuerna_lvsigtestdfmeta <- data.frame(row.names = colnames(crosstissuerna_lvsigtestdf),
                                             "Tissue" = c(rep("Adipose",6),
                                                          rep("Blood",10),
                                                          rep("Muscle",6)),
                                             "Timepoint" = c("P15-45M","P3.5/4H","P24H",
                                                             "P15-45M","P3.5/4H","P24H",
                                                             "D20M","D40M","P10M",
                                                             "P15-45M","P3.5/4H","P24H",
                                                             "P10M","P15-45M","P3.5/4H","P24H",
                                                             "P15-45M","P3.5/4H","P24H",
                                                             "P15-45M","P3.5/4H","P24H"),
                                             "Modality" = c("EE","EE","EE",
                                                            "RE","RE","RE",
                                                            "EE","EE","EE",
                                                            "EE","EE","EE",
                                                            "RE","RE","RE","RE",
                                                            "EE","EE","EE",
                                                            "RE","RE","RE"))

crosstissuerna_lvsigtestdfmeta$Tissue <- factor(crosstissuerna_lvsigtestdfmeta$Tissue,levels = c("Adipose","Blood","Muscle"))
crosstissuerna_lvsigtestdfmeta$Timepoint <- factor(crosstissuerna_lvsigtestdfmeta$Timepoint,levels = c("D20M","D40M","P10M","P15-45M","P3.5/4H","P24H"))
crosstissuerna_lvsigtestdfmeta$Modality <- factor(crosstissuerna_lvsigtestdfmeta$Modality,levels = c("EE","RE"))

crosstissuerna_lvsigtestdfmeta <- crosstissuerna_lvsigtestdfmeta[order(crosstissuerna_lvsigtestdfmeta$Tissue,crosstissuerna_lvsigtestdfmeta$Modality,crosstissuerna_lvsigtestdfmeta$Timepoint),]
crosstissuerna_lvsigtestdf <- crosstissuerna_lvsigtestdf[,rownames(crosstissuerna_lvsigtestdfmeta)]

#####
# Similar to above, but now we generate a finalized form with bubbleheatmap
####

crosstissuerna_lvsigtest_bhsize_df <- matrix(0L,nrow = 100,ncol = 22)
rownames(crosstissuerna_lvsigtest_bhsize_df) <- paste("LV",c(1:100),sep = "")
colnames(crosstissuerna_lvsigtest_bhsize_df) <- c("Adipose_45min_Endurance","Adipose_4hr_Endurance","Adipose_24hr_Endurance",
                                                  "Adipose_45min_Resistance","Adipose_4hr_Resistance","Adipose_24hr_Resistance",
                                                  "Blood_20minDuring_Endurance","Blood_40minDuring_Endurance","Blood_10min_Endurance",
                                                  "Blood_30min_Endurance","Blood_4hr_Endurance","Blood_24hr_Endurance",
                                                  "Blood_10min_Resistance","Blood_30min_Resistance","Blood_4hr_Resistance","Blood_24hr_Resistance",
                                                  "Muscle_15min_Endurance","Muscle_3_5hr_Endurance","Muscle_24hr_Endurance",
                                                  "Muscle_15min_Resistance","Muscle_3_5hr_Resistance","Muscle_24hr_Resistance")

crosstissuerna_lvsigtest_bhcolor_df <- matrix(0L,nrow = 100,ncol = 22)
rownames(crosstissuerna_lvsigtest_bhcolor_df) <- paste("LV",c(1:100),sep = "")
colnames(crosstissuerna_lvsigtest_bhcolor_df) <- c("Adipose_45min_Endurance","Adipose_4hr_Endurance","Adipose_24hr_Endurance",
                                                   "Adipose_45min_Resistance","Adipose_4hr_Resistance","Adipose_24hr_Resistance",
                                                   "Blood_20minDuring_Endurance","Blood_40minDuring_Endurance","Blood_10min_Endurance",
                                                   "Blood_30min_Endurance","Blood_4hr_Endurance","Blood_24hr_Endurance",
                                                   "Blood_10min_Resistance","Blood_30min_Resistance","Blood_4hr_Resistance","Blood_24hr_Resistance",
                                                   "Muscle_15min_Endurance","Muscle_3_5hr_Endurance","Muscle_24hr_Endurance",
                                                   "Muscle_15min_Resistance","Muscle_3_5hr_Resistance","Muscle_24hr_Resistance")

for(i in 1:dim(crosstissuerna_lvsigtestdf)[1]){

  metamerge <- trimcrosstissuernametadata
  metamerge$B <- crosstissuernacombofin.plierResult.all.7$B[i,rownames(trimcrosstissuernametadata)]
  metamerge$pid <- factor(metamerge$pid,levels = unique(metamerge$pid))
  metamerge$Diff <- 0

  metamergeadipose <- metamerge[metamerge$Tissue %in% "Adipose",]
  metamergemuscle <- metamerge[metamerge$Tissue %in% "Muscle",]
  metamergeblood <- metamerge[metamerge$Tissue %in% "Blood",]

  metamergeadipose <- metamergeadipose[metamergeadipose$pid %in% adiposernasubjects,]
  metamergemuscle <- metamergemuscle[metamergemuscle$pid %in% musclernasubjects,]
  metamergeblood <- metamergeblood[metamergeblood$pid %in% bloodrnasubjects,]

  for(j in 1:dim(metamergeadipose)[1]){
    ourpid <- metamergeadipose[j,"pid"]
    ourtissue <- metamergeadipose[j,"Tissue"]
    ourbaseline <- metamergeadipose[metamergeadipose$pid %in% ourpid & metamergeadipose$Timepoint %in% "Pre" & metamergeadipose$Tissue %in% ourtissue,"B"]
    metamergeadipose[j,"Diff"] <- metamergeadipose[j,"B"] - ourbaseline
  }

  for(j in 1:dim(metamergemuscle)[1]){
    ourpid <- metamergemuscle[j,"pid"]
    ourtissue <- metamergemuscle[j,"Tissue"]
    ourbaseline <- metamergemuscle[metamergemuscle$pid %in% ourpid & metamergemuscle$Timepoint %in% "Pre" & metamergemuscle$Tissue %in% ourtissue,"B"]
    metamergemuscle[j,"Diff"] <- metamergemuscle[j,"B"] - ourbaseline
  }

  for(j in 1:dim(metamergeblood)[1]){
    ourpid <- metamergeblood[j,"pid"]
    ourtissue <- metamergeblood[j,"Tissue"]
    ourbaseline <- metamergeblood[metamergeblood$pid %in% ourpid & metamergeblood$Timepoint %in% "Pre" & metamergeblood$Tissue %in% ourtissue,"B"]
    metamergeblood[j,"Diff"] <- metamergeblood[j,"B"] - ourbaseline
  }

  adiposetestout <- compare_means(Diff ~ Modality,data = metamergeadipose[!(metamergeadipose$Timepoint %in% c("Pre")),],ref.group = "CON",group.by = "Timepoint")
  muscletestout <- compare_means(Diff ~ Modality,data = metamergemuscle[!(metamergemuscle$Timepoint %in% c("Pre")),],ref.group = "CON",group.by = "Timepoint")
  bloodtestout <- compare_means(Diff ~ Modality,data = metamergeblood[!(metamergeblood$Timepoint %in% c("Pre")),],ref.group = "CON",group.by = "Timepoint")

  crosstissuerna_lvsigtest_bhsize_df[i,"Adipose_45min_Endurance"] <- -log10(adiposetestout[1,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Adipose_45min_Resistance"] <- -log10(adiposetestout[2,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Adipose_4hr_Endurance"] <- -log10(adiposetestout[3,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Adipose_4hr_Resistance"] <- -log10(adiposetestout[4,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Adipose_24hr_Endurance"] <- -log10(adiposetestout[5,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Adipose_24hr_Resistance"] <- -log10(adiposetestout[6,"p.adj"][[1]])

  crosstissuerna_lvsigtest_bhsize_df[i,"Muscle_15min_Endurance"] <- -log10(muscletestout[1,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Muscle_15min_Resistance"] <- -log10(muscletestout[2,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Muscle_3_5hr_Endurance"] <- -log10(muscletestout[3,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Muscle_3_5hr_Resistance"] <- -log10(muscletestout[4,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Muscle_24hr_Endurance"] <- -log10(muscletestout[5,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Muscle_24hr_Resistance"] <- -log10(muscletestout[6,"p.adj"][[1]])

  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_20minDuring_Endurance"] <- -log10(bloodtestout[1,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_40minDuring_Endurance"] <- -log10(bloodtestout[2,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_10min_Endurance"] <- -log10(bloodtestout[3,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_10min_Resistance"] <- -log10(bloodtestout[4,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_30min_Endurance"] <- -log10(bloodtestout[5,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_30min_Resistance"] <- -log10(bloodtestout[6,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_4hr_Endurance"] <- -log10(bloodtestout[7,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_4hr_Resistance"] <- -log10(bloodtestout[8,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_24hr_Endurance"] <- -log10(bloodtestout[9,"p.adj"][[1]])
  crosstissuerna_lvsigtest_bhsize_df[i,"Blood_24hr_Resistance"] <- -log10(bloodtestout[10,"p.adj"][[1]])

  crosstissuerna_lvsigtest_bhcolor_df[i,"Adipose_45min_Endurance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "EE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Adipose_45min_Resistance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "RE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Adipose_4hr_Endurance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "EE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Adipose_4hr_Resistance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "RE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Adipose_24hr_Endurance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "EE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Adipose_24hr_Resistance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "RE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"])

  crosstissuerna_lvsigtest_bhcolor_df[i,"Muscle_15min_Endurance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "EE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Muscle_15min_Resistance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "RE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Muscle_3_5hr_Endurance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "EE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Muscle_3_5hr_Resistance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "RE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Muscle_24hr_Endurance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "EE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Muscle_24hr_Resistance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "RE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"])

  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_20minDuring_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_40minDuring_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_10min_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_10min_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_30min_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_30min_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_4hr_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_4hr_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_24hr_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuerna_lvsigtest_bhcolor_df[i,"Blood_24hr_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"])
}

#testout <- bubbleHeatmap(colorMat = crosstissuerna_lvsigtest_bhcolor_df[names(crosstissuernalvmarkerlist),],sizeMat = crosstissuerna_lvsigtest_bhsize_df[names(crosstissuernalvmarkerlist),],legendTitles = c("Significance","log2FC"),colorSeq = c("blue","white","red"),colorLim = c(-0.5,0.5),sizeLim = c(0,5))

testout <- bubbleHeatmap(colorMat = crosstissuerna_lvsigtest_bhcolor_df[levels(crosstissuernalvmembershipbarplotdf$LV)[levels(crosstissuernalvmembershipbarplotdf$LV) %in% names(crosstissuernalvmarkerlist)],],
                         sizeMat = crosstissuerna_lvsigtest_bhsize_df[levels(crosstissuernalvmembershipbarplotdf$LV)[levels(crosstissuernalvmembershipbarplotdf$LV) %in% names(crosstissuernalvmarkerlist)],],
                         legendTitles = c("Significance","log2FC"),colorSeq = c("blue","white","red"),colorLim = c(-0.5,0.5),sizeLim = c(0,5))


#--------------------------------------------------------------------------------------
#Chris: This code is to replace the above bubbleheatmap code to keep consistent
#with all of the other plots and to avoid having to combine different figures
#in illustrator. The functions are largely similar to those from 6B/S6.
# In the future, i should trim these to avoid adding too much code and make a pliers_function.R file i source.
# Figure 5B - Nov 17th, 2025.

keep_rows = levels(crosstissuernalvmembershipbarplotdf$LV)[levels(crosstissuernalvmembershipbarplotdf$LV) %in% names(crosstissuernalvmarkerlist)]
colorMat = crosstissuerna_lvsigtest_bhcolor_df[keep_rows,]
sizeMat  = crosstissuerna_lvsigtest_bhsize_df[keep_rows,]

#here greg names things differently compared to figure 6, so i have to change the way i parse text.
#the timepoint grep part is changed.
annotate_matrix = function(mat) {
  mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("size_info") %>%
    tidyr::pivot_longer(
      cols = -size_info,
      names_to = "tissue_mod_timepoint",
      values_to = "value"
    )  %>%
    dplyr::mutate(
      Tissue = dplyr::case_when(
        grepl("Adipose", tissue_mod_timepoint, ignore.case = TRUE) ~ "adipose",
        grepl("Blood", tissue_mod_timepoint, ignore.case = TRUE) ~ "blood",
        grepl("Muscle", tissue_mod_timepoint, ignore.case = TRUE) ~ "muscle",
        TRUE ~ NA_character_
      ),
      Modality = dplyr::case_when(
        grepl("Endur", tissue_mod_timepoint, ignore.case = TRUE) ~ "EE",
        grepl("Resist", tissue_mod_timepoint, ignore.case = TRUE) ~ "RE",
        TRUE ~ NA_character_
      ),
      Timepoint = dplyr::case_when(
        grepl("20", tissue_mod_timepoint, ignore.case = TRUE) ~ "D20M",
        grepl("40", tissue_mod_timepoint, ignore.case = TRUE) ~ "D40M",
        grepl("10", tissue_mod_timepoint, ignore.case = TRUE) ~ "P10M",
        grepl("15|30|45", tissue_mod_timepoint, ignore.case = TRUE) ~ "P15-45M",
        grepl("_3_5hr|_4hr_", tissue_mod_timepoint, ignore.case = TRUE) ~ "P3.5/4H",
        grepl("24", tissue_mod_timepoint, ignore.case = TRUE) ~ "P24H",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(Tissue, Modality, Timepoint, tissue_mod_timepoint, value)
}

anno_df <- annotate_matrix(sizeMat) %>%
  dplyr::select(Tissue, Modality, Timepoint) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D20M",
                                                  "D40M",
                                                  "P10M",
                                                  "P15-45M",
                                                  "P3.5/4H",
                                                  "P24H"))) %>%
  distinct()

.contrast_colors <- function() {
  structure(
    c("#fde725",
      "#bad071",
      "#d1bbd7",
      "#ae76a3",
      "#882e72",
      "#61194f"),
    names = c("D20M",
              "D40M",
              "P10M",
              "P15-45M",
              "P3.5/4H",
              "P24H")
  )
}

anno_col <- list(
  "Tissue" = MotrpacHumanPreSuspension::HUMAN_TISSUE_COLORS,
  "Modality" = setNames(c("#d95f02", "#1b9e77"),
                        c("EE", "RE")),
  "Timepoint" = .contrast_colors()[levels(anno_df$Timepoint)]
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = anno_df,
  col = anno_col,
  which = "column",
  border = TRUE,
  gap = unit(2, "pt"),
  annotation_name_gp = gpar(fontsize = 0.9 * 14),
  annotation_legend_param = list(
    border = TRUE,
    title_gp = gpar(fontsize = 0.9 * unit(14, "pt"),
                    fontface = "bold"),
    labels_gp = gpar(fontsize = 0.9 * unit(14, "pt"))
  ),
  annotation_name_side = "left"
)

col_fun = circlize::colorRamp2(
  c(-1, 0, 1),
  c("blue", "white", "red")
)

color_lg = ComplexHeatmap::Legend(
  title = "log2FC",
  col_fun = col_fun,
  at = c(-1, 0, 1),
  labels = c("-1", "0", "1")
)

size_lg = ComplexHeatmap::Legend(
  title = "-log10(adj_p)",
  type = "points",
  pch = 16, #circle code
  size = grid::unit(c(0, 1, 2, 3, 4, 5), "mm"), #this is diameter
  labels = c("0", "1", "2", "3", "4", "5+")
)
#replace bubble heatmap with Complexheatmap to keep consistent w other figures
cell_size = unit(6, "mm")
make_ht = function(colorMat_sub, sizeMat_sub) {
  ComplexHeatmap::Heatmap(
    matrix = matrix(NA_real_,
                    nrow = nrow(colorMat_sub),
                    ncol = ncol(colorMat_sub),
                    dimnames = dimnames(colorMat_sub)),

    width  = cell_size * ncol(colorMat_sub),
    height = cell_size * nrow(colorMat_sub),

    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(fill = NA, col = "black"),

    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_heatmap_legend = FALSE,

    top_annotation = top_annotation,

    cell_fun = function(j, i, x, y, width, height, fill) {

      fc = colorMat_sub[i, j]
      pv = sizeMat_sub[i, j]
      if (pv > 5) pv = 5
      pv = unit(pv, "mm")

      if (!is.na(fc) && !is.na(pv)) {
        bubble_color = col_fun(fc)
        bubble_radius = pv / 2

        grid::grid.circle(
          x = x, y = y,
          r = bubble_radius,
          gp = grid::gpar(fill = bubble_color, col = NA)
        )
      }
    }
  )
}

ht = make_ht(colorMat, sizeMat)

pdf(file = "5B_muscle_TF_protph_bubbleheatmap_chris.pdf", height = 8, width = 12)
ComplexHeatmap::draw(
  ht,
  heatmap_legend_list = list(color_lg, size_lg),
  merge_legend = FALSE
)
dev.off()
#-----------------------------


#####
# Generating a pathway enrichment heatmap for the LVs of interest in the dataset
####

crosstissuernatoppathways <- Reduce(union,list(rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,2]),])[1:1],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,6]),])[1:5],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,14]),])[1:5],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,16]),])[1:5],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,24]),])[1:5],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,31]),])[1:5],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,35]),])[1:5],
                                               rownames(crosstissuernacombofin.plierResult.all.7$U[order(-crosstissuernacombofin.plierResult.all.7$U[,89]),])[1:5]))

ourindexCol = c(2,6,14,16,24,31,35,89)

crosstissuernatoppathwaysdisplay <- c("DNA Binding Transcription Factor Pathway",
                                      "RNA Splicing Regulation",
                                      "Capped Intron Containing Pre mRNA Processing",
                                      "Nuclear Transcribed mRNA Catabolic Process",
                                      "Transcription Regulator Activity",
                                      "Sequence Specific DNA Binding",
                                      "Preribosome",
                                      "Ribosomal Large Subunit Biogenesis",
                                      "Translation",
                                      "Protein Folding Chaperone Complex",
                                      "Cytosolic tRNA Aminoacylation",
                                      "Lymphoid Cell Immunoregulatory Interactions",
                                      "Alpha Beta T Cell Activation",
                                      "G Protein Coupled Receptor Activity",
                                      "T Cell Receptor Signaling Pathway",
                                      "Copy Number 17Q12 Variation Syndrome",
                                      "NGF Stimulated Transcription",
                                      "Glucocorticoid Receptor Pathway",
                                      "ATF2 Pathway",
                                      "ncRNA Transcription Regulation",
                                      "Fibroblast Growth Factor Response",
                                      "Phagophore Assembly Site",
                                      "Chromosome Segregation Regulation",
                                      "WNT Signaling Pathway Regulation",
                                      "Attenuation Phase",
                                      "HSF1 Activation",
                                      "Actin Filament Binding",
                                      "HSP90 Chaperone Cycle for Steroid Hormone Receptors",
                                      "MAPK6 MAPK4 Signaling")

#Figure 5E
png(file = "Figure5E_Cross Tissue RNA PLIER Top 5 Pathways for Top LVs_102925.png",width = 8,height = 8,units = "in",res = 600)
pheatmap(crosstissuernacombofin.plierResult.all.7$Uauc[crosstissuernatoppathways,as.numeric(gsub("LV","",levels(crosstissuernalvmembershipbarplotdf$LV)[levels(crosstissuernalvmembershipbarplotdf$LV) %in% names(crosstissuernalvmarkerlist)]))],breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,angle_col = 0,labels_row = crosstissuernatoppathwaysdisplay)
dev.off()



#####
# Making the RNAseq portion of the Supplemental Table S5
####

ourindexCol <- as.numeric(gsub("LV","",names(crosstissuernalvmarkerlist)))

toplvZmat <- crosstissuernacombofin.plierResult.all.7$Z[Reduce(union,crosstissuernalvmarkerlist),ourindexCol]
colnames(toplvZmat) <- paste("LV",ourindexCol,sep = "")
toplvUmat <- crosstissuernacombofin.plierResult.all.7$Uauc[,ourindexCol]
toplvUmat <- toplvUmat[apply(toplvUmat,1,max) > 0,]

RNAPLIER_supplementaltable <- rbind(toplvUmat,toplvZmat)
write.csv(RNAPLIER_supplementaltable,file = "Supplemental Table S5_RNAseq_102925.csv",row.names = T)

# In case you want to save progress
#save.image("precovid_rnaplieranalysis_102925.RData")
