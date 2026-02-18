library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacBicQC)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(Mfuzz)
library(gprofiler2)
library(ggpubr)
library(PLIER)
library(gridExtra)
library(biomaRt)
library(bubbleHeatmap)
library(ggsankey)
library(TMSig)
library(ggrepel)
library(geomtextpath)
library(here)

#note: this script is currently undergoing refactoring.

setwd(file.path(here(), "figures", "landscape", "figure_6"))


# First thing is to identify the DEGs for each tissue comparison for homer enrichment analysis
full_degs = MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(selected_omes = "transcript-rna-seq",
                                                                      single_matrix = TRUE,
                                                                      combine_with_featgene = TRUE)
group_vars = c("tissue", "randomGroupCode", "Timepoint")
unique_groups = full_degs %>%
  dplyr::distinct(dplyr::across(all_of(group_vars)))

for(i in seq_len(nrow(unique_groups))) {
  this_tissue = unique_groups$tissue[i]
  this_code = unique_groups$randomGroupCode[i]
  this_tp = unique_groups$Timepoint[i]

  subset_df = full_degs %>%
    dplyr::filter(p_value < 0.05,
                  tissue == this_tissue,
                  randomGroupCode == this_code,
                  Timepoint == this_tp)

  outfile = paste0("homer_", this_tissue, "_",this_code, "_",this_tp,".txt")
  # We write the files that are used as input for HOMER. HOMER is run outside of R using the findMotifs.pl command with default parameters
  write.table(subset_df, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
}


# These DEG lists are used as inputs in homer with the command findMotifs.pl bloodrna_immpost_resist_sig.txt human bloodrna_immpost_resist_sig_output/
# Then we load the enrichment results of all TFs for each DEG set here for further analysis
saved_homer_files_path = file.path(here(), "data", "tmp", "Precovid_DEG_HOMER_KnownTF_Results")
#ok pause here.
musclerna_early_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/muscle_early_resist_sig_knownResults.txt",header = T,sep = "\t")
musclerna_early_resist_sig_tfout <- musclerna_early_resist_sig_tfout[!duplicated(musclerna_early_resist_sig_tfout$Motif.Name),]
rownames(musclerna_early_resist_sig_tfout) <- musclerna_early_resist_sig_tfout$Motif.Name

musclerna_mid_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/muscle_mid_resist_sig_knownResults.txt",header = T,sep = "\t")
musclerna_mid_resist_sig_tfout <- musclerna_mid_resist_sig_tfout[!duplicated(musclerna_mid_resist_sig_tfout$Motif.Name),]
rownames(musclerna_mid_resist_sig_tfout) <- musclerna_mid_resist_sig_tfout$Motif.Name

musclerna_late_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/muscle_late_resist_sig_knownResults.txt",header = T,sep = "\t")
musclerna_late_resist_sig_tfout <- musclerna_late_resist_sig_tfout[!duplicated(musclerna_late_resist_sig_tfout$Motif.Name),]
rownames(musclerna_late_resist_sig_tfout) <- musclerna_late_resist_sig_tfout$Motif.Name

musclerna_early_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/muscle_early_endur_sig_knownResults.txt",header = T,sep = "\t")
musclerna_early_endur_sig_tfout <- musclerna_early_endur_sig_tfout[!duplicated(musclerna_early_endur_sig_tfout$Motif.Name),]
rownames(musclerna_early_endur_sig_tfout) <- musclerna_early_endur_sig_tfout$Motif.Name

musclerna_mid_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/muscle_mid_endur_sig_knownResults.txt",header = T,sep = "\t")
musclerna_mid_endur_sig_tfout <- musclerna_mid_endur_sig_tfout[!duplicated(musclerna_mid_endur_sig_tfout$Motif.Name),]
rownames(musclerna_mid_endur_sig_tfout) <- musclerna_mid_endur_sig_tfout$Motif.Name

musclerna_late_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/muscle_late_endur_sig_knownResults.txt",header = T,sep = "\t")
musclerna_late_endur_sig_tfout <- musclerna_late_endur_sig_tfout[!duplicated(musclerna_late_endur_sig_tfout$Motif.Name),]
rownames(musclerna_late_endur_sig_tfout) <- musclerna_late_endur_sig_tfout$Motif.Name


# blood
bloodrna_immpost_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_immpost_resist_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_immpost_resist_sig_tfout <- bloodrna_immpost_resist_sig_tfout[!duplicated(bloodrna_immpost_resist_sig_tfout$Motif.Name),]
rownames(bloodrna_immpost_resist_sig_tfout) <- bloodrna_immpost_resist_sig_tfout$Motif.Name

bloodrna_early_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_early_resist_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_early_resist_sig_tfout <- bloodrna_early_resist_sig_tfout[!duplicated(bloodrna_early_resist_sig_tfout$Motif.Name),]
rownames(bloodrna_early_resist_sig_tfout) <- bloodrna_early_resist_sig_tfout$Motif.Name

bloodrna_mid_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_mid_resist_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_mid_resist_sig_tfout <- bloodrna_mid_resist_sig_tfout[!duplicated(bloodrna_mid_resist_sig_tfout$Motif.Name),]
rownames(bloodrna_mid_resist_sig_tfout) <- bloodrna_mid_resist_sig_tfout$Motif.Name

bloodrna_late_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_late_resist_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_late_resist_sig_tfout <- bloodrna_late_resist_sig_tfout[!duplicated(bloodrna_late_resist_sig_tfout$Motif.Name),]
rownames(bloodrna_late_resist_sig_tfout) <- bloodrna_late_resist_sig_tfout$Motif.Name

bloodrna_during20_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_during20_endur_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_during20_endur_sig_tfout <- bloodrna_during20_endur_sig_tfout[!duplicated(bloodrna_during20_endur_sig_tfout$Motif.Name),]
rownames(bloodrna_during20_endur_sig_tfout) <- bloodrna_during20_endur_sig_tfout$Motif.Name

bloodrna_during40_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_during40_endur_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_during40_endur_sig_tfout <- bloodrna_during40_endur_sig_tfout[!duplicated(bloodrna_during40_endur_sig_tfout$Motif.Name),]
rownames(bloodrna_during40_endur_sig_tfout) <- bloodrna_during40_endur_sig_tfout$Motif.Name

bloodrna_immpost_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_immpost_endur_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_immpost_endur_sig_tfout <- bloodrna_immpost_endur_sig_tfout[!duplicated(bloodrna_immpost_endur_sig_tfout$Motif.Name),]
rownames(bloodrna_immpost_endur_sig_tfout) <- bloodrna_immpost_endur_sig_tfout$Motif.Name

bloodrna_early_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_early_endur_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_early_endur_sig_tfout <- bloodrna_early_endur_sig_tfout[!duplicated(bloodrna_early_endur_sig_tfout$Motif.Name),]
rownames(bloodrna_early_endur_sig_tfout) <- bloodrna_early_endur_sig_tfout$Motif.Name

bloodrna_mid_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_mid_endur_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_mid_endur_sig_tfout <- bloodrna_mid_endur_sig_tfout[!duplicated(bloodrna_mid_endur_sig_tfout$Motif.Name),]
rownames(bloodrna_mid_endur_sig_tfout) <- bloodrna_mid_endur_sig_tfout$Motif.Name

bloodrna_late_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/blood_late_endur_sig_knownResults.txt",header = T,sep = "\t")
bloodrna_late_endur_sig_tfout <- bloodrna_late_endur_sig_tfout[!duplicated(bloodrna_late_endur_sig_tfout$Motif.Name),]
rownames(bloodrna_late_endur_sig_tfout) <- bloodrna_late_endur_sig_tfout$Motif.Name


# adipose
adiposerna_early_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/adipose_early_resist_sig_knownResults.txt",header = T,sep = "\t")
adiposerna_early_resist_sig_tfout <- adiposerna_early_resist_sig_tfout[!duplicated(adiposerna_early_resist_sig_tfout$Motif.Name),]
rownames(adiposerna_early_resist_sig_tfout) <- adiposerna_early_resist_sig_tfout$Motif.Name

adiposerna_mid_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/adipose_mid_resist_sig_knownResults.txt",header = T,sep = "\t")
adiposerna_mid_resist_sig_tfout <- adiposerna_mid_resist_sig_tfout[!duplicated(adiposerna_mid_resist_sig_tfout$Motif.Name),]
rownames(adiposerna_mid_resist_sig_tfout) <- adiposerna_mid_resist_sig_tfout$Motif.Name

adiposerna_late_resist_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/adipose_late_resist_sig_knownResults.txt",header = T,sep = "\t")
adiposerna_late_resist_sig_tfout <- adiposerna_late_resist_sig_tfout[!duplicated(adiposerna_late_resist_sig_tfout$Motif.Name),]
rownames(adiposerna_late_resist_sig_tfout) <- adiposerna_late_resist_sig_tfout$Motif.Name

adiposerna_early_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/adipose_early_endur_sig_knownResults.txt",header = T,sep = "\t")
adiposerna_early_endur_sig_tfout <- adiposerna_early_endur_sig_tfout[!duplicated(adiposerna_early_endur_sig_tfout$Motif.Name),]
rownames(adiposerna_early_endur_sig_tfout) <- adiposerna_early_endur_sig_tfout$Motif.Name

adiposerna_mid_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/adipose_mid_endur_sig_knownResults.txt",header = T,sep = "\t")
adiposerna_mid_endur_sig_tfout <- adiposerna_mid_endur_sig_tfout[!duplicated(adiposerna_mid_endur_sig_tfout$Motif.Name),]
rownames(adiposerna_mid_endur_sig_tfout) <- adiposerna_mid_endur_sig_tfout$Motif.Name

adiposerna_late_endur_sig_tfout <- read.delim("Precovid_DEG_HOMER_KnownTF_Results/adipose_late_endur_sig_knownResults.txt",header = T,sep = "\t")
adiposerna_late_endur_sig_tfout <- adiposerna_late_endur_sig_tfout[!duplicated(adiposerna_late_endur_sig_tfout$Motif.Name),]
rownames(adiposerna_late_endur_sig_tfout) <- adiposerna_late_endur_sig_tfout$Motif.Name


# The total list of known TFs generated by homer
tflist <- rownames(musclerna_early_resist_sig_tfout)

# Matrices of adjusted p-values (q-value) of TF enrichment for each tissue
musclerna_tfenrich_qvalmat <- cbind(musclerna_early_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                    musclerna_mid_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                    musclerna_late_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                    musclerna_early_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                    musclerna_mid_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                    musclerna_late_endur_sig_tfout[tflist,"q.value..Benjamini."])
colnames(musclerna_tfenrich_qvalmat) <- c("Resist 15/30/45 Min Post","Resist 3.5/4 Hr Post","Resist 24 Hr Post",
                                          "Endur 15/30/45 Min Post","Endur 3.5/4 Hr Post","Endur 24 Hr Post")
rownames(musclerna_tfenrich_qvalmat) <- tflist

adiposerna_tfenrich_qvalmat <- cbind(adiposerna_early_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                     adiposerna_mid_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                     adiposerna_late_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                     adiposerna_early_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                     adiposerna_mid_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                     adiposerna_late_endur_sig_tfout[tflist,"q.value..Benjamini."])
colnames(adiposerna_tfenrich_qvalmat) <- c("Resist 15/30/45 Min Post","Resist 3.5/4 Hr Post","Resist 24 Hr Post",
                                           "Endur 15/30/45 Min Post","Endur 3.5/4 Hr Post","Endur 24 Hr Post")
rownames(adiposerna_tfenrich_qvalmat) <- tflist

bloodrna_tfenrich_qvalmat <- cbind(bloodrna_immpost_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_early_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_mid_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_late_resist_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_during20_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_during40_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_immpost_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_early_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_mid_endur_sig_tfout[tflist,"q.value..Benjamini."],
                                   bloodrna_late_endur_sig_tfout[tflist,"q.value..Benjamini."])
colnames(bloodrna_tfenrich_qvalmat) <- c("Resist 10 Min Post","Resist 15/30/45 Min Post","Resist 3.5/4 Hr Post","Resist 24 Hr Post",
                                         "Endur 20 Min During","Endur 40 Min During","Endur 10 Min Post",
                                         "Endur 15/30/45 Min Post","Endur 3.5/4 Hr Post","Endur 24 Hr Post")
rownames(bloodrna_tfenrich_qvalmat) <- tflist

# Defining the colors for the figures
tf_ann_cols <- list("Tissue" = MotrpacHumanPreSuspension::HUMAN_TISSUE_COLORS,
                    "Sex" = MotrpacHumanPreSuspension::HUMAN_SEX_COLORS[c("Male","Female")],
                    "Timepoint" = MotrpacHumanPreSuspension::HUMAN_ACUTE_TIMEPOINT_COLORS[c("Pre_exercise","during_20_min","during_40_min","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr")],
                    "Modality" = MotrpacHumanPreSuspension::HUMAN_EXERCISE_GROUP_COLORS)
names(tf_ann_cols$Tissue) <- c("Blood","Muscle","Adipose")
names(tf_ann_cols$Timepoint) <- c("Pre","D20M","D40M","P10M","P15-45M","P3.5/4H","P24H")
names(tf_ann_cols$Modality) <- c("RE","EE","CON")

# The minimum non-zero muscle q value is 1E-4 so I will set the zero muscle q values at 1E-5

musclerna_tfenrich_qvalmat[musclerna_tfenrich_qvalmat == 0] <- 1e-05

bloodrna_tfenrich_qvalmat[bloodrna_tfenrich_qvalmat == 0] <- 1e-05

alltissue_tfenrich_qvalmat <- cbind(musclerna_tfenrich_qvalmat,adiposerna_tfenrich_qvalmat,bloodrna_tfenrich_qvalmat)

colnames(alltissue_tfenrich_qvalmat) <- c(paste("Muscle_",colnames(musclerna_tfenrich_qvalmat),sep = ""),
                                          paste("Adipose_",colnames(adiposerna_tfenrich_qvalmat),sep = ""),
                                          paste("Blood_",colnames(bloodrna_tfenrich_qvalmat),sep = ""))
alltissue_tfenrich_qvalmatcolumndf <- data.frame(row.names = colnames(alltissue_tfenrich_qvalmat),
                                                 "Tissue" = c(rep("Muscle",6),
                                                              rep("Adipose",6),
                                                              rep("Blood",10)),
                                                 "Modality" = c(rep("RE",3),
                                                                rep("EE",3),
                                                                rep("RE",3),
                                                                rep("EE",3),
                                                                rep("RE",4),
                                                                rep("EE",6)),
                                                 "Timepoint" = c("P15-45M",
                                                                 "P3.5/4H",
                                                                 "P24H",
                                                                 "P15-45M",
                                                                 "P3.5/4H",
                                                                 "P24H",
                                                                 "P15-45M",
                                                                 "P3.5/4H",
                                                                 "P24H",
                                                                 "P15-45M",
                                                                 "P3.5/4H",
                                                                 "P24H",
                                                                 "P10M",
                                                                 "P15-45M",
                                                                 "P3.5/4H",
                                                                 "P24H",
                                                                 "D20M",
                                                                 "D40M",
                                                                 "P10M",
                                                                 "P15-45M",
                                                                 "P3.5/4H",
                                                                 "P24H"))
alltissue_tfenrich_qvalmatcolumndf$Time.Point <- factor(alltissue_tfenrich_qvalmatcolumndf$Timepoint,
                                                        levels = c("D20M","D40M","P10M",
                                                                   "P15-45M","P3.5/4H","P24H"))
alltissue_tfenrich_qvalmatcolumndf$Tissue <- factor(alltissue_tfenrich_qvalmatcolumndf$Tissue,levels = c("Muscle","Adipose","Blood"))

#png(file = "All Tissue TF Enrich log10 Q Val Heatmap_062025.png",width = 8,height = 12,units = "in",res = 600)
#pheatmap(-log10(alltissue_tfenrich_qvalmat[apply(alltissue_tfenrich_qvalmat,1,min) < 0.01,]),labels_row = gsub("\\(.*","",rownames(alltissue_tfenrich_qvalmat[apply(alltissue_tfenrich_qvalmat,1,min) < 0.01,])),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 315,annotation_col = alltissue_tfenrich_qvalmatcolumndf[,c("Timepoint","Modality","Tissue")],annotation_colors = tf_ann_cols,show_colnames = F)
#dev.off()

colnames(musclerna_tfenrich_qvalmat) <- paste("Muscle_",colnames(musclerna_tfenrich_qvalmat),sep = "")
colnames(adiposerna_tfenrich_qvalmat) <- paste("Adipose_",colnames(adiposerna_tfenrich_qvalmat),sep = "")
colnames(bloodrna_tfenrich_qvalmat) <- paste("Blood_",colnames(bloodrna_tfenrich_qvalmat),sep = "")

# We want some statistics for our different TF enrichments
# First how many genes are taken as inputs for the different analyses

tfgeneinputcompdf <- data.frame("Sig.Count" = c(length(musclerna_early_resist_sig),
                                                length(musclerna_mid_resist_sig),
                                                length(musclerna_late_resist_sig),
                                                length(musclerna_early_endur_sig),
                                                length(musclerna_mid_endur_sig),
                                                length(musclerna_late_endur_sig),
                                                length(adiposerna_early_resist_sig),
                                                length(adiposerna_mid_resist_sig),
                                                length(adiposerna_late_resist_sig),
                                                length(adiposerna_early_endur_sig),
                                                length(adiposerna_mid_endur_sig),
                                                length(adiposerna_late_endur_sig),
                                                NA,
                                                NA,
                                                length(bloodrna_immpost_resist_sig),
                                                length(bloodrna_early_resist_sig),
                                                length(bloodrna_mid_resist_sig),
                                                length(bloodrna_late_resist_sig),
                                                length(bloodrna_during20_endur_sig),
                                                length(bloodrna_during40_endur_sig),
                                                length(bloodrna_immpost_endur_sig),
                                                length(bloodrna_early_endur_sig),
                                                length(bloodrna_mid_endur_sig),
                                                length(bloodrna_late_endur_sig)),
                                "Tissue" = c(rep("Muscle",6),
                                             rep("Adipose",6),
                                             rep("Blood",12)),
                                "Modality" = c(rep("RE",3),
                                               rep("EE",3),
                                               rep("RE",3),
                                               rep("EE",3),
                                               rep("RE",6),
                                               rep("EE",6)),
                                "Timepoint" = c("P15-45M",
                                                "P3.5/4H",
                                                "P24H",
                                                "P15-45M",
                                                "P3.5/4H",
                                                "P24H",
                                                "P15-45M",
                                                "P3.5/4H",
                                                "P24H",
                                                "P15-45M",
                                                "P3.5/4H",
                                                "P24H",
                                                "D20M",
                                                "D40M",
                                                "P10M",
                                                "P15-45M",
                                                "P3.5/4H",
                                                "P24H",
                                                "D20M",
                                                "D40M",
                                                "P10M",
                                                "P15-45M",
                                                "P3.5/4H",
                                                "P24H"))


tfgeneinputcompdf$Timepoint <- factor(tfgeneinputcompdf$Timepoint,levels = c("D20M","D40M","P10M","P15-45M","P3.5/4H","P24H"))
tfgeneinputcompdf$Tissue <- factor(tfgeneinputcompdf$Tissue,levels = c("Muscle","Adipose","Blood"))

# FIGURE 6A
png(file = "Figure6A_Counts of Sig Inputs for TF Analysis_102925.png",width = 5,height = 7,units = "in",res = 600)
ggbarplot(data = tfgeneinputcompdf,x = "Timepoint",y = "Sig.Count",fill = "Modality",color = "Modality",position = position_dodge(0.9)) + facet_wrap(~Tissue,ncol = 1) + theme(axis.text = element_text(angle = 45,hjust = 1)) + scale_fill_manual(values = tf_ann_cols$Modality) + ylab("Sig Count") + xlab("Time Point")
dev.off()

siglist <- list(musclerna_early_resist_sig,
                musclerna_mid_resist_sig,
                musclerna_late_resist_sig,
                musclerna_early_endur_sig,
                musclerna_mid_endur_sig,
                musclerna_late_endur_sig,
                adiposerna_early_resist_sig,
                adiposerna_mid_resist_sig,
                adiposerna_late_resist_sig,
                adiposerna_early_endur_sig,
                adiposerna_mid_endur_sig,
                adiposerna_late_endur_sig,
                bloodrna_immpost_resist_sig,
                bloodrna_early_resist_sig,
                bloodrna_mid_resist_sig,
                bloodrna_late_resist_sig,
                bloodrna_during20_endur_sig,
                bloodrna_during40_endur_sig,
                bloodrna_immpost_endur_sig,
                bloodrna_early_endur_sig,
                bloodrna_mid_endur_sig,
                bloodrna_late_endur_sig)


# Now we want to identify the top TFs for each comparison.
# 22 comparisons. Let's identify 5 per comparison. See how effective that is

mostsigtflist <- Reduce(union,list(musclerna_early_resist_sig_tfout[order(musclerna_early_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   musclerna_mid_resist_sig_tfout[order(musclerna_mid_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   musclerna_late_resist_sig_tfout[order(musclerna_late_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   musclerna_early_endur_sig_tfout[order(musclerna_early_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   musclerna_mid_endur_sig_tfout[order(musclerna_mid_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   musclerna_late_endur_sig_tfout[order(musclerna_late_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   adiposerna_early_resist_sig_tfout[order(adiposerna_early_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   adiposerna_mid_resist_sig_tfout[order(adiposerna_mid_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   adiposerna_late_resist_sig_tfout[order(adiposerna_late_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   adiposerna_early_endur_sig_tfout[order(adiposerna_early_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   adiposerna_mid_endur_sig_tfout[order(adiposerna_mid_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   adiposerna_late_endur_sig_tfout[order(adiposerna_late_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_immpost_resist_sig_tfout[order(bloodrna_immpost_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_early_resist_sig_tfout[order(bloodrna_early_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_mid_resist_sig_tfout[order(bloodrna_mid_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_late_resist_sig_tfout[order(bloodrna_late_resist_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_during20_endur_sig_tfout[order(bloodrna_during20_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_during40_endur_sig_tfout[order(bloodrna_during40_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_immpost_endur_sig_tfout[order(bloodrna_immpost_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_early_endur_sig_tfout[order(bloodrna_early_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_mid_endur_sig_tfout[order(bloodrna_mid_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5],
                                   bloodrna_late_endur_sig_tfout[order(bloodrna_late_endur_sig_tfout$q.value..Benjamini.),"Motif.Name"][1:5]))
# While this could have been a list of 22*5 = 110 TFs, there is considerable overlap of top enriched TFs so this is only 46 TFs long.

# Figure 6C
png(file = "Figure6C_All Tissue Select TF Enrich log10 Q Val Heatmap_102925.png",width = 5,height = 7,units = "in",res = 600)
pheatmap(-log10(alltissue_tfenrich_qvalmat[mostsigtflist,c(10,11,12,7,8,9,17,18,19,20,21,22,13,14,15,16,4,5,6,1,2,3)]),labels_row = gsub("\\(.*","",mostsigtflist),cluster_cols = F,color = colorpanel(101,"white","firebrick"),angle_col = 315,annotation_col = alltissue_tfenrich_qvalmatcolumndf[,c("Timepoint","Modality","Tissue")],annotation_colors = tf_ann_cols,show_colnames = F)
dev.off()

# Now we want to do a series of dotplots to show the variability in TF enrichment across each comparison

alltissue_tfenrich_qvaldf <- data.frame("qval" = as.vector(alltissue_tfenrich_qvalmat),
                                        "Comparison" = rep(colnames(alltissue_tfenrich_qvalmat),each = dim(alltissue_tfenrich_qvalmat)[1]))
alltissue_tfenrich_qvaldf$Tissue <- alltissue_tfenrich_qvalmatcolumndf[alltissue_tfenrich_qvaldf$Comparison,"Tissue"]
alltissue_tfenrich_qvaldf$Modality <- alltissue_tfenrich_qvalmatcolumndf[alltissue_tfenrich_qvaldf$Comparison,"Modality"]
alltissue_tfenrich_qvaldf$Timepoint <- alltissue_tfenrich_qvalmatcolumndf[alltissue_tfenrich_qvaldf$Comparison,"Timepoint"]

alltissue_tfenrich_qvaldf$log10qval <- -log10(alltissue_tfenrich_qvaldf$qval)

# We generate a more thorough Homer TF annotation list from a previous list used for PASS1B
human_feature_to_gene <- as.data.frame(MotrpacHumanPreSuspensionData::HUMAN_FEATURE_TO_GENE)

tfproanno <- readRDS("tfproanno.RDS")

tfproanno <- tfproanno[,c("Gene.Name","Ensembl")]
tfproanno$Gene.Name <- toupper(tfproanno$Gene.Name)

tfproanno[11,"Gene.Name"] <- "JUN"
tfproanno[42,"Gene.Name"] <- "LHX6"
tfproanno[23,"Gene.Name"] <- "ETS1"
tfproanno[33,"Gene.Name"] <- "EWSR1"
tfproanno[104,"Gene.Name"] <- "IRF8"
tfproanno[137,"Gene.Name"] <- "STAT3"
tfproanno[175,"Gene.Name"] <- "NANOG"
tfproanno[190,"Gene.Name"] <- "POU3F1"
tfproanno[205,"Gene.Name"] <- "POU5F1"
tfproanno[220,"Gene.Name"] <- "RUNX1"
tfproanno[320,"Gene.Name"] <- "BATF"
tfproanno[337,"Gene.Name"] <- "NFKB2"
tfproanno[379,"Gene.Name"] <- "DMC1"
tfproanno[406,"Gene.Name"] <- "AR"
tfproanno[435,"Gene.Name"] <- "CTCF"


newtfproanno <- data.frame(row.names = rownames(alltissue_tfenrich_qvalmat),
                           "Gene.Name" = toupper(gsub("\\(.*","",rownames(alltissue_tfenrich_qvalmat))),
                           "Ensembl" = rep("",length(rownames(alltissue_tfenrich_qvalmat))))
for(i in 1:dim(newtfproanno)[1]){
  if(rownames(newtfproanno)[i] %in% rownames(tfproanno)){
    newtfproanno[i,"Gene.Name"] <- tfproanno[rownames(newtfproanno)[i],"Gene.Name"]
  }
  if(newtfproanno[i,"Gene.Name"] %in% human_feature_to_gene$gene_symbol){
    newtfproanno[i,"Ensembl"] <- as.character(human_feature_to_gene[human_feature_to_gene$gene_symbol %in% newtfproanno[i,"Gene.Name"],"ensembl_gene"][1])
  }
}
newtfproanno$prot_pr_id <- ""
for(i in 1:dim(newtfproanno)[1]){
  ourtfgene <- newtfproanno[i,"Gene.Name"]
  if(ourtfgene %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-pr","gene_symbol"]){
    newtfproanno[i,"prot_pr_id"] <- as.character(human_feature_to_gene[human_feature_to_gene$gene_symbol %in% ourtfgene & human_feature_to_gene$assay %in% "prot-pr","feature_id"][[1]])
  }
}
tfproanno <- newtfproanno
rm(newtfproanno)

tfproanno$prot_ol_id <- ""
for(i in 1:dim(tfproanno)[1]){
  ourtfgene <- tfproanno[i,"Gene.Name"]
  if(ourtfgene %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ol","gene_symbol"]){
    tfproanno[i,"prot_ol_id"] <- as.character(human_feature_to_gene[human_feature_to_gene$gene_symbol %in% ourtfgene & human_feature_to_gene$assay %in% "prot-ol","feature_id"][[1]])
  }
}

# We identify a subset of 381 TFs in the tflist where their corresponding gene symbol is present
# in the human data

tftrim <- intersect(tflist,rownames(tfproanno))
tftrimmed <- rownames(tfproanno)[rownames(tfproanno) %in% tftrim & toupper(tfproanno$Gene.Name) %in% human_feature_to_gene$gene_symbol]
tftrimmedgenesym <- toupper(tfproanno[tftrimmed,"Gene.Name"])

# We can identify the TFs present in the protein abundance, phosphorylation and olink data for the different tissues

tfmuscleprotpr <- tftrimmedgenesym[tftrimmedgenesym %in% MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_QC$feature_metadata$gene_symbol]
tfmuscleprotph <- tftrimmedgenesym[tftrimmedgenesym %in% MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_QC$feature_metadata$gene_symbol]
tfadiposeprotpr <- tftrimmedgenesym[tftrimmedgenesym %in% MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_QC$feature_metadata$gene_symbol]
tfadiposeprotph <- tftrimmedgenesym[tftrimmedgenesym %in% MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_QC$feature_metadata$gene_symbol]
tfbloodprotol <- tftrimmedgenesym[tftrimmedgenesym %in% MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_QC$feature_metadata$assay]

tfmuscleprotprids <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_QC$feature_metadata[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_QC$feature_metadata$gene_symbol %in% tfmuscleprotpr,"protein_id"]
tfmuscleprotphids <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_QC$feature_metadata[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_QC$feature_metadata$gene_symbol %in% tfmuscleprotph,"id"]

tfmusclernaids <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC$feature_metadata[gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC$feature_metadata$feature_id) %in% tfproanno$Ensembl,"feature_id"]
tfadiposernaids <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC$feature_metadata[gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC$feature_metadata$feature_id) %in% tfproanno$Ensembl,"feature_id"]
tfbloodrnaids <- MotrpacHumanPreSuspensionData::blood_transcript_1$feature_metadata[gsub("\\..*","",MotrpacHumanPreSuspensionData::blood_transcript_1$feature_metadata$feature_id) %in% tfproanno$Ensembl,"feature_id"]


tfadiposeprotprids <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_QC$feature_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_QC$feature_metadata$gene_symbol %in% tfadiposeprotpr,"protein_id"]
tfadiposeprotprids <- intersect(tfadiposeprotprids,MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id)
tfadiposeprotphids <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_QC$feature_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_QC$feature_metadata$gene_symbol %in% tfadiposeprotph,"id"]
tfadiposeprotphids <- intersect(tfadiposeprotphids,MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$feature_id)

tfbloodprotolids <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_QC$feature_metadata[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_QC$feature_metadata$assay %in% tfbloodprotol,"feature_id"]


# We want to list the exercise responses for each of the TFs in our database - First at a L2FC and pval level

# Muscle
muscletfprotprpval <- matrix(0L,nrow = length(tfmuscleprotprids),ncol = 6)
rownames(muscletfprotprpval) <- tfmuscleprotprids
colnames(muscletfprotprpval) <- colnames(musclerna_tfenrich_qvalmat)

muscletfprotprl2fc <- matrix(0L,nrow = length(tfmuscleprotprids),ncol = 6)
rownames(muscletfprotprl2fc) <- tfmuscleprotprids
colnames(muscletfprotprl2fc) <- colnames(musclerna_tfenrich_qvalmat)

muscletfprotprmeta <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_QC$feature_metadata[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_QC$feature_metadata$protein_id %in% tfmuscleprotprids,]
rownames(muscletfprotprmeta) <- muscletfprotprmeta$protein_id
muscletfprotprmeta <- muscletfprotprmeta[tfmuscleprotprids,]

for(i in 1:dim(muscletfprotprpval)[1]){
  ourid <- rownames(muscletfprotprpval)[i]
  muscletfprotprpval[i,"Muscle_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                                MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotprl2fc[i,"Muscle_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                                MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotprpval[i,"Muscle_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotprl2fc[i,"Muscle_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotprpval[i,"Muscle_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotprl2fc[i,"Muscle_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotprpval[i,"Muscle_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotprl2fc[i,"Muscle_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotprpval[i,"Muscle_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                           MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotprl2fc[i,"Muscle_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                           MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotprpval[i,"Muscle_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotprl2fc[i,"Muscle_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]

}


muscletfprotprpvaldisplay <- muscletfprotprpval
rownames(muscletfprotprpvaldisplay) <- muscletfprotprmeta$gene_symbol

muscletfprotprl2fcdisplay <- muscletfprotprl2fc
rownames(muscletfprotprl2fcdisplay) <- muscletfprotprmeta$gene_symbol

tfmuscleprotphidstrim <- tfmuscleprotphids[tfmuscleprotphids %in% MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id]

muscletfprotphpval <- matrix(0L,nrow = length(tfmuscleprotphidstrim),ncol = 6)
rownames(muscletfprotphpval) <- tfmuscleprotphidstrim
colnames(muscletfprotphpval) <- colnames(musclerna_tfenrich_qvalmat)

muscletfprotphl2fc <- matrix(0L,nrow = length(tfmuscleprotphidstrim),ncol = 6)
rownames(muscletfprotphl2fc) <- tfmuscleprotphidstrim
colnames(muscletfprotphl2fc) <- colnames(musclerna_tfenrich_qvalmat)

muscletfprotphmeta <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_QC$feature_metadata[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_QC$feature_metadata$id %in% tfmuscleprotphids,]
rownames(muscletfprotphmeta) <- muscletfprotphmeta$id
muscletfprotphmeta <- muscletfprotphmeta[tfmuscleprotphidstrim,]

for(i in 1:dim(muscletfprotphpval)[1]){
  ourid <- rownames(muscletfprotphpval)[i]
  muscletfprotphpval[i,"Muscle_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                                MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotphl2fc[i,"Muscle_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                                MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotphpval[i,"Muscle_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotphl2fc[i,"Muscle_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotphpval[i,"Muscle_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotphl2fc[i,"Muscle_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotphpval[i,"Muscle_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotphl2fc[i,"Muscle_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotphpval[i,"Muscle_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                           MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotphl2fc[i,"Muscle_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                           MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfprotphpval[i,"Muscle_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfprotphl2fc[i,"Muscle_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]

}

# For phosphorylation, we need to merge the gene_symbol with the phosphosite to make the appropriate name

muscletfprotphpvaldisplay <- muscletfprotphpval
rownames(muscletfprotphpvaldisplay) <- paste(muscletfprotphmeta$gene_symbol,gsub(".*_","",rownames(muscletfprotphpval)),sep = "_")

muscletfprotphl2fcdisplay <- muscletfprotphl2fc
rownames(muscletfprotphl2fcdisplay) <- paste(muscletfprotphmeta$gene_symbol,gsub(".*_","",rownames(muscletfprotphl2fc)),sep = "_")

# Muscle RNA
muscletfrnapval <- matrix(0L,nrow = length(tfmusclernaids ),ncol = 6)
rownames(muscletfrnapval) <- tfmusclernaids
colnames(muscletfrnapval) <- colnames(musclerna_tfenrich_qvalmat)

muscletfrnal2fc <- matrix(0L,nrow = length(tfmusclernaids ),ncol = 6)
rownames(muscletfrnal2fc) <- tfmusclernaids
colnames(muscletfrnal2fc) <- colnames(musclerna_tfenrich_qvalmat)

muscletfrnameta <- data.frame(row.names = MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC$feature_metadata[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC$feature_metadata$feature_id %in% tfmusclernaids,],
                              "feature_id" = MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC$feature_metadata[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC$feature_metadata$feature_id %in% tfmusclernaids,])



for(i in 1:dim(muscletfrnapval)[1]){
  ourid <- rownames(muscletfrnapval)[i]
  muscletfrnapval[i,"Muscle_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                              MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfrnal2fc[i,"Muscle_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                              MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfrnapval[i,"Muscle_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                          MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfrnal2fc[i,"Muscle_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                          MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfrnapval[i,"Muscle_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                       MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfrnal2fc[i,"Muscle_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                       MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfrnapval[i,"Muscle_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                             MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfrnal2fc[i,"Muscle_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                             MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfrnapval[i,"Muscle_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfrnal2fc[i,"Muscle_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  muscletfrnapval[i,"Muscle_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  muscletfrnal2fc[i,"Muscle_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]

}

muscletfrnameta$gene_symbol <- ""
for(i in 1:dim(muscletfrnameta)[1]){
  ourtf <- muscletfrnameta[i,"feature_id"]
  muscletfrnameta[i,"gene_symbol"] <- as.character(human_feature_to_gene[human_feature_to_gene$feature_id %in% ourtf,"gene_symbol"][1])
}

muscletfrnapvaldisplay <- muscletfrnapval
rownames(muscletfrnapvaldisplay) <- muscletfrnameta$gene_symbol

muscletfrnal2fcdisplay <- muscletfrnal2fc
rownames(muscletfrnal2fcdisplay) <- muscletfrnameta$gene_symbol


# Adipose RNA
adiposetfrnapval <- matrix(0L,nrow = length(tfadiposernaids ),ncol = 6)
rownames(adiposetfrnapval) <- tfadiposernaids
colnames(adiposetfrnapval) <- colnames(adiposerna_tfenrich_qvalmat)

adiposetfrnal2fc <- matrix(0L,nrow = length(tfadiposernaids ),ncol = 6)
rownames(adiposetfrnal2fc) <- tfadiposernaids
colnames(adiposetfrnal2fc) <- colnames(adiposerna_tfenrich_qvalmat)

adiposetfrnameta <- data.frame(row.names = MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC$feature_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC$feature_metadata$feature_id %in% tfadiposernaids,],
                               "feature_id" = MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC$feature_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC$feature_metadata$feature_id %in% tfadiposernaids,])



for(i in 1:dim(adiposetfrnapval)[1]){
  ourid <- rownames(adiposetfrnapval)[i]
  adiposetfrnapval[i,"Adipose_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                                 MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfrnal2fc[i,"Adipose_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                                 MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfrnapval[i,"Adipose_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                             MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfrnal2fc[i,"Adipose_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                             MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfrnapval[i,"Adipose_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                          MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfrnal2fc[i,"Adipose_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                          MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfrnapval[i,"Adipose_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                                MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfrnal2fc[i,"Adipose_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                                MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfrnapval[i,"Adipose_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfrnal2fc[i,"Adipose_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfrnapval[i,"Adipose_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfrnal2fc[i,"Adipose_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]

}

adiposetfrnameta$gene_symbol <- ""
for(i in 1:dim(adiposetfrnameta)[1]){
  ourtf <- adiposetfrnameta[i,"feature_id"]
  adiposetfrnameta[i,"gene_symbol"] <- as.character(human_feature_to_gene[human_feature_to_gene$feature_id %in% ourtf,"gene_symbol"][1])
}

adiposetfrnapvaldisplay <- adiposetfrnapval
rownames(adiposetfrnapvaldisplay) <- adiposetfrnameta$gene_symbol

adiposetfrnal2fcdisplay <- adiposetfrnal2fc
rownames(adiposetfrnal2fcdisplay) <- adiposetfrnameta$gene_symbol


# Blood RNA
bloodtfrnapval <- matrix(0L,nrow = length(tfbloodrnaids ),ncol = 10)
rownames(bloodtfrnapval) <- tfbloodrnaids
colnames(bloodtfrnapval) <- colnames(bloodrna_tfenrich_qvalmat)

bloodtfrnal2fc <- matrix(0L,nrow = length(tfbloodrnaids ),ncol = 10)
rownames(bloodtfrnal2fc) <- tfbloodrnaids
colnames(bloodtfrnal2fc) <- colnames(bloodrna_tfenrich_qvalmat)

bloodtfrnameta <- data.frame(row.names = MotrpacHumanPreSuspensionData::blood_transcript_1$feature_metadata[MotrpacHumanPreSuspensionData::blood_transcript_1$feature_metadata$feature_id %in% tfbloodrnaids,],
                             "feature_id" = MotrpacHumanPreSuspensionData::blood_transcript_1$feature_metadata[MotrpacHumanPreSuspensionData::blood_transcript_1$feature_metadata$feature_id %in% tfbloodrnaids,])



for(i in 1:dim(bloodtfrnapval)[1]){
  ourid <- rownames(bloodtfrnapval)[i]
  bloodtfrnapval[i,"Blood_Resist 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                     MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_10_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Resist 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                     MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_10_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                           MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                           MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                       MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                       MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                    MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                    MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Endur 20 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_20_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Endur 20 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_20_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Endur 40 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.during_40_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_40_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Endur 40 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.during_40_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_40_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Endur 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                    MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_10_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Endur 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                    MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                          MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                          MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfrnapval[i,"Blood_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                   MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfrnal2fc[i,"Blood_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id %in% ourid &
                                                                                                   MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]

}

bloodtfrnameta$gene_symbol <- ""
for(i in 1:dim(bloodtfrnameta)[1]){
  ourtf <- bloodtfrnameta[i,"feature_id"]
  bloodtfrnameta[i,"gene_symbol"] <- as.character(human_feature_to_gene[human_feature_to_gene$feature_id %in% ourtf,"gene_symbol"][1])
}

bloodtfrnapvaldisplay <- bloodtfrnapval
rownames(bloodtfrnapvaldisplay) <- bloodtfrnameta$gene_symbol

bloodtfrnal2fcdisplay <- bloodtfrnal2fc
rownames(bloodtfrnal2fcdisplay) <- bloodtfrnameta$gene_symbol

# Adipose Prot-pr
adiposetfprotprpval <- matrix(0L,nrow = length(tfadiposeprotprids),ncol = 2)
rownames(adiposetfprotprpval) <- tfadiposeprotprids
colnames(adiposetfprotprpval) <- colnames(adiposerna_tfenrich_qvalmat)[c(2,5)]

adiposetfprotprl2fc <- matrix(0L,nrow = length(tfadiposeprotprids),ncol = 2)
rownames(adiposetfprotprl2fc) <- tfadiposeprotprids
colnames(adiposetfprotprl2fc) <- colnames(adiposerna_tfenrich_qvalmat)[c(2,5)]

adiposetfprotprmeta <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_QC$feature_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_QC$feature_metadata$protein_id %in% tfadiposeprotprids,]
rownames(adiposetfprotprmeta) <- adiposetfprotprmeta$protein_id
adiposetfprotprmeta <- adiposetfprotprmeta[tfadiposeprotprids,]

for(i in 1:dim(adiposetfprotprpval)[1]){
  ourid <- rownames(adiposetfprotprpval)[i]
  adiposetfprotprpval[i,"Adipose_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfprotprl2fc[i,"Adipose_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfprotprpval[i,"Adipose_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                              MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfprotprl2fc[i,"Adipose_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id %in% ourid &
                                                                                                              MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
}


adiposetfprotprpvaldisplay <- adiposetfprotprpval
rownames(adiposetfprotprpvaldisplay) <- adiposetfprotprmeta$gene_symbol

adiposetfprotprl2fcdisplay <- adiposetfprotprl2fc
rownames(adiposetfprotprl2fcdisplay) <- adiposetfprotprmeta$gene_symbol

# Adipose Prot-ph
tfadiposeprotphidstrim <- tfadiposeprotphids[tfadiposeprotphids %in% MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$feature_id]

adiposetfprotphpval <- matrix(0L,nrow = length(tfadiposeprotphidstrim),ncol = 2)
rownames(adiposetfprotphpval) <- tfadiposeprotphidstrim
colnames(adiposetfprotphpval) <- colnames(adiposerna_tfenrich_qvalmat)[c(2,5)]

adiposetfprotphl2fc <- matrix(0L,nrow = length(tfadiposeprotphidstrim),ncol = 2)
rownames(adiposetfprotphl2fc) <- tfadiposeprotphidstrim
colnames(adiposetfprotphl2fc) <- colnames(adiposerna_tfenrich_qvalmat)[c(2,5)]

adiposetfprotphmeta <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_QC$feature_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_QC$feature_metadata$id %in% tfadiposeprotphids,]
rownames(adiposetfprotphmeta) <- adiposetfprotphmeta$id
adiposetfprotphmeta <- adiposetfprotphmeta[tfadiposeprotphidstrim,]

for(i in 1:dim(adiposetfprotphpval)[1]){
  ourid <- rownames(adiposetfprotphpval)[i]
  adiposetfprotphpval[i,"Adipose_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfprotphl2fc[i,"Adipose_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                               MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  adiposetfprotphpval[i,"Adipose_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                              MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  adiposetfprotphl2fc[i,"Adipose_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$feature_id %in% ourid &
                                                                                                              MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
}

# For phosphorylation, we need to merge the gene_symbol with the phosphosite to make the appropriate name

adiposetfprotphpvaldisplay <- adiposetfprotphpval
rownames(adiposetfprotphpvaldisplay) <- paste(adiposetfprotphmeta$gene_symbol,gsub(".*_","",rownames(adiposetfprotphpval)),sep = "_")

adiposetfprotphl2fcdisplay <- adiposetfprotphl2fc
rownames(adiposetfprotphl2fcdisplay) <- paste(adiposetfprotphmeta$gene_symbol,gsub(".*_","",rownames(adiposetfprotphl2fc)),sep = "_")


# Blood Prot-ol

bloodtfprotolpval <- matrix(0L,nrow = length(tfbloodprotolids),ncol = 10)
rownames(bloodtfprotolpval) <- tfbloodprotolids
colnames(bloodtfprotolpval) <- colnames(bloodrna_tfenrich_qvalmat)

bloodtfprotoll2fc <- matrix(0L,nrow = length(tfbloodprotolids),ncol = 10)
rownames(bloodtfprotoll2fc) <- tfbloodprotolids
colnames(bloodtfprotoll2fc) <- colnames(bloodrna_tfenrich_qvalmat)

bloodtfprotolmeta <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_QC$feature_metadata[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_QC$feature_metadata$feature_id %in% tfbloodprotolids,]
rownames(bloodtfprotolmeta) <- bloodtfprotolmeta$feature_id
bloodtfprotolmeta <- bloodtfprotolmeta[tfbloodprotolids,]

for(i in 1:dim(bloodtfprotolpval)[1]){
  ourid <- rownames(bloodtfprotolpval)[i]
  bloodtfprotolpval[i,"Blood_Resist 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                       MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_10_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Resist 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                       MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_10_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                             MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                             MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                         MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Endur 20 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_20_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Endur 20 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_20_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Endur 40 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.during_40_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_40_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Endur 40 Min During"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.during_40_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_40_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Endur 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_10_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Endur 10 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                      MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_10_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                            MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]
  bloodtfprotolpval[i,"Blood_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                     MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  bloodtfprotoll2fc[i,"Blood_Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id %in% ourid &
                                                                                                     MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","logFC"][[1]]

}

bloodtfprotolpvaldisplay <- bloodtfprotolpval
rownames(bloodtfprotolpvaldisplay) <- bloodtfprotolmeta$assay

bloodtfprotoll2fcdisplay <- bloodtfprotoll2fc
rownames(bloodtfprotoll2fcdisplay) <- bloodtfprotolmeta$assay

muscletfprotprmeta$TF <- ""
for(i in 1:dim(muscletfprotprmeta)[1]){
  ourid <- rownames(muscletfprotprmeta)[i]
  muscletfprotprmeta[i,"TF"] <- rownames(tfproanno[toupper(tfproanno$Gene.Name) %in% muscletfprotprmeta[ourid,"gene_symbol"],])[1]
}
muscletfprotprmeta$TFtrim <- gsub("\\(.*","",muscletfprotprmeta$TF)


#####
# We next want to create a dataframe where we connect a TF's enrichment, differential gene expression, phosphorylation
####

# Let's start with muscle

muscle_tf_summarymat <- matrix(1L,nrow = dim(alltissue_tfenrich_qvalmat)[1],ncol = 18)
rownames(muscle_tf_summarymat) <- rownames(alltissue_tfenrich_qvalmat)
colnames(muscle_tf_summarymat) <- c("Enrich Endur 15/30/45 Min Post",
                                    "Enrich Endur 3.5/4 Hr Post",
                                    "Enrich Endur 24 Hr Post",
                                    "Enrich Resist 15/30/45 Min Post",
                                    "Enrich Resist 3.5/4 Hr Post",
                                    "Enrich Resist 24 Hr Post",
                                    "TF Gene Endur 15/30/45 Min Post",
                                    "TF Gene Endur 3.5/4 Hr Post",
                                    "TF Gene Endur 24 Hr Post",
                                    "TF Gene Resist 15/30/45 Min Post",
                                    "TF Gene Resist 3.5/4 Hr Post",
                                    "TF Gene Resist 24 Hr Post",
                                    "TF Phos Endur 15/30/45 Min Post",
                                    "TF Phos Endur 3.5/4 Hr Post",
                                    "TF Phos Endur 24 Hr Post",
                                    "TF Phos Resist 15/30/45 Min Post",
                                    "TF Phos Resist 3.5/4 Hr Post",
                                    "TF Phos Resist 24 Hr Post")
muscle_tf_summarymat <- -1*muscle_tf_summarymat
for(i in 1:dim(muscle_tf_summarymat)[1]){
  if(i%%20 == 0){
    print(i)
  }

  ourtf <- rownames(muscle_tf_summarymat)[i]

  muscle_tf_summarymat[i,"Enrich Endur 15/30/45 Min Post"] <- alltissue_tfenrich_qvalmat[ourtf,"Muscle_Endur 15/30/45 Min Post"]
  muscle_tf_summarymat[i,"Enrich Endur 3.5/4 Hr Post"] <- alltissue_tfenrich_qvalmat[ourtf,"Muscle_Endur 3.5/4 Hr Post"]
  muscle_tf_summarymat[i,"Enrich Endur 24 Hr Post"] <- alltissue_tfenrich_qvalmat[ourtf,"Muscle_Endur 24 Hr Post"]
  muscle_tf_summarymat[i,"Enrich Resist 15/30/45 Min Post"] <- alltissue_tfenrich_qvalmat[ourtf,"Muscle_Resist 15/30/45 Min Post"]
  muscle_tf_summarymat[i,"Enrich Resist 3.5/4 Hr Post"] <- alltissue_tfenrich_qvalmat[ourtf,"Muscle_Resist 3.5/4 Hr Post"]
  muscle_tf_summarymat[i,"Enrich Resist 24 Hr Post"] <- alltissue_tfenrich_qvalmat[ourtf,"Muscle_Resist 24 Hr Post"]

  if(tfproanno[ourtf,"Ensembl"] %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id)){
    ourens <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id[grep(tfproanno[ourtf,"Ensembl"],MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id)][1]
    muscle_tf_summarymat[i,"TF Gene Endur 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourens &
                                                                                                                     MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
    muscle_tf_summarymat[i,"TF Gene Endur 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourens &
                                                                                                                 MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
    muscle_tf_summarymat[i,"TF Gene Endur 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourens &
                                                                                                              MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
    muscle_tf_summarymat[i,"TF Gene Resist 15/30/45 Min Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourens &
                                                                                                                      MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
    muscle_tf_summarymat[i,"TF Gene Resist 3.5/4 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourens &
                                                                                                                  MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
    muscle_tf_summarymat[i,"TF Gene Resist 24 Hr Post"] <- MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id %in% ourens &
                                                                                                               MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"][[1]]
  }
  if("prot-ph" %in% names(table(human_feature_to_gene[human_feature_to_gene$gene_symbol %in% tfproanno[ourtf,"Gene.Name"],"assay"]))){

    ourphossites <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph" & human_feature_to_gene$gene_symbol %in% tfproanno[ourtf,"Gene.Name"],"feature_id"]
    if(length(intersect(ourphossites,MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id)) > 0){

      muscle_tf_summarymat[i,"TF Phos Endur 15/30/45 Min Post"] <- min(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourphossites &
                                                                                                                          MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"])
      muscle_tf_summarymat[i,"TF Phos Endur 3.5/4 Hr Post"] <- min(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourphossites &
                                                                                                                      MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"])
      muscle_tf_summarymat[i,"TF Phos Endur 24 Hr Post"] <- min(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourphossites &
                                                                                                                   MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"])
      muscle_tf_summarymat[i,"TF Phos Resist 15/30/45 Min Post"] <- min(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourphossites &
                                                                                                                           MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise","adj_p_value"])
      muscle_tf_summarymat[i,"TF Phos Resist 3.5/4 Hr Post"] <- min(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourphossites &
                                                                                                                       MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise","adj_p_value"])
      muscle_tf_summarymat[i,"TF Phos Resist 24 Hr Post"] <- min(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$feature_id %in% ourphossites &
                                                                                                                    MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise","adj_p_value"])

    }
  }
}

muscle_tf_summarymatlog <- muscle_tf_summarymat
for(i in 1:dim(muscle_tf_summarymatlog)[1]){
  for(j in 1:dim(muscle_tf_summarymatlog)[2]){
    if(muscle_tf_summarymatlog[i,j] == -1){
      muscle_tf_summarymatlog[i,j] <- -5
    }
    if(muscle_tf_summarymatlog[i,j] > -1){
      muscle_tf_summarymatlog[i,j] <- -1*log(muscle_tf_summarymatlog[i,j])
    }
  }
}


muscle_tf_summarymatlog_sig <- muscle_tf_summarymatlog[apply(muscle_tf_summarymatlog[,c(1:6)],1,max) > 2.995732,]

muscle_tf_summarymatlog_sig_phossig <- muscle_tf_summarymatlog_sig[apply(muscle_tf_summarymatlog_sig[,c(13:18)],1,max) >= 2.995732,]

rna_human_feature_to_gene <- human_feature_to_gene[human_feature_to_gene$platform %in% "transcript-rna-seq",]
rownames(rna_human_feature_to_gene) <- rna_human_feature_to_gene$feature_id

# Figure 6H
png("Figure6H_Muscle TFs enriched for DEGs and showing phospho changes Reaarranged_trimmedlegend_102925.png",width = 6,height = 6,units = "in",res = 600)
pheatmap(muscle_tf_summarymatlog_sig[apply(muscle_tf_summarymatlog_sig[,c(13:18)],1,max) >= 2.995732,c(13:18,1:6)],breaks = seq(0,5,length.out = 101),color = colorpanel(101,"white","firebrick"),angle_col = 315,labels_row = toupper(gsub("\\(.*","",rownames(muscle_tf_summarymatlog_sig[apply(muscle_tf_summarymatlog_sig[,c(13:18)],1,max) >= 2.995732,]))),cluster_cols = F)
dev.off()


#####
# Now we generate a bubbleHeatmap of TFs responding to exercise in the muscle phosphoproteomics data
####


# I generate a tissue-specific heatmap for the color annotation on top of the bubbleheatmap

png(filename = "muscleheatmap_topannotationtemplateforbubble_102925.png",width = 8,height = 8,units = "in",res = 600)
pheatmap(-log10(muscletfprotphpvaldisplay[apply(muscletfprotphpvaldisplay,1,min) < 0.05,]),cluster_cols = F,breaks = seq(0,3,length.out = 101),color = colorpanel(101,"white","firebrick"),display_numbers = T,show_colnames = F,annotation_col = alltissue_tfenrich_qvalmatcolumndf[,c("Timepoint","Modality")],annotation_colors = tf_ann_cols,cellwidth = 10)
dev.off()

rownames(muscletfprotphl2fcdisplay) <- gsub("s|t|y","",rownames(muscletfprotphl2fcdisplay))
muscletfprotphl2fcdisplay <- muscletfprotphl2fcdisplay[order(rownames(muscletfprotphl2fcdisplay)),]

rownames(muscletfprotphpvaldisplay) <- gsub("s|t|y","",rownames(muscletfprotphpvaldisplay))
muscletfprotphpvaldisplay <- muscletfprotphpvaldisplay[order(rownames(muscletfprotphpvaldisplay)),]

testout <- bubbleHeatmap(colorMat = muscletfprotphl2fcdisplay[apply(muscletfprotphpvaldisplay,1,min) < 0.05,c(4,5,6,1,2,3)],sizeMat = -log10(muscletfprotphpvaldisplay[apply(muscletfprotphpvaldisplay,1,min) < 0.05,c(4,5,6,1,2,3)]),legendTitles = c("Significance","log2FC"),colorSeq = c("blue","white","red"),colorLim = c(-1,1),sizeLim = c(0,5))


# Filter rows

#--------------------------------------------------------------------------------------
#Chris: This code is to replace the above bubbleheatmap code to keep consistent
#with all of the other plots and to avoid having to combine different figures
#in illustrator. - Nov 17th, 2025.
# Figure 6D

keep_rows = apply(muscletfprotphpvaldisplay, 1, min) < 0.05
colorMat = muscletfprotphl2fcdisplay[keep_rows, c(4,5,6,1,2,3)]
sizeMat  = -log10(muscletfprotphpvaldisplay[keep_rows, c(4,5,6,1,2,3)])

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
        grepl("3.5", tissue_mod_timepoint, ignore.case = TRUE) ~ "P3.5/4H",
        grepl("24", tissue_mod_timepoint, ignore.case = TRUE) ~ "P24H",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(Tissue, Modality, Timepoint, tissue_mod_timepoint, value)
}

anno_df <- annotate_matrix(sizeMat) %>%
  dplyr::select(Modality, Timepoint) %>%
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
  )
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

pdf(file = "6D_muscle_TF_protph_bubbleheatmap_chris.pdf", height = 12, width = 5)
ComplexHeatmap::draw(
  ht,
  annotation_legend_list = NULL,
  heatmap_legend_list = list(color_lg, size_lg),
  merge_legend = TRUE
)
dev.off()
#--------------------------------------------------------------------------------------

# We next make a combined TF RNA bubbleheatmap from the three tissues
# First we want to make a combined matrix for the three tissues
tissuecombornal2fcdisplay <- matrix(0L,nrow = length(Reduce(union,list(rownames(bloodtfrnal2fcdisplay),rownames(adiposetfrnal2fcdisplay),rownames(muscletfrnal2fcdisplay)))),ncol = 22)
rownames(tissuecombornal2fcdisplay) <- Reduce(union,list(rownames(bloodtfrnal2fcdisplay),rownames(adiposetfrnal2fcdisplay),rownames(muscletfrnal2fcdisplay)))
colnames(tissuecombornal2fcdisplay) <- c(colnames(muscletfrnal2fcdisplay),colnames(adiposetfrnal2fcdisplay),colnames(bloodtfrnal2fcdisplay))

tissuecombornapvaldisplay <- matrix(1L,nrow = length(Reduce(union,list(rownames(bloodtfrnapvaldisplay),rownames(adiposetfrnapvaldisplay),rownames(muscletfrnapvaldisplay)))),ncol = 22)
rownames(tissuecombornapvaldisplay) <- Reduce(union,list(rownames(bloodtfrnapvaldisplay),rownames(adiposetfrnapvaldisplay),rownames(muscletfrnapvaldisplay)))
colnames(tissuecombornapvaldisplay) <- c(colnames(muscletfrnapvaldisplay),colnames(adiposetfrnapvaldisplay),colnames(bloodtfrnapvaldisplay))

for(i in 1:dim(tissuecombornal2fcdisplay)[1]){
  ourtf <- rownames(tissuecombornal2fcdisplay)[i]
  if(ourtf %in% rownames(muscletfrnal2fcdisplay)){
    tissuecombornal2fcdisplay[i,c(1:6)] <- muscletfrnal2fcdisplay[ourtf,]
  }
  if(ourtf %in% rownames(adiposetfrnal2fcdisplay)){
    tissuecombornal2fcdisplay[i,c(7:12)] <- adiposetfrnal2fcdisplay[ourtf,]
  }
  if(ourtf %in% rownames(bloodtfrnal2fcdisplay)){
    tissuecombornal2fcdisplay[i,c(13:22)] <- bloodtfrnal2fcdisplay[ourtf,]
  }
  if(ourtf %in% rownames(muscletfrnapvaldisplay)){
    tissuecombornapvaldisplay[i,c(1:6)] <- muscletfrnapvaldisplay[ourtf,]
  }
  if(ourtf %in% rownames(adiposetfrnapvaldisplay)){
    tissuecombornapvaldisplay[i,c(7:12)] <- adiposetfrnapvaldisplay[ourtf,]
  }
  if(ourtf %in% rownames(bloodtfrnapvaldisplay)){
    tissuecombornapvaldisplay[i,c(13:22)] <- bloodtfrnapvaldisplay[ourtf,]
  }
}

tissuecombornal2fcdisplay <- tissuecombornal2fcdisplay[order(rownames(tissuecombornal2fcdisplay)),]
tissuecombornapvaldisplay <- tissuecombornapvaldisplay[order(rownames(tissuecombornapvaldisplay)),]

testout <- bubbleHeatmap(colorMat = tissuecombornal2fcdisplay[apply(tissuecombornapvaldisplay,1,min) < 0.00001,c(10,11,12,7,8,9,17,18,19,20,21,22,13,14,15,16,4,5,6,1,2,3)],sizeMat = -log10(tissuecombornapvaldisplay[apply(tissuecombornapvaldisplay,1,min) < 0.00001,c(10,11,12,7,8,9,17,18,19,20,21,22,13,14,15,16,4,5,6,1,2,3)]),legendTitles = c("Significance","log2FC"),colorSeq = c("blue","white","red"),colorLim = c(-1,1),sizeLim = c(0,5),)

# Supplemental Figure S6
png(file = "FigureS6_tissuecombo_TF_RNA_bubbleheatmap_102925.png",width = 8,height = 25,units = "in",res = 600)
grid.draw(testout)
dev.off()

png(filename = "combornaheatmap_topannotationtemplateforbubble_102925.png",width = 8,height = 8,units = "in",res = 600)
pheatmap(-log10(tissuecombornapvaldisplay[apply(tissuecombornapvaldisplay,1,min) < 0.00001,c(10,11,12,7,8,9,17,18,19,20,21,22,13,14,15,16,4,5,6,1,2,3)]),cluster_cols = F,breaks = seq(0,3,length.out = 101),color = colorpanel(101,"white","firebrick"),display_numbers = T,show_colnames = F,annotation_col = alltissue_tfenrich_qvalmatcolumndf[,c("Timepoint","Modality","Tissue")],annotation_colors = tf_ann_cols,cellwidth = 10)
dev.off()

#--------------------------------------------------------------------------------------
#Chris: This code is to replace the above bubbleheatmap code to keep consistent
#with all of the other plots and to avoid having to combine different figures
#in illustrator. - Nov 17th, 2025.
# Figure S6

keep_rows = apply(tissuecombornapvaldisplay, 1, min) < 0.00001
colorMat = tissuecombornal2fcdisplay[keep_rows, c(10,11,12,7,8,9,17,18,19,20,21,22,13,14,15,16,4,5,6,1,2,3)]
sizeMat  = -log10(tissuecombornapvaldisplay[keep_rows, c(10,11,12,7,8,9,17,18,19,20,21,22,13,14,15,16,4,5,6,1,2,3)])
# First half bc the new tissue combo is too long
n_total = nrow(colorMat)
half = n_total / 2
colorMat1 = colorMat[1:half, , drop = FALSE]
sizeMat1  = sizeMat[1:half, , drop = FALSE]

# Second half
colorMat2 = colorMat[(half + 1):n_total, , drop = FALSE]
sizeMat2  = sizeMat[(half + 1):n_total, , drop = FALSE]

#we copy this stuff from above, the only difference is we now want to add a tissue
#part of the legend, so we have to remake anno_df and add a color to anno_col
anno_df <- annotate_matrix(sizeMat) %>%
  dplyr::select(Tissue, Modality, Timepoint) %>%
  mutate(Timepoint = factor(Timepoint, levels = c("D20M",
                                                  "D40M",
                                                  "P10M",
                                                  "P15-45M",
                                                  "P3.5/4H",
                                                  "P24H"))) %>%
  distinct()
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
  )
)
ht1 = make_ht(colorMat1, sizeMat1)
pdf(file = "S6_firsthalf_tissue_combo_chris.pdf", height = 16, width = 8)
ComplexHeatmap::draw(
  ht1,
  annotation_legend_list = NULL,
  heatmap_legend_list = list(color_lg, size_lg),
  merge_legend = TRUE
)
dev.off()

ht2 = make_ht(colorMat2, sizeMat2)
pdf(file = "S6_secondhalf_tissue_combo_chris.pdf", height = 16, width = 8)
ComplexHeatmap::draw(
  ht2,
  annotation_legend_list = NULL,
  heatmap_legend_list = list(color_lg, size_lg),
  merge_legend = TRUE
)
dev.off()


#--------------------------------------------------------------------------------------




#####
# I next display the overlap between TFs that are enriched at given time points and showing responses
# or other time points. First up, I want to show the TFs that are generally enriched in a tissue and how
# that relates to the TFs that are responding at a gene, protein or phosphoprotein level
####


# Let's add some annotations to our boxes

alltissue_tfenrich_qvaldf$Timepoint <- factor(alltissue_tfenrich_qvaldf$Timepoint,levels = c("D20M","D40M","P10M","P15-45M","P3.5/4H","P24H"))

alltissue_tfenrich_qvaldf2 = distinct(alltissue_tfenrich_qvaldf, Timepoint, Tissue) %>%
  arrange(Tissue,Timepoint)
alltissue_tfenrich_qvaldf2$SigCount <- c(paste(sum(alltissue_tfenrich_qvalmat[,"Muscle_Endur 15/30/45 Min Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Muscle_Resist 15/30/45 Min Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Muscle_Endur 3.5/4 Hr Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Muscle_Resist 3.5/4 Hr Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Muscle_Endur 24 Hr Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Muscle_Resist 24 Hr Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Adipose_Endur 15/30/45 Min Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Adipose_Resist 15/30/45 Min Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Adipose_Endur 3.5/4 Hr Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Adipose_Resist 3.5/4 Hr Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Adipose_Endur 24 Hr Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Adipose_Resist 24 Hr Post"] < 0.05),sep = ""),
                                         sum(alltissue_tfenrich_qvalmat[,"Blood_Endur 20 Min During"] < 0.05),
                                         sum(alltissue_tfenrich_qvalmat[,"Blood_Endur 40 Min During"] < 0.05),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Blood_Endur 10 Min Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Blood_Resist 10 Min Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Blood_Endur 15/30/45 Min Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Blood_Resist 15/30/45 Min Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Blood_Endur 3.5/4 Hr Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Blood_Resist 3.5/4 Hr Post"] < 0.05),sep = ""),
                                         paste(sum(alltissue_tfenrich_qvalmat[,"Blood_Endur 24 Hr Post"] < 0.05),";",sum(alltissue_tfenrich_qvalmat[,"Blood_Resist 24 Hr Post"] < 0.05),sep = ""))

# Now let's just do the Sig counts without the qval distribution
alltissue_tfenrich_qvaldf2$EE.Sig <- gsub(";.*","",alltissue_tfenrich_qvaldf2$SigCount)
alltissue_tfenrich_qvaldf2$RE.Sig <- gsub(".*;","",alltissue_tfenrich_qvaldf2$SigCount)
alltissue_tfenrich_qvaldf2[7,"RE.Sig"] <- "0"
alltissue_tfenrich_qvaldf2[8,"RE.Sig"] <- "0"

alltissue_tfenrich_qvaldf3 <- rbind(alltissue_tfenrich_qvaldf2,alltissue_tfenrich_qvaldf2)
alltissue_tfenrich_qvaldf3[,"Sig"] <- c(alltissue_tfenrich_qvaldf2$EE.Sig,alltissue_tfenrich_qvaldf2$RE.Sig)
alltissue_tfenrich_qvaldf3[,"Modality"] <- c(rep("EE",dim(alltissue_tfenrich_qvaldf2)[1]),rep("RE",dim(alltissue_tfenrich_qvaldf2)[1]))
alltissue_tfenrich_qvaldf3$Sig <- as.numeric(alltissue_tfenrich_qvaldf3$Sig)

alltissue_tfenrich_qvaldf3$Timepoint <- factor(alltissue_tfenrich_qvaldf3$Timepoint,levels = c("D20M","D40M","P10M","P15-45M","P3.5/4H","P24H"))

alltissue_tfenrich_qvaldf3[19,"Sig"] <- NA
alltissue_tfenrich_qvaldf3[20,"Sig"] <- NA

# Figure 6B
png(file = "Figure6B_Sig Enriched TF Barplot_102925.png",width = 5,height = 7,units = "in",res = 600)
ggbarplot(alltissue_tfenrich_qvaldf3,x = "Timepoint",y = "Sig",fill = "Modality",color = "Modality",position = position_dodge(0.75)) + facet_wrap(facets = "Tissue",ncol = 1) + theme(axis.text = element_text(angle = 45,hjust = 1)) + scale_fill_manual(values = tf_ann_cols$Modality) + ylab("Significantly Enriched TF Count") + xlab("Time Point")
dev.off()

#####
# We next create a series of bar plots describing whether, in each tissue and for each ome, if TFs have
# significant responses to exercise and if they are enriched for DEGs in a given tissue
####

muscleenrichtflist <- rownames(alltissue_tfenrich_qvalmat)[apply(alltissue_tfenrich_qvalmat[,c(1:6)],1,min) < 0.05]
adiposeenrichtflist <- rownames(alltissue_tfenrich_qvalmat)[apply(alltissue_tfenrich_qvalmat[,c(7:12)],1,min) < 0.05]
bloodenrichtflist <- rownames(alltissue_tfenrich_qvalmat)[apply(alltissue_tfenrich_qvalmat[,c(13:22)],1,min) < 0.05]

muscleenrichtflist_ens <- human_feature_to_gene[human_feature_to_gene$assay %in% "transcript-rna-seq" & human_feature_to_gene$gene_symbol %in% tfproanno[muscleenrichtflist,"Gene.Name"],"feature_id"]
adiposeenrichtflist_ens <- human_feature_to_gene[human_feature_to_gene$assay %in% "transcript-rna-seq" & human_feature_to_gene$gene_symbol %in% tfproanno[adiposeenrichtflist,"Gene.Name"],"feature_id"]
bloodenrichtflist_ens <- human_feature_to_gene[human_feature_to_gene$assay %in% "transcript-rna-seq" & human_feature_to_gene$gene_symbol %in% tfproanno[bloodenrichtflist,"Gene.Name"],"feature_id"]

muscleenrichtflist_pro <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-pr" & human_feature_to_gene$gene_symbol %in% tfproanno[muscleenrichtflist,"Gene.Name"],"feature_id"]
adiposeenrichtflist_pro <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-pr" & human_feature_to_gene$gene_symbol %in% tfproanno[adiposeenrichtflist,"Gene.Name"],"feature_id"]
bloodenrichtflist_pro <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-pr" & human_feature_to_gene$gene_symbol %in% tfproanno[bloodenrichtflist,"Gene.Name"],"feature_id"]

muscleenrichtflist_proph <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph" & human_feature_to_gene$gene_symbol %in% tfproanno[muscleenrichtflist,"Gene.Name"],"feature_id"]
adiposeenrichtflist_proph <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph" & human_feature_to_gene$gene_symbol %in% tfproanno[adiposeenrichtflist,"Gene.Name"],"feature_id"]
bloodenrichtflist_proph <- human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph" & human_feature_to_gene$gene_symbol %in% tfproanno[bloodenrichtflist,"Gene.Name"],"feature_id"]


musclernasig <- unique(MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$adj_p_value < 0.05 &
                                                                           MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]

adiposernasig <- unique(MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$adj_p_value < 0.05 &
                                                                             MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]

bloodrnasig <- unique(MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA[MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$adj_p_value < 0.05 &
                                                                         MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUEndur.post_10_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUEndur.during_40_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_40_min + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_20_min + group_timepointADUControl.pre_exercise",
                                                                                                                                          "group_timepointADUResist.post_10_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]

muscleprosig <- unique(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$adj_p_value < 0.05 &
                                                                          MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                           "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                           "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                           "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                           "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                           "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]

adiposeprosig <- unique(MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$adj_p_value < 0.05 &
                                                                            MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                              "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                              "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                              "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                              "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                              "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]

bloodprosig <- unique(MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA[MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$adj_p_value < 0.05 &
                                                                        MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUEndur.post_10_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUEndur.during_40_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_40_min + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUEndur.during_20_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.during_20_min + group_timepointADUControl.pre_exercise",
                                                                                                                                        "group_timepointADUResist.post_10_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_10_min + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]


muscleprophsig <- unique(MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA[MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$adj_p_value < 0.05 &
                                                                            MotrpacHumanPreSuspensionData::MUSCLE_PROT_PH_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                             "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]

adiposeprophsig <- unique(MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA[MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$adj_p_value < 0.05 &
                                                                              MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PH_DA$contrast %in% c("group_timepointADUEndur.post_15_30_45_min - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUEndur.post_3.5_4_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUResist.post_15_30_45_min - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_15_30_45_min + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUResist.post_3.5_4_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_3.5_4_hr + group_timepointADUControl.pre_exercise",
                                                                                                                                                "group_timepointADUResist.post_24_hr - group_timepointADUResist.pre_exercise - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"),"feature_id"])[[1]]



muscleprosig_ens <- as.data.frame(HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$feature_id %in% muscleprosig,"ensembl_gene"])
muscleprosig_ens <- unique(muscleprosig_ens$ensembl_gene)

adiposeprosig_ens <- as.data.frame(HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$feature_id %in% adiposeprosig,"ensembl_gene"])
adiposeprosig_ens <- unique(adiposeprosig_ens$ensembl_gene)

bloodprosig_ens <- as.data.frame(HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$feature_id %in% bloodprosig,"ensembl_gene"])
bloodprosig_ens <- unique(bloodprosig_ens$ensembl_gene)

muscleprophsig_ens <- as.data.frame(HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$feature_id %in% muscleprophsig,"ensembl_gene"])
muscleprophsig_ens <- unique(muscleprophsig_ens$ensembl_gene)

adiposeprophsig_ens <- as.data.frame(HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$feature_id %in% adiposeprophsig,"ensembl_gene"])
adiposeprophsig_ens <- unique(adiposeprophsig_ens$ensembl_gene)


muscletfcompbarplot_stackeddf <- data.frame("Ome" = c(rep("Transcriptomics",4),rep("Proteomics (MS)",4),rep("Phosphoproteomics",4)),
                                            "Classification" = c(rep(c("Not Sig, Not Enriched","Not Sig, Enriched","Sig, Not Enriched","Sig, Enriched"),3)),
                                            "Count" = rep(0,12))

muscletfcompbarplot_stackeddf[1,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])[1] - dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",musclernasig),])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% gsub("\\..*","",musclernasig)),]),rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])))
muscletfcompbarplot_stackeddf[2,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% gsub("\\..*","",musclernasig)),]),rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])))
muscletfcompbarplot_stackeddf[3,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",musclernasig),])[1] - length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",musclernasig),]),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])))
muscletfcompbarplot_stackeddf[4,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",musclernasig),]),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_DA$feature_id),])))

muscletfcompbarplot_stackeddf[5,"Count"] <- dim(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])[1] - dim(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",muscleprosig),])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$prot_pr_id %in% gsub("\\..*","",muscleprosig)),]),rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])))
muscletfcompbarplot_stackeddf[6,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$prot_pr_id %in% gsub("\\..*","",muscleprosig)),]),rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])))
muscletfcompbarplot_stackeddf[7,"Count"] <- dim(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",muscleprosig),])[1] - length(intersect(rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",muscleprosig),]),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])))
muscletfcompbarplot_stackeddf[8,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",muscleprosig),]),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::MUSCLE_PROT_PR_DA$feature_id),])))

muscletfcompbarplot_stackeddf[9,"Count"] <- dim(tfproanno[tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])[1] - dim(tfproanno[tfproanno$Ensembl %in% muscleprophsig_ens,])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% muscleprophsig_ens),]),rownames(tfproanno[tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])))
muscletfcompbarplot_stackeddf[10,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% muscleprophsig_ens),]),rownames(tfproanno[tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])),rownames(tfproanno[rownames(tfproanno) %in% muscleenrichtflist & tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])))
muscletfcompbarplot_stackeddf[11,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% muscleprophsig_ens,])[1] - length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% muscleprophsig_ens,]),muscleenrichtflist))
muscletfcompbarplot_stackeddf[12,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% muscleprophsig_ens,]),muscleenrichtflist))

muscletfcompbarplot_stackeddf$Classification <- factor(muscletfcompbarplot_stackeddf$Classification,levels = c("Not Sig, Not Enriched","Not Sig, Enriched","Sig, Not Enriched","Sig, Enriched"))
muscletfcompbarplot_stackeddf$Ome <- factor(muscletfcompbarplot_stackeddf$Ome,levels = c("Transcriptomics","Proteomics (MS)","Phosphoproteomics"))

# FIGURE 6G.III
png(file = "Figure6GIII_Muscle Transcription Factor Activity Barplot_102925.png",width = 6, height = 6,units = "in",res = 600)
ggplot(muscletfcompbarplot_stackeddf, aes(fill=Classification, y=Count, x=Ome)) +
  geom_bar(position="stack", stat="identity") + theme_classic() + scale_fill_manual(values = c("#97c976","#dba24c","#fd717a","#ce69cd")) + ggtitle("Muscle Transcription Factor Activity")
dev.off()


adiposetfcompbarplot_stackeddf <- data.frame("Ome" = c(rep("Transcriptomics",4),rep("Proteomics (MS)",4),rep("Phosphoproteomics",4)),
                                             "Classification" = c(rep(c("Not Sig, Not Enriched","Not Sig, Enriched","Sig, Not Enriched","Sig, Enriched"),3)),
                                             "Count" = rep(0,12))

adiposetfcompbarplot_stackeddf[1,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])[1] - dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",adiposernasig),])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% gsub("\\..*","",adiposernasig)),]),rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])))
adiposetfcompbarplot_stackeddf[2,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% gsub("\\..*","",adiposernasig)),]),rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])))
adiposetfcompbarplot_stackeddf[3,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",adiposernasig),])[1] - length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",adiposernasig),]),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])))
adiposetfcompbarplot_stackeddf[4,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",adiposernasig),]),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_DA$feature_id),])))

adiposetfcompbarplot_stackeddf[5,"Count"] <- dim(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])[1] - dim(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",adiposeprosig),])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$prot_pr_id %in% gsub("\\..*","",adiposeprosig)),]),rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])))
adiposetfcompbarplot_stackeddf[6,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$prot_pr_id %in% gsub("\\..*","",adiposeprosig)),]),rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])))
adiposetfcompbarplot_stackeddf[7,"Count"] <- dim(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",adiposeprosig),])[1] - length(intersect(rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",adiposeprosig),]),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])))
adiposetfcompbarplot_stackeddf[8,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$prot_pr_id %in% gsub("\\..*","",adiposeprosig),]),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$prot_pr_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::ADIPOSE_PROT_PR_DA$feature_id),])))

adiposetfcompbarplot_stackeddf[9,"Count"] <- dim(tfproanno[tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])[1] - dim(tfproanno[tfproanno$Ensembl %in% adiposeprophsig_ens,])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% adiposeprophsig_ens),]),rownames(tfproanno[tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])))
adiposetfcompbarplot_stackeddf[10,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% adiposeprophsig_ens),]),rownames(tfproanno[tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])),rownames(tfproanno[rownames(tfproanno) %in% adiposeenrichtflist & tfproanno$Gene.Name %in% human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph","gene_symbol"],])))
adiposetfcompbarplot_stackeddf[11,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% adiposeprophsig_ens,])[1] - length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% adiposeprophsig_ens,]),adiposeenrichtflist))
adiposetfcompbarplot_stackeddf[12,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% adiposeprophsig_ens,]),adiposeenrichtflist))

adiposetfcompbarplot_stackeddf$Classification <- factor(adiposetfcompbarplot_stackeddf$Classification,levels = c("Not Sig, Not Enriched","Not Sig, Enriched","Sig, Not Enriched","Sig, Enriched"))
adiposetfcompbarplot_stackeddf$Ome <- factor(adiposetfcompbarplot_stackeddf$Ome,levels = c("Transcriptomics","Proteomics (MS)","Phosphoproteomics"))

# FIGURE 6G.I
png(file = "Figure6GI_Adipose Transcription Factor Activity Barplot_102925.png",width = 6, height = 6,units = "in",res = 600)
ggplot(adiposetfcompbarplot_stackeddf, aes(fill=Classification, y=Count, x=Ome)) +
  geom_bar(position="stack", stat="identity") + theme_classic() + scale_fill_manual(values = c("#97c976","#dba24c","#fd717a","#ce69cd")) + ggtitle("Adipose Transcription Factor Activity")
dev.off()


bloodtfcompbarplot_stackeddf <- data.frame("Ome" = c(rep("Transcriptomics",4),rep("Proteomics (Olink)",4)),
                                           "Classification" = c(rep(c("Not Sig, Not Enriched","Not Sig, Enriched","Sig, Not Enriched","Sig, Enriched"),2)),
                                           "Count" = rep(0,8))

bloodtfcompbarplot_stackeddf[1,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])[1] - dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",bloodrnasig),])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% gsub("\\..*","",bloodrnasig)),]),rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])))
bloodtfcompbarplot_stackeddf[2,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$Ensembl %in% gsub("\\..*","",bloodrnasig)),]),rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])))
bloodtfcompbarplot_stackeddf[3,"Count"] <- dim(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",bloodrnasig),])[1] - length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",bloodrnasig),]),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])))
bloodtfcompbarplot_stackeddf[4,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$Ensembl %in% gsub("\\..*","",bloodrnasig),]),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$Ensembl %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_TRNSCRPT_DA$feature_id),])))

bloodtfcompbarplot_stackeddf[5,"Count"] <- dim(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])[1] - dim(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",bloodprosig),])[1] - length(intersect(intersect(rownames(tfproanno[!(tfproanno$prot_ol_id %in% gsub("\\..*","",bloodprosig)),]),rownames(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])))
bloodtfcompbarplot_stackeddf[6,"Count"] <- length(intersect(intersect(rownames(tfproanno[!(tfproanno$prot_ol_id %in% gsub("\\..*","",bloodprosig)),]),rownames(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])))
bloodtfcompbarplot_stackeddf[7,"Count"] <- dim(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",bloodprosig),])[1] - length(intersect(rownames(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",bloodprosig),]),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])))
bloodtfcompbarplot_stackeddf[8,"Count"] <- length(intersect(rownames(tfproanno[tfproanno$prot_ol_id %in% gsub("\\..*","",bloodprosig),]),rownames(tfproanno[rownames(tfproanno) %in% bloodenrichtflist & tfproanno$prot_ol_id %in% gsub("\\..*","",MotrpacHumanPreSuspensionData::BLOOD_PROT_OL_DA$feature_id),])))

bloodtfcompbarplot_stackeddf$Classification <- factor(bloodtfcompbarplot_stackeddf$Classification,levels = c("Not Sig, Not Enriched","Not Sig, Enriched","Sig, Not Enriched","Sig, Enriched"))
bloodtfcompbarplot_stackeddf$Ome <- factor(bloodtfcompbarplot_stackeddf$Ome,levels = c("Transcriptomics","Proteomics (Olink)"))

# FIGURE 6G.II
png(file = "Figure6GII_Blood Transcription Factor Activity Barplot_102925.png",width = 5, height = 6,units = "in",res = 600)
ggplot(bloodtfcompbarplot_stackeddf, aes(fill=Classification, y=Count, x=Ome)) +
  geom_bar(position="stack", stat="identity") + theme_classic() + scale_fill_manual(values = c("#97c976","#dba24c","#fd717a","#ce69cd")) + ggtitle("Blood Transcription Factor Activity")
dev.off()

#####
# Now we generate Supplemental Table S6 describing the TF enrichments
# and exercise responses across all three omes.
####

tfenrichexpresstable <- matrix(1L,nrow = length(tflist),ncol = 70)
rownames(tfenrichexpresstable) <- tflist
colnames(tfenrichexpresstable) <- c("Muscle Enrich EE P15-45M","Muscle Enrich EE P3.5/4H","Muscle Enrich EE P24H",
                                    "Muscle Enrich RE P15-45M","Muscle Enrich RE P3.5/4H","Muscle Enrich RE P24H",
                                    "Adipose Enrich EE P15-45M","Adipose Enrich EE P3.5/4H","Adipose Enrich EE P24H",
                                    "Adipose Enrich RE P15-45M","Adipose Enrich RE P3.5/4H","Adipose Enrich RE P24H",
                                    "Blood Enrich EE D20M","Blood Enrich EE D40M","Blood Enrich EE P10M",
                                    "Blood Enrich EE P15-45M","Blood Enrich EE P3.5/4H","Blood Enrich EE P24H",
                                    "Blood Enrich RE P10M","Blood Enrich RE P15-45M","Blood Enrich RE P3.5/4H","Blood Enrich RE P24H",
                                    "Muscle RNA Expression EE P15-45M","Muscle RNA Expression EE P3.5/4H","Muscle RNA Expression EE P24H",
                                    "Muscle RNA Expression RE P15-45M","Muscle RNA Expression RE P3.5/4H","Muscle RNA Expression RE P24H",
                                    "Adipose RNA Expression EE P15-45M","Adipose RNA Expression EE P3.5/4H","Adipose RNA Expression EE P24H",
                                    "Adipose RNA Expression RE P15-45M","Adipose RNA Expression RE P3.5/4H","Adipose RNA Expression RE P24H",
                                    "Blood RNA Expression EE D20M","Blood RNA Expression EE D40M","Blood RNA Expression EE P10M",
                                    "Blood RNA Expression EE P15-45M","Blood RNA Expression EE P3.5/4H","Blood RNA Expression EE P24H",
                                    "Blood RNA Expression RE P10M","Blood RNA Expression RE P15-45M","Blood RNA Expression RE P3.5/4H","Blood RNA Expression RE P24H",
                                    "Muscle Protein Abundance EE P15-45M","Muscle Protein Abundance EE P3.5/4H","Muscle Protein Abundance EE P24H",
                                    "Muscle Protein Abundance RE P15-45M","Muscle Protein Abundance RE P3.5/4H","Muscle Protein Abundance RE P24H",
                                    "Adipose Protein Abundance EE P3.5/4H","Adipose Protein Abundance RE P3.5/4H",
                                    "Blood Protein Abundance EE D20M","Blood Protein Abundance EE D40M","Blood Protein Abundance EE P10M",
                                    "Blood Protein Abundance EE P15-45M","Blood Protein Abundance EE P3.5/4H","Blood Protein Abundance EE P24H",
                                    "Blood Protein Abundance RE P10M","Blood Protein Abundance RE P15-45M","Blood Protein Abundance RE P3.5/4H","Blood Protein Abundance RE P24H",
                                    "Muscle Protein Phosphorylation EE P15-45M","Muscle Protein Phosphorylation EE P3.5/4H","Muscle Protein Phosphorylation EE P24H",
                                    "Muscle Protein Phosphorylation RE P15-45M","Muscle Protein Phosphorylation RE P3.5/4H","Muscle Protein Phosphorylation RE P24H",
                                    "Adipose Protein Phosphorylation EE P3.5/4H","Adipose Protein Phosphorylation RE P3.5/4H")

tfenrichexpresstable[tflist,c(1:6)] <- musclerna_tfenrich_qvalmat[tflist,]
tfenrichexpresstable[tflist,c(7:12)] <- adiposerna_tfenrich_qvalmat[tflist,]
tfenrichexpresstable[tflist,c(13:22)] <- bloodrna_tfenrich_qvalmat[tflist,]

for(i in 1:length(tflist)){

  if(i%%25 == 0){
    print(i)
  }

  ourtf <- tflist[i]
  if(tfproanno[ourtf,"Ensembl"] %in% gsub("\\..*","",rownames(muscletfrnapval))){
    tfenrichexpresstable[i,"Muscle RNA Expression RE P15-45M"] <- muscletfrnapval[rownames(muscletfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(muscletfrnapval))],"Muscle_Resist 15/30/45 Min Post"]
    tfenrichexpresstable[i,"Muscle RNA Expression RE P3.5/4H"] <- muscletfrnapval[rownames(muscletfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(muscletfrnapval))],"Muscle_Endur 3.5/4 Hr Post"]
    tfenrichexpresstable[i,"Muscle RNA Expression RE P24H"] <- muscletfrnapval[rownames(muscletfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(muscletfrnapval))],"Muscle_Resist 24 Hr Post"]
    tfenrichexpresstable[i,"Muscle RNA Expression EE P15-45M"] <- muscletfrnapval[rownames(muscletfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(muscletfrnapval))],"Muscle_Endur 15/30/45 Min Post"]
    tfenrichexpresstable[i,"Muscle RNA Expression EE P3.5/4H"] <- muscletfrnapval[rownames(muscletfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(muscletfrnapval))],"Muscle_Endur 3.5/4 Hr Post"]
    tfenrichexpresstable[i,"Muscle RNA Expression EE P24H"] <- muscletfrnapval[rownames(muscletfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(muscletfrnapval))],"Muscle_Endur 24 Hr Post"]
  }
  if(tfproanno[ourtf,"Ensembl"] %in% gsub("\\..*","",rownames(adiposetfrnapval))){
    tfenrichexpresstable[i,"Adipose RNA Expression RE P15-45M"] <- adiposetfrnapval[rownames(adiposetfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(adiposetfrnapval))],"Adipose_Resist 15/30/45 Min Post"]
    tfenrichexpresstable[i,"Adipose RNA Expression RE P3.5/4H"] <- adiposetfrnapval[rownames(adiposetfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(adiposetfrnapval))],"Adipose_Endur 3.5/4 Hr Post"]
    tfenrichexpresstable[i,"Adipose RNA Expression RE P24H"] <- adiposetfrnapval[rownames(adiposetfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(adiposetfrnapval))],"Adipose_Resist 24 Hr Post"]
    tfenrichexpresstable[i,"Adipose RNA Expression EE P15-45M"] <- adiposetfrnapval[rownames(adiposetfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(adiposetfrnapval))],"Adipose_Endur 15/30/45 Min Post"]
    tfenrichexpresstable[i,"Adipose RNA Expression EE P3.5/4H"] <- adiposetfrnapval[rownames(adiposetfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(adiposetfrnapval))],"Adipose_Endur 3.5/4 Hr Post"]
    tfenrichexpresstable[i,"Adipose RNA Expression EE P24H"] <- adiposetfrnapval[rownames(adiposetfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(adiposetfrnapval))],"Adipose_Endur 24 Hr Post"]
  }
  if(tfproanno[ourtf,"Ensembl"] %in% gsub("\\..*","",rownames(bloodtfrnapval))){
    tfenrichexpresstable[i,"Blood RNA Expression RE P10M"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Resist 10 Min Post"]
    tfenrichexpresstable[i,"Blood RNA Expression RE P15-45M"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Resist 15/30/45 Min Post"]
    tfenrichexpresstable[i,"Blood RNA Expression RE P3.5/4H"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 3.5/4 Hr Post"]
    tfenrichexpresstable[i,"Blood RNA Expression RE P24H"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Resist 24 Hr Post"]
    tfenrichexpresstable[i,"Blood RNA Expression EE D20M"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 20 Min During"]
    tfenrichexpresstable[i,"Blood RNA Expression EE D40M"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 40 Min During"]
    tfenrichexpresstable[i,"Blood RNA Expression EE P10M"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 10 Min Post"]
    tfenrichexpresstable[i,"Blood RNA Expression EE P15-45M"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 15/30/45 Min Post"]
    tfenrichexpresstable[i,"Blood RNA Expression EE P3.5/4H"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 3.5/4 Hr Post"]
    tfenrichexpresstable[i,"Blood RNA Expression EE P24H"] <- bloodtfrnapval[rownames(bloodtfrnapval)[grep(tfproanno[ourtf,"Ensembl"],rownames(bloodtfrnapval))],"Blood_Endur 24 Hr Post"]
  }


  if(tfproanno[ourtf,"prot_pr_id"] %in% gsub("\\..*","",rownames(muscletfprotprpval))){
    tfenrichexpresstable[i,"Muscle Protein Abundance RE P15-45M"] <- min(muscletfprotprpval[rownames(muscletfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(muscletfprotprpval))],"Muscle_Resist 15/30/45 Min Post"])
    tfenrichexpresstable[i,"Muscle Protein Abundance RE P3.5/4H"] <- min(muscletfprotprpval[rownames(muscletfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(muscletfprotprpval))],"Muscle_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Muscle Protein Abundance RE P24H"] <- min(muscletfprotprpval[rownames(muscletfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(muscletfprotprpval))],"Muscle_Resist 24 Hr Post"])
    tfenrichexpresstable[i,"Muscle Protein Abundance EE P15-45M"] <- min(muscletfprotprpval[rownames(muscletfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(muscletfprotprpval))],"Muscle_Endur 15/30/45 Min Post"])
    tfenrichexpresstable[i,"Muscle Protein Abundance EE P3.5/4H"] <- min(muscletfprotprpval[rownames(muscletfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(muscletfprotprpval))],"Muscle_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Muscle Protein Abundance EE P24H"] <- min(muscletfprotprpval[rownames(muscletfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(muscletfprotprpval))],"Muscle_Endur 24 Hr Post"])
  }
  if(tfproanno[ourtf,"prot_pr_id"] %in% gsub("\\..*","",rownames(adiposetfprotprpval))){
    tfenrichexpresstable[i,"Adipose Protein Abundance RE P3.5/4H"] <- min(adiposetfprotprpval[rownames(adiposetfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(adiposetfprotprpval))],"Adipose_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Adipose Protein Abundance EE P3.5/4H"] <- min(adiposetfprotprpval[rownames(adiposetfprotprpval)[grep(tfproanno[ourtf,"prot_pr_id"],rownames(adiposetfprotprpval))],"Adipose_Endur 3.5/4 Hr Post"])
  }
  if(tfproanno[ourtf,"prot_ol_id"] %in% gsub("\\..*","",rownames(bloodtfprotolpval))){
    tfenrichexpresstable[i,"Blood Protein Abundance RE P10M"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Resist 10 Min Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance RE P15-45M"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Resist 15/30/45 Min Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance RE P3.5/4H"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance RE P24H"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Resist 24 Hr Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance EE D20M"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 20 Min During"])
    tfenrichexpresstable[i,"Blood Protein Abundance EE D40M"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 40 Min During"])
    tfenrichexpresstable[i,"Blood Protein Abundance EE P10M"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 10 Min Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance EE P15-45M"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 15/30/45 Min Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance EE P3.5/4H"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Blood Protein Abundance EE P24H"] <- min(bloodtfprotolpval[rownames(bloodtfprotolpval)[grep(tfproanno[ourtf,"prot_ol_id"],rownames(bloodtfprotolpval))],"Blood_Endur 24 Hr Post"])
  }


  ourtfgene <- tfproanno[ourtf,"Gene.Name"]
  ourtfphospho <- as.character(human_feature_to_gene[human_feature_to_gene$assay %in% "prot-ph" & human_feature_to_gene$gene_symbol %in% ourtfgene,"feature_id"])

  if(length(intersect(ourtfphospho,rownames(muscletfprotphpval))) > 0){
    tfenrichexpresstable[i,"Muscle Protein Phosphorylation RE P15-45M"] <- min(muscletfprotphpval[intersect(ourtfphospho,rownames(muscletfprotphpval)),"Muscle_Resist 15/30/45 Min Post"])
    tfenrichexpresstable[i,"Muscle Protein Phosphorylation RE P3.5/4H"] <- min(muscletfprotphpval[intersect(ourtfphospho,rownames(muscletfprotphpval)),"Muscle_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Muscle Protein Phosphorylation RE P24H"] <- min(muscletfprotphpval[intersect(ourtfphospho,rownames(muscletfprotphpval)),"Muscle_Resist 24 Hr Post"])
    tfenrichexpresstable[i,"Muscle Protein Phosphorylation EE P15-45M"] <- min(muscletfprotphpval[intersect(ourtfphospho,rownames(muscletfprotphpval)),"Muscle_Endur 15/30/45 Min Post"])
    tfenrichexpresstable[i,"Muscle Protein Phosphorylation EE P3.5/4H"] <- min(muscletfprotphpval[intersect(ourtfphospho,rownames(muscletfprotphpval)),"Muscle_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Muscle Protein Phosphorylation EE P24H"] <- min(muscletfprotphpval[intersect(ourtfphospho,rownames(muscletfprotphpval)),"Muscle_Endur 24 Hr Post"])
  }
  if(length(intersect(ourtfphospho,rownames(adiposetfprotphpval))) > 0){
    tfenrichexpresstable[i,"Adipose Protein Phosphorylation RE P3.5/4H"] <- min(adiposetfprotphpval[intersect(ourtfphospho,rownames(adiposetfprotphpval)),"Adipose_Endur 3.5/4 Hr Post"])
    tfenrichexpresstable[i,"Adipose Protein Phosphorylation EE P3.5/4H"] <- min(adiposetfprotphpval[intersect(ourtfphospho,rownames(adiposetfprotphpval)),"Adipose_Endur 3.5/4 Hr Post"])
  }
}

write.csv(tfenrichexpresstable,file = "Pre-covid Landscape Supplemental Table S6.csv",row.names = T)

#####
# Lastly, we generate scatter plots comparing the EE and RE responses to acute exercise by tissue for both
# Transcriptomics and Phosphoproteomics
####

musclevsadiposeTFprophdf <- data.frame("TF" = c(rownames(muscletfprotphl2fcdisplay),
                                                rownames(adiposetfprotphl2fcdisplay)),
                                       "RE.L2FC" = c(muscletfprotphl2fcdisplay[,"Muscle_Resist 3.5/4 Hr Post"],
                                                     adiposetfprotphl2fcdisplay[,"Adipose_Resist 3.5/4 Hr Post"]),
                                       "EE.L2FC" = c(muscletfprotphl2fcdisplay[,"Muscle_Endur 3.5/4 Hr Post"],
                                                     adiposetfprotphl2fcdisplay[,"Adipose_Endur 3.5/4 Hr Post"]),
                                       "Tissue" = c(rep("Muscle",length(rownames(muscletfprotphl2fcdisplay))),
                                                    rep("Adipose",length(rownames(adiposetfprotphl2fcdisplay)))))

musclevsadiposeTFprophdf$TF <- gsub("s|t|y","",musclevsadiposeTFprophdf$TF)

# Figure 6F
png(file = "Figure6F_Muscle vs Adipose Phospho EE vs RE 4 Hour Response Scatter Plot_102925.png",width = 5,height = 5,units = "in",res = 600)
ggplot(musclevsadiposeTFprophdf, aes(x = RE.L2FC, y = EE.L2FC, color = Tissue,label = TF)) + geom_abline(1,intercept = 0,colour = "grey50",linetype = "dashed") + geom_hline(yintercept = 0,colour = "grey40") + geom_vline(xintercept = 0,colour = "grey40") + geom_point() +
  geom_text_repel(max.overlaps = 15) +
  geom_labelsmooth(aes(label = Tissue), fill = "white",
                   method = "lm", formula = y ~ x,
                   size = 5, linewidth = 1, boxlinewidth = 0.8) +
  theme_bw() + guides(color = 'none') + xlim(-0.8,1.2) + ylim(-0.8,1.2) + ggtitle("Comparing EE vs RE in Phosphoproteomics P3.5/4H") + scale_color_manual(values = list("Muscle" = "darkblue","Adipose" = "goldenrod","Blood" = "firebrick")) + theme_classic()
dev.off()



musclevsbloodvsadiposeTFrnadf <- data.frame("TF" = c(rownames(muscletfrnal2fcdisplay),
                                                     rownames(adiposetfrnal2fcdisplay),
                                                     rownames(bloodtfrnal2fcdisplay)),
                                            "RE.L2FC" = c(muscletfrnal2fcdisplay[,"Muscle_Resist 3.5/4 Hr Post"],
                                                          adiposetfrnal2fcdisplay[,"Adipose_Resist 3.5/4 Hr Post"],
                                                          bloodtfrnal2fcdisplay[,"Blood_Resist 3.5/4 Hr Post"]),
                                            "EE.L2FC" = c(muscletfrnal2fcdisplay[,"Muscle_Endur 3.5/4 Hr Post"],
                                                          adiposetfrnal2fcdisplay[,"Adipose_Endur 3.5/4 Hr Post"],
                                                          bloodtfrnal2fcdisplay[,"Blood_Endur 3.5/4 Hr Post"]),
                                            "Tissue" = c(rep("Muscle",length(rownames(muscletfrnal2fcdisplay))),
                                                         rep("Adipose",length(rownames(adiposetfrnal2fcdisplay))),
                                                         rep("Blood",length(rownames(bloodtfrnal2fcdisplay)))))

# Figure 6E
png(file = "Figure6E_EE vs RE 4 Hour Response Scatter Plot_102925.png",width = 5,height = 5,units = "in",res = 600)
ggplot(musclevsbloodvsadiposeTFrnadf, aes(x = RE.L2FC, y = EE.L2FC, color = Tissue,label = TF)) + geom_abline(1,intercept = 0,colour = "grey50",linetype = "dashed") + geom_hline(yintercept = 0,linetype = "dashed",colour = "grey40") + geom_vline(xintercept = 0,linetype = "dashed",colour = "grey40") + geom_point() +
  geom_text_repel(max.overlaps = 15) +
  geom_labelsmooth(aes(label = Tissue), fill = "white",
                   method = "lm", formula = y ~ x,
                   size = 5, linewidth = 1, boxlinewidth = 0.8) +
  theme_bw() + guides(color = 'none') + xlim(-2,5.75) + ylim(-2,5.75) + ggtitle("Comparing EE vs RE in Transcriptomics P3.5/4H") + scale_color_manual(values = list("Muscle" = "darkblue","Adipose" = "goldenrod","Blood" = "firebrick")) + theme_classic()
dev.off()

# In case you want to save progress
# save.image("precovid_homer_analysis_102925.RData")
