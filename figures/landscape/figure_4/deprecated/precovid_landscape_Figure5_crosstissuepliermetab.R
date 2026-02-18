library(MotrpacHumanPreSuspension)
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
library(circlize)
library(umap)
library(bubbleHeatmap)

# setwd("C:/Users/gsmit/Documents/Mount Sinai/Sealfon Laboratory/MoTrPAC/PreCovid Human")

# Loading the refmet metabolite classifications for use in PLIER's prior knowledge
refmetmetabolites <- Reduce(union,MOLECULAR_SIGNATURES$REFMET)

refmetmat <- matrix(0L,nrow = length(refmetmetabolites),ncol = length(names(MOLECULAR_SIGNATURES$REFMET)))
rownames(refmetmat) <- refmetmetabolites
colnames(refmetmat) <- names(MOLECULAR_SIGNATURES$REFMET)
for(i in 1:length(names(MOLECULAR_SIGNATURES$REFMET))){
  refmetmat[MOLECULAR_SIGNATURES$REFMET[i][[1]],i] <- 1
}

# Collecting all of the qcnorm matrices from the metabolomics data in each tissue
musclemetab_t_amine_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$qc_norm
musclemetab_t_nuc_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_NUC_QC$qc_norm
musclemetab_t_oxylipneg_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_OXYLIPNEG_QC$qc_norm
musclemetab_t_tca_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_TCA_QC$qc_norm
musclemetab_u_hilicpos_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_HILICPOS_QC$qc_norm
musclemetab_u_ionpneg_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_IONPNEG_QC$qc_norm
musclemetab_u_rpneg_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_RPNEG_QC$qc_norm
musclemetab_u_rppos_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_RPPOS_QC$qc_norm
musclemetab_u_lrpneg_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_LRPNEG_QC$qc_norm
musclemetab_u_lrppos_data <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_LRPPOS_QC$qc_norm

bloodmetab_t_amine_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$qc_norm
bloodmetab_t_oxylipneg_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_T_OXYLIPNEG_QC$qc_norm
bloodmetab_t_tca_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_T_TCA_QC$qc_norm
bloodmetab_t_conv_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_T_CONV_QC$qc_norm
bloodmetab_t_imm_crt_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_T_IMM_CRT_QC$qc_norm
bloodmetab_u_hilicpos_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_U_HILICPOS_QC$qc_norm
bloodmetab_u_ionpneg_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_U_IONPNEG_QC$qc_norm
bloodmetab_u_rpneg_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_U_RPNEG_QC$qc_norm
bloodmetab_u_rppos_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_U_RPPOS_QC$qc_norm
bloodmetab_u_lrpneg_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_U_LRPNEG_QC$qc_norm
bloodmetab_u_lrppos_data <- MotrpacHumanPreSuspensionData::BLOOD_METAB_U_LRPPOS_QC$qc_norm

adiposemetab_t_amine_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$qc_norm
adiposemetab_t_nuc_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_NUC_QC$qc_norm
adiposemetab_t_oxylipneg_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_OXYLIPNEG_QC$qc_norm
adiposemetab_t_tca_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_TCA_QC$qc_norm
adiposemetab_t_acoa_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_ACOA_QC$qc_norm
adiposemetab_t_ka_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_KA_QC$qc_norm
adiposemetab_u_hilicpos_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_HILICPOS_QC$qc_norm
adiposemetab_u_ionpneg_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_IONPNEG_QC$qc_norm
adiposemetab_u_rpneg_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_RPNEG_QC$qc_norm
adiposemetab_u_rppos_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_RPPOS_QC$qc_norm
adiposemetab_u_lrpneg_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_LRPNEG_QC$qc_norm
adiposemetab_u_lrppos_data <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_LRPPOS_QC$qc_norm

# We identify the samples that overlap across each metabolite designation per tissue so we can integrate all matrices into a
colnames(musclemetab_t_amine_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_t_oxylipneg_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_OXYLIPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_OXYLIPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_t_tca_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_TCA_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_TCA_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_u_hilicpos_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_HILICPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_HILICPOS_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_u_ionpneg_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_IONPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_IONPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_u_rpneg_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_RPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_RPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_u_rppos_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_RPPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_RPPOS_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_u_lrpneg_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_LRPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_LRPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(musclemetab_u_lrppos_data) <- paste(MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_LRPPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::MUSCLE_METAB_U_LRPPOS_QC$sample_metadata$Timepoint,sep = "")

musclesamplelist <- Reduce(intersect,list(colnames(musclemetab_t_amine_data),
                                          colnames(musclemetab_t_oxylipneg_data),
                                          colnames(musclemetab_t_tca_data),
                                          colnames(musclemetab_u_hilicpos_data),
                                          colnames(musclemetab_u_ionpneg_data),
                                          colnames(musclemetab_u_rpneg_data),
                                          colnames(musclemetab_u_rppos_data),
                                          colnames(musclemetab_u_lrpneg_data),
                                          colnames(musclemetab_u_lrppos_data)))

colnames(adiposemetab_t_amine_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_t_oxylipneg_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_OXYLIPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_OXYLIPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_t_tca_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_TCA_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_TCA_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_u_hilicpos_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_HILICPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_HILICPOS_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_u_ionpneg_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_IONPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_IONPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_u_rpneg_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_RPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_RPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_u_rppos_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_RPPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_RPPOS_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_u_lrpneg_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_LRPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_LRPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(adiposemetab_u_lrppos_data) <- paste(MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_LRPPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::ADIPOSE_METAB_U_LRPPOS_QC$sample_metadata$Timepoint,sep = "")

adiposesamplelist <- Reduce(intersect,list(colnames(adiposemetab_t_amine_data),
                                           colnames(adiposemetab_t_oxylipneg_data),
                                           colnames(adiposemetab_t_tca_data),
                                           colnames(adiposemetab_u_hilicpos_data),
                                           colnames(adiposemetab_u_ionpneg_data),
                                           colnames(adiposemetab_u_rpneg_data),
                                           colnames(adiposemetab_u_rppos_data),
                                           colnames(adiposemetab_u_lrpneg_data),
                                           colnames(adiposemetab_u_lrppos_data)))

colnames(bloodmetab_t_amine_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_t_oxylipneg_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_T_OXYLIPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_T_OXYLIPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_t_tca_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_T_TCA_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_T_TCA_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_u_hilicpos_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_U_HILICPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_U_HILICPOS_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_u_ionpneg_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_U_IONPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_U_IONPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_u_rpneg_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_U_RPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_U_RPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_u_rppos_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_U_RPPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_U_RPPOS_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_u_lrpneg_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_U_LRPNEG_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_U_LRPNEG_QC$sample_metadata$Timepoint,sep = "")
colnames(bloodmetab_u_lrppos_data) <- paste(MotrpacHumanPreSuspensionData::BLOOD_METAB_U_LRPPOS_QC$sample_metadata$pid,"_",MotrpacHumanPreSuspensionData::BLOOD_METAB_U_LRPPOS_QC$sample_metadata$Timepoint,sep = "")

bloodsamplelist <- Reduce(intersect,list(colnames(bloodmetab_t_amine_data),
                                         colnames(bloodmetab_t_oxylipneg_data),
                                         colnames(bloodmetab_t_tca_data),
                                         colnames(bloodmetab_u_hilicpos_data),
                                         colnames(bloodmetab_u_ionpneg_data),
                                         colnames(bloodmetab_u_rpneg_data),
                                         colnames(bloodmetab_u_rppos_data),
                                         colnames(bloodmetab_u_lrpneg_data),
                                         colnames(bloodmetab_u_lrppos_data)))

# Now we merge each qcnorm matrix together into a single matrix per tissue
muscle_combo_rownames <- c(rownames(musclemetab_t_amine_data),
                           rownames(musclemetab_t_nuc_data),
                           rownames(musclemetab_t_oxylipneg_data),
                           rownames(musclemetab_t_tca_data),
                           rownames(musclemetab_u_hilicpos_data),
                           rownames(musclemetab_u_ionpneg_data),
                           rownames(musclemetab_u_rpneg_data),
                           rownames(musclemetab_u_rppos_data),
                           rownames(musclemetab_u_lrpneg_data),
                           rownames(musclemetab_u_lrppos_data))

adipose_combo_rownames <- c(rownames(adiposemetab_t_amine_data),
                           rownames(adiposemetab_t_nuc_data),
                           rownames(adiposemetab_t_oxylipneg_data),
                           rownames(adiposemetab_t_tca_data),
                           rownames(adiposemetab_u_hilicpos_data),
                           rownames(adiposemetab_u_ionpneg_data),
                           rownames(adiposemetab_u_rpneg_data),
                           rownames(adiposemetab_u_rppos_data),
                           rownames(adiposemetab_u_lrpneg_data),
                           rownames(adiposemetab_u_lrppos_data))

blood_combo_rownames <- c(rownames(bloodmetab_t_amine_data),
                           rownames(bloodmetab_t_oxylipneg_data),
                           rownames(bloodmetab_t_tca_data),
                           rownames(bloodmetab_t_conv_data),
                           rownames(bloodmetab_t_imm_crt_data),
                           rownames(bloodmetab_u_hilicpos_data),
                           rownames(bloodmetab_u_ionpneg_data),
                           rownames(bloodmetab_u_rpneg_data),
                           rownames(bloodmetab_u_rppos_data),
                           rownames(bloodmetab_u_lrpneg_data),
                           rownames(bloodmetab_u_lrppos_data))

musclemetab_combo_data <- rbind(musclemetab_t_amine_data[,musclesamplelist],
                                musclemetab_t_oxylipneg_data[,musclesamplelist],
                                musclemetab_t_tca_data[,musclesamplelist],
                                musclemetab_u_hilicpos_data[,musclesamplelist],
                                musclemetab_u_ionpneg_data[,musclesamplelist],
                                musclemetab_u_rpneg_data[,musclesamplelist],
                                musclemetab_u_rppos_data[,musclesamplelist],
                                musclemetab_u_lrpneg_data[,musclesamplelist],
                                musclemetab_u_lrppos_data[,musclesamplelist])

bloodmetab_combo_data <- rbind(bloodmetab_t_amine_data[,bloodsamplelist],
                               bloodmetab_t_oxylipneg_data[,bloodsamplelist],
                               bloodmetab_t_tca_data[,bloodsamplelist],
                               bloodmetab_u_hilicpos_data[,bloodsamplelist],
                               bloodmetab_u_ionpneg_data[,bloodsamplelist],
                               bloodmetab_u_rpneg_data[,bloodsamplelist],
                               bloodmetab_u_rppos_data[,bloodsamplelist],
                               bloodmetab_u_lrpneg_data[,bloodsamplelist],
                               bloodmetab_u_lrppos_data[,bloodsamplelist])

adiposemetab_combo_data <- rbind(adiposemetab_t_amine_data[,adiposesamplelist],
                                 adiposemetab_t_oxylipneg_data[,adiposesamplelist],
                                 adiposemetab_t_tca_data[,adiposesamplelist],
                                 adiposemetab_u_hilicpos_data[,adiposesamplelist],
                                 adiposemetab_u_ionpneg_data[,adiposesamplelist],
                                 adiposemetab_u_rpneg_data[,adiposesamplelist],
                                 adiposemetab_u_rppos_data[,adiposesamplelist],
                                 adiposemetab_u_lrpneg_data[,adiposesamplelist],
                                 adiposemetab_u_lrppos_data[,adiposesamplelist])

# We now identify metabolites that are measured in all three tissues so we can generate a cross-tissue comparison matrix
musclemetab_combo_data <- musclemetab_combo_data[intersect(rownames(musclemetab_combo_data),Reduce(intersect,list(muscle_combo_rownames,adipose_combo_rownames,blood_combo_rownames))),]
bloodmetab_combo_data <- bloodmetab_combo_data[intersect(rownames(bloodmetab_combo_data),Reduce(intersect,list(muscle_combo_rownames,adipose_combo_rownames,blood_combo_rownames))),]
adiposemetab_combo_data <- adiposemetab_combo_data[intersect(rownames(adiposemetab_combo_data),Reduce(intersect,list(muscle_combo_rownames,adipose_combo_rownames,blood_combo_rownames))),]

musclemetab_combo_data <- musclemetab_combo_data[Reduce(intersect,list(rownames(musclemetab_combo_data),
                                                                       rownames(bloodmetab_combo_data),
                                                                       rownames(adiposemetab_combo_data))),]
adiposemetab_combo_data <- adiposemetab_combo_data[Reduce(intersect,list(rownames(musclemetab_combo_data),
                                                                       rownames(bloodmetab_combo_data),
                                                                       rownames(adiposemetab_combo_data))),]
bloodmetab_combo_data <- bloodmetab_combo_data[Reduce(intersect,list(rownames(musclemetab_combo_data),
                                                                       rownames(bloodmetab_combo_data),
                                                                       rownames(adiposemetab_combo_data))),]

# We individually z-score each tissue's normalized data matrix prior to concatenation
musclemetab_combo_dataz <- t(scale(t(musclemetab_combo_data)))
bloodmetab_combo_dataz <- t(scale(t(bloodmetab_combo_data)))
adiposemetab_combo_dataz <- t(scale(t(adiposemetab_combo_data)))

colnames(musclemetab_combo_dataz) <- paste("muscle_",colnames(musclemetab_combo_dataz),sep = "")
colnames(bloodmetab_combo_dataz) <- paste("blood_",colnames(bloodmetab_combo_dataz),sep = "")
colnames(adiposemetab_combo_dataz) <- paste("adipose_",colnames(adiposemetab_combo_dataz),sep = "")
crosstissue_metab_dataz <- cbind(musclemetab_combo_dataz,bloodmetab_combo_dataz,adiposemetab_combo_dataz)

# Generating sample metadata matricesfor each tissue to combine
blood_metab_samplemetadata <- data.frame(row.names = bloodsamplelist,
                                         "Tissue" = rep("Blood",length(bloodsamplelist)),
                                         "Timepoint" = gsub(".*post","post",gsub(".*during","during",gsub(".*pre","pre",bloodsamplelist))),
                                         "pid" = gsub("_.*","",bloodsamplelist))
blood_metab_samplemetadata$Sex = "Female"
blood_metab_samplemetadata$randomGroupCode = "Control"
for(i in 1:dim(blood_metab_samplemetadata)[1]){
  ourpid <- blood_metab_samplemetadata$pid[i]
  blood_metab_samplemetadata[i,"Sex"] <- c("Female","Male")[MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$sample_metadata[MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$sample_metadata$pid %in% ourpid,"Sex"][1]]
  blood_metab_samplemetadata[i,"randomGroupCode"] <- MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$sample_metadata[MotrpacHumanPreSuspensionData::BLOOD_METAB_T_AMINES_QC$sample_metadata$pid %in% ourpid,"randomGroupCode"][1]
}

adipose_metab_samplemetadata <- data.frame(row.names = adiposesamplelist,
                                           "Tissue" = rep("Adipose",length(adiposesamplelist)),
                                           "Timepoint" = gsub(".*post","post",gsub(".*during","during",gsub(".*pre","pre",adiposesamplelist))),
                                           "pid" = gsub("_.*","",adiposesamplelist))
adipose_metab_samplemetadata$Sex = "Female"
adipose_metab_samplemetadata$randomGroupCode = "Control"
for(i in 1:dim(adipose_metab_samplemetadata)[1]){
  ourpid <- adipose_metab_samplemetadata$pid[i]
  adipose_metab_samplemetadata[i,"Sex"] <- c("Female","Male")[MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$sample_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$sample_metadata$pid %in% ourpid,"Sex"][1]]
  adipose_metab_samplemetadata[i,"randomGroupCode"] <- MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$sample_metadata[MotrpacHumanPreSuspensionData::ADIPOSE_METAB_T_AMINES_QC$sample_metadata$pid %in% ourpid,"randomGroupCode"][1]
}

muscle_metab_samplemetadata <- data.frame(row.names = musclesamplelist,
                                          "Tissue" = rep("Muscle",length(musclesamplelist)),
                                          "Timepoint" = gsub(".*post","post",gsub(".*during","during",gsub(".*pre","pre",musclesamplelist))),
                                          "pid" = gsub("_.*","",musclesamplelist))
muscle_metab_samplemetadata$Sex = "Female"
muscle_metab_samplemetadata$randomGroupCode = "Control"
for(i in 1:dim(muscle_metab_samplemetadata)[1]){
  ourpid <- muscle_metab_samplemetadata$pid[i]
  muscle_metab_samplemetadata[i,"Sex"] <- c("Female","Male")[MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$sample_metadata[MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$sample_metadata$pid %in% ourpid,"Sex"][1]]
  muscle_metab_samplemetadata[i,"randomGroupCode"] <- MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$sample_metadata[MotrpacHumanPreSuspensionData::MUSCLE_METAB_T_AMINES_QC$sample_metadata$pid %in% ourpid,"randomGroupCode"][1]
}

rownames(blood_metab_samplemetadata) <- paste("blood_",rownames(blood_metab_samplemetadata),sep = "")
rownames(adipose_metab_samplemetadata) <- paste("adipose_",rownames(adipose_metab_samplemetadata),sep = "")
rownames(muscle_metab_samplemetadata) <- paste("muscle_",rownames(muscle_metab_samplemetadata),sep = "")

crosstissue_metab_samplemetadata <- rbind(muscle_metab_samplemetadata[intersect(rownames(muscle_metab_samplemetadata),colnames(crosstissue_metab_dataz)),],
                                          blood_metab_samplemetadata[intersect(rownames(blood_metab_samplemetadata),colnames(crosstissue_metab_dataz)),],
                                          adipose_metab_samplemetadata[intersect(rownames(adipose_metab_samplemetadata),colnames(crosstissue_metab_dataz)),])

# Editing metadata names to be in keeping with desired figure-ready labels
crosstissue_metab_samplemetadata$Timepoint <- as.character(crosstissue_metab_samplemetadata$Timepoint)

crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "pre_exercise","Timepoint"] <- "Pre"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "during_20_min","Timepoint"] <- "D20M"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "during_40_min","Timepoint"] <- "D40M"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "post_10_min","Timepoint"] <- "P10M"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "post_15_30_45_min","Timepoint"] <- "P15-45M"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "post_3.5_4_hr","Timepoint"] <- "P3.5/4H"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$Timepoint %in% "post_24_hr","Timepoint"] <- "P24H"

crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$randomGroupCode %in% "ADUResist","randomGroupCode"] <- "RE"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$randomGroupCode %in% "ADUEndur","randomGroupCode"] <- "EE"
crosstissue_metab_samplemetadata[crosstissue_metab_samplemetadata$randomGroupCode %in% "ADUControl","randomGroupCode"] <- "CON"

names(crosstissue_metab_samplemetadata)[5] <- "Modality"


#metab_featuremetadata <- data.frame(row.names = rownames(bloodmetab_combo_data),
#                                    "Assay" = c(rep("t_amine",length(t_amine_crosstissue)),
#                                                rep("t_oxylipneg",length(t_oxylipneg_crosstissue)),
#                                                rep("t_tca",length(t_tca_crosstissue)),
#                                                rep("u_hilicpos",length(u_hilicpos_crosstissue)),
#                                                rep("u_ionpneg",length(u_ionpneg_crosstissue)),
#                                                rep("u_rpneg",length(u_rpneg_crosstissue)),
#                                                rep("u_rppos",length(u_rppos_crosstissue)),
#                                                rep("u_lrpneg",length(u_lrpneg_crosstissue)),
#                                                rep("u_lrppos",length(u_lrppos_crosstissue))))
#metab_featurematrix <- matrix(0L,nrow = length(rownames(bloodmetab_combo_data)),ncol = 9)
#rownames(metab_featurematrix) <- rownames(bloodmetab_combo_data)
#colnames(metab_featurematrix) <- unique(metab_featuremetadata$Assay)
#metab_featurematrix[t_amine_crosstissue,"t_amine"] <- 1
#metab_featurematrix[t_oxylipneg_crosstissue,"t_oxylipneg"] <- 1
#metab_featurematrix[t_tca_crosstissue,"t_tca"] <- 1
#metab_featurematrix[u_hilicpos_crosstissue,"u_hilicpos"] <- 1
#metab_featurematrix[u_ionpneg_crosstissue,"u_ionpneg"] <- 1
#metab_featurematrix[u_rpneg_crosstissue,"u_rpneg"] <- 1
#metab_featurematrix[u_rppos_crosstissue,"u_rppos"] <- 1
#metab_featurematrix[u_lrpneg_crosstissue,"u_lrpneg"] <- 1
#metab_featurematrix[u_lrppos_crosstissue,"u_lrppos"] <- 1

# Converting the refmetmatrix into a prior knowledge format for running PLIER
refmetmatrix <- matrix(0L,nrow = length(rownames(bloodmetab_combo_data)),ncol = dim(refmetmat)[2])
rownames(refmetmatrix) <- rownames(bloodmetab_combo_data)
colnames(refmetmatrix) <- colnames(refmetmat)
for(i in 1:dim(refmetmatrix)[1]){
  ourmetab <- rownames(refmetmatrix)[i]
  if(ourmetab %in% rownames(refmetmat)){
    refmetmatrix[i,] <- refmetmat[ourmetab,]
  }
}

metab_finalpaths <- refmetmatrix

# Specifying the number of LVs we desire from PLIER. We typically use as input double the calculated number of PCs from the num.pc function
crosstissuemetabcombovalue <- num.pc(as.data.frame(crosstissue_metab_dataz), seed = 1)

# Running PLIER with default parameters
crosstissuemetabcombofin.plierResult.all.7=PLIER(as.matrix(crosstissue_metab_dataz),metab_finalpaths,k=2*crosstissuemetabcombovalue ,trace=F, frac=0.7, scale=T, seed=1)
# Saving the PLIER output
saveRDS(crosstissuemetabcombofin.plierResult.all.7,file = "crosstissuemetabplier_102925.RDS")

# Collecting color annotations for figure generation
precawg_ann_cols <- list("Tissue" = MotrpacHumanPreSuspension::HUMAN_TISSUE_COLORS,
                         "Sex" = MotrpacHumanPreSuspension::HUMAN_SEX_COLORS[c("Male","Female")],
                         "Timepoint" = MotrpacHumanPreSuspension::HUMAN_ACUTE_TIMEPOINT_COLORS[c("Pre_exercise","during_20_min","during_40_min","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr")],
                         "Modality" = MotrpacHumanPreSuspension::HUMAN_EXERCISE_GROUP_COLORS)
names(precawg_ann_cols$Tissue) <- c("Blood","Muscle","Adipose")
names(precawg_ann_cols$Timepoint) <- c("Pre","D20M","D40M","P10M","P15-45M","P3.5/4H","P24H")
names(precawg_ann_cols$Modality) <- c("RE","EE","CON")

crosstissue_metab_samplemetadata$Timepoint <- factor(crosstissue_metab_samplemetadata$Timepoint,levels = c("Pre",
                                                                                                           "D20M",
                                                                                                           "D40M",
                                                                                                           "P10M",
                                                                                                           "P15-45M",
                                                                                                           "P3.5/4H",
                                                                                                           "P24H"))
crosstissue_metab_samplemetadata$Tissue <- factor(crosstissue_metab_samplemetadata$Tissue,levels = c("Adipose","Blood","Muscle"))
crosstissue_metab_samplemetadata <- crosstissue_metab_samplemetadata[order(crosstissue_metab_samplemetadata$Tissue,crosstissue_metab_samplemetadata$Modality,crosstissue_metab_samplemetadata$Timepoint,crosstissue_metab_samplemetadata$Sex),]
crosstissue_metab_dataz <- crosstissue_metab_dataz[,rownames(crosstissue_metab_samplemetadata)]


# Generating heatmaps of z-scored expression across all samples for top metabolites associated with each LV
# Heatmaps of select LVs of interest are part of Supplemental Figure S5
for(i in 1:dim(crosstissuemetabcombofin.plierResult.all.7$Z)[2]){

  top10rows <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,i]),][1:10,])
  datatoplot <- crosstissue_metab_dataz[top10rows,]

  png(filename=paste("crosstissuemetab_PLIER_LV",toString(i),"Plot102925.png",sep=""),width = 8,height = 3,units = "in",res = 600)

  pheatmap(datatoplot,breaks = seq(-5,5, length.out = 100),color = colorpanel(101,"blue", "white", "red"),show_colnames=F,cluster_cols = FALSE,cluster_rows = FALSE,annotation_col = crosstissue_metab_samplemetadata[,c("Sex","Timepoint","Modality","Tissue")],fontsize=14,annotation_colors = precawg_ann_cols,annotation_legend = FALSE)

  dev.off()
}


crosstissue_metab_samplemetadata$pid <- factor(crosstissue_metab_samplemetadata$pid,levels = unique(crosstissue_metab_samplemetadata$pid))

trimcrosstissuemetabmetadata <- crosstissue_metab_samplemetadata
#trimcrosstissuemetabmetadata <- crosstissuemetabmetadata[!(crosstissuemetabmetadata$pid %in% c("11989394","11326794","10480749","13098039","15783162","16196910","10484493","11258187","12279390","")),]

# Identifying the pids for subjects with data in the different tissues
adiposemetabsubjects <- rownames(table(trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Adipose","pid"],trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Adipose","Timepoint"]))[table(trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Adipose","pid"],trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Adipose","Timepoint"])[,"Pre"] == 1]
musclemetabsubjects <- rownames(table(trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Muscle","pid"],trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Muscle","Timepoint"]))[table(trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Muscle","pid"],trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Muscle","Timepoint"])[,"Pre"] == 1]
bloodmetabsubjects <- rownames(table(trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Blood","pid"],trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Blood","Timepoint"]))[table(trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Blood","pid"],trimcrosstissuemetabmetadata[trimcrosstissuemetabmetadata$Tissue %in% "Blood","Timepoint"])[,"Pre"] == 1]

#crosstissuemetabBmat <- crosstissuemetabcombofin.plierResult.all.7$B
#rownames(crosstissuemetabBmat) <- paste("LV",c(1:dim(crosstissuemetabBmat)[1]),sep = "")

#####
# We want to make similar figures to what was made for the cross tissue RNA
####

# First we need to update how marker lists are generated.

crosstissuemetabcombofin_Zmatscore <- t(scale(t(crosstissuemetabcombofin.plierResult.all.7$Z)))
crosstissuemetablv11markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,11]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[11],])
crosstissuemetablv12markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,12]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[12],])
crosstissuemetablv17markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,17]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[17],])
crosstissuemetablv27markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,27]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[27],])
crosstissuemetablv38markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,38]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[38],])
crosstissuemetablv40markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,40]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[40],])
crosstissuemetablv87markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,87]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[87],])

crosstissuemetablvmarkerlist <- list("LV11" = crosstissuemetablv11markers,
                                     "LV12" = crosstissuemetablv12markers,
                                     "LV17" = crosstissuemetablv17markers,
                                     "LV27" = crosstissuemetablv27markers,
                                     "LV38" = crosstissuemetablv38markers,
                                     "LV40" = crosstissuemetablv40markers,
                                     "LV87" = crosstissuemetablv87markers)

#####
# Now, we want to plot how similar/different the responses are in the three tissues
# We will get a score for relationship/correlation between each tissue
# First we need to divide the B matrix into three matrices for the tissues
# Then we need to rename the matrix columns to reflect pid, timepoint
####

crosstissuemetabbmat <- crosstissuemetabcombofin.plierResult.all.7$B
rownames(crosstissuemetabbmat) <- paste("LV",c(1:88),sep = "")
crosstissuemetabbmat_muscle <- crosstissuemetabbmat[,grep("muscle",colnames(crosstissuemetabbmat))]
crosstissuemetabbmat_adipose <- crosstissuemetabbmat[,grep("adipose",colnames(crosstissuemetabbmat))]
crosstissuemetabbmat_blood <- crosstissuemetabbmat[,grep("blood",colnames(crosstissuemetabbmat))]

muscle_metab_samplemetadata$pid_timepoint <- paste(muscle_metab_samplemetadata$pid,muscle_metab_samplemetadata$Timepoint,sep = "_")
colnames(crosstissuemetabbmat_muscle) <- muscle_metab_samplemetadata[colnames(crosstissuemetabbmat_muscle),"pid_timepoint"]

adipose_metab_samplemetadata$pid_timepoint <- paste(adipose_metab_samplemetadata$pid,adipose_metab_samplemetadata$Timepoint,sep = "_")
colnames(crosstissuemetabbmat_adipose) <- adipose_metab_samplemetadata[colnames(crosstissuemetabbmat_adipose),"pid_timepoint"]

blood_metab_samplemetadata$pid_timepoint <- paste(blood_metab_samplemetadata$pid,blood_metab_samplemetadata$Timepoint,sep = "_")
colnames(crosstissuemetabbmat_blood) <- blood_metab_samplemetadata[colnames(crosstissuemetabbmat_blood),"pid_timepoint"]

crosstissuemetabbmat_tissuecordf <- data.frame(row.names = rownames(crosstissuemetabbmat),
                                             "Adipose_v_Blood" = rep(0,length(rownames(crosstissuemetabbmat))),
                                             "Adipose_v_Muscle" = rep(0,length(rownames(crosstissuemetabbmat))),
                                             "Blood_v_Muscle" = rep(0,length(rownames(crosstissuemetabbmat))))
for(i in 1:dim(crosstissuemetabbmat_tissuecordf)[1]){
  crosstissuemetabbmat_tissuecordf[i,"Adipose_v_Blood"] <- cor(crosstissuemetabbmat_adipose[i,intersect(colnames(crosstissuemetabbmat_adipose),colnames(crosstissuemetabbmat_blood))],
                                                             crosstissuemetabbmat_blood[i,intersect(colnames(crosstissuemetabbmat_adipose),colnames(crosstissuemetabbmat_blood))])
  crosstissuemetabbmat_tissuecordf[i,"Adipose_v_Muscle"] <- cor(crosstissuemetabbmat_adipose[i,intersect(colnames(crosstissuemetabbmat_adipose),colnames(crosstissuemetabbmat_muscle))],
                                                              crosstissuemetabbmat_muscle[i,intersect(colnames(crosstissuemetabbmat_adipose),colnames(crosstissuemetabbmat_muscle))])
  crosstissuemetabbmat_tissuecordf[i,"Blood_v_Muscle"] <- cor(crosstissuemetabbmat_blood[i,intersect(colnames(crosstissuemetabbmat_blood),colnames(crosstissuemetabbmat_muscle))],
                                                            crosstissuemetabbmat_muscle[i,intersect(colnames(crosstissuemetabbmat_blood),colnames(crosstissuemetabbmat_muscle))])
}
crosstissuemetabbmat_tissuecordf$LV <- rownames(crosstissuemetabbmat_tissuecordf)
crosstissuemetabbmat_tissuecordf$SigLV <- 1 + 1*(crosstissuemetabbmat_tissuecordf$LV %in% names(crosstissuemetablvmarkerlist))

#with(crosstissuemetabbmat_tissuecordf,plot3d(Adipose_v_Blood,Adipose_v_Muscle,Blood_v_Muscle))
#with(crosstissuemetabbmat_tissuecordf,text3d(x = Adipose_v_Blood, y = Adipose_v_Muscle,z = Blood_v_Muscle,texts = LV,col = SigLV))

crosstissuemetabbmat_tissuecordf_forggplot <- data.frame("Correlation" = c(crosstissuemetabbmat_tissuecordf$Adipose_v_Blood,
                                                                         crosstissuemetabbmat_tissuecordf$Adipose_v_Muscle,
                                                                         crosstissuemetabbmat_tissuecordf$Blood_v_Muscle),
                                                       "Comparison" = c(rep("Adipose vs Blood",length(crosstissuemetabbmat_tissuecordf$Adipose_v_Blood)),
                                                                        rep("Adipose vs Muscle",length(crosstissuemetabbmat_tissuecordf$Adipose_v_Muscle)),
                                                                        rep("Blood vs Muscle",length(crosstissuemetabbmat_tissuecordf$Blood_v_Muscle))),
                                                       "LV" = rep(rownames(crosstissuemetabbmat_tissuecordf),3))
crosstissuemetabbmat_tissuecordf_forggplot$LV <- factor(crosstissuemetabbmat_tissuecordf_forggplot$LV,levels = crosstissuemetabbmat_tissuecordf_forggplot$LV[order(crosstissuemetabbmat_tissuecordf$Blood_v_Muscle)])


crosstissuemetablvmembershipbarplotdf <- data.frame("LV" = c(rep("LV11",length(crosstissuemetablvmarkerlist$LV11)),
                                                           rep("LV12",length(crosstissuemetablvmarkerlist$LV12)),
                                                           rep("LV17",length(crosstissuemetablvmarkerlist$LV17)),
                                                           rep("LV27",length(crosstissuemetablvmarkerlist$LV27)),
                                                           rep("LV38",length(crosstissuemetablvmarkerlist$LV38)),
                                                           rep("LV40",length(crosstissuemetablvmarkerlist$LV40)),
                                                           rep("LV87",length(crosstissuemetablvmarkerlist$LV87))))

crosstissuemetablvmembershipbarplotdf$LV <- factor(crosstissuemetablvmembershipbarplotdf$LV,levels = crosstissuemetabbmat_tissuecordf_forggplot$LV[order(crosstissuemetabbmat_tissuecordf$Blood_v_Muscle)])

# Figure 5G
png(file = "Figure5G_crosstissuemetab_lvtissuecorrelation_102925.png",width = 4,height = 4,units = "in",res = 600)
ggplot(crosstissuemetabbmat_tissuecordf_forggplot[crosstissuemetabbmat_tissuecordf_forggplot$LV %in% names(crosstissuemetablvmarkerlist),]) + geom_hline(yintercept = 0,size = 1,linetype = "dashed") + geom_point(aes(x=LV,y=Correlation,color=Comparison),size = 3) + theme_classic() + scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7")) + theme(axis.title.x=element_blank(),legend.position = "top",text = element_text(size=10),legend.key.spacing = unit(0,"cm")) + guides(col = guide_legend(ncol = 1))
dev.off()

# Figure 5H
png(file = "Figure5H_crosstissuemetab_lvmembershipbarplot_102925.png",width = 4,height = 3,units = "in",res = 600)
ggplot(crosstissuemetablvmembershipbarplotdf,aes(LV)) + geom_bar() + theme_classic() + ylab("Membership") + theme(axis.title.x=element_blank(),text = element_text(size=10))
dev.off()

crosstissuemetabcombofinUaucplottrim <- crosstissuemetabcombofin.plierResult.all.7$Uauc[,c(11,12,17,27,38,40,87)]
crosstissuemetabcombofinUaucplottrim <- crosstissuemetabcombofinUaucplottrim[apply(crosstissuemetabcombofinUaucplottrim,1,max) > 0,]

# Figure 5I
png("Figure5I_CrossTissueMetab_UaucPlot_102925.png",width = 6,height = 4, units = "in",res = 600)
pheatmap(crosstissuemetabcombofinUaucplottrim[,levels(crosstissuemetabbmat_tissuecordf_forggplot$LV)[levels(crosstissuemetabbmat_tissuecordf_forggplot$LV) %in% names(crosstissuemetablvmarkerlist)]],breaks = seq(0,1,length.out = 101),color = colorpanel(101,"white","firebrick"),cluster_cols = F,angle_col = 0,labels_row = c("Acyl Carnitines","Amino Acids","Saturated FA","Unsaturated FA"))
dev.off()




#####
# Now, we want to make a plot of the significant (up or down) responses by each LV
# So we will need to use compare_means like we do in the boxplots to generate a p-val for each LV
# for each comparison - then create a simpler heatmap with a column for each Tissue/Time-Point/Endurance+Resistance Groups
####

crosstissuemetab_lvsigtestdf <- matrix(0L,nrow = 88,ncol = 22)
rownames(crosstissuemetab_lvsigtestdf) <- paste("LV",c(1:88),sep = "")
colnames(crosstissuemetab_lvsigtestdf) <- c("Adipose_45min_Endurance","Adipose_4hr_Endurance","Adipose_24hr_Endurance",
                                            "Adipose_45min_Resistance","Adipose_4hr_Resistance","Adipose_24hr_Resistance",
                                            "Blood_20minDuring_Endurance","Blood_40minDuring_Endurance","Blood_10min_Endurance",
                                            "Blood_30min_Endurance","Blood_4hr_Endurance","Blood_24hr_Endurance",
                                            "Blood_10min_Resistance","Blood_30min_Resistance","Blood_4hr_Resistance","Blood_24hr_Resistance",
                                            "Muscle_15min_Endurance","Muscle_3_5hr_Endurance","Muscle_24hr_Endurance",
                                            "Muscle_15min_Resistance","Muscle_3_5hr_Resistance","Muscle_24hr_Resistance")

for(i in 1:dim(crosstissuemetab_lvsigtestdf)[1]){

  metamerge <- trimcrosstissuemetabmetadata
  metamerge$B <- crosstissuemetabcombofin.plierResult.all.7$B[i,rownames(trimcrosstissuemetabmetadata)]
  metamerge$pid <- factor(metamerge$pid,levels = unique(metamerge$pid))
  metamerge$Diff <- 0

  metamergeadipose <- metamerge[metamerge$Tissue %in% "Adipose",]
  metamergemuscle <- metamerge[metamerge$Tissue %in% "Muscle",]
  metamergeblood <- metamerge[metamerge$Tissue %in% "Blood",]

  metamergeadipose <- metamergeadipose[metamergeadipose$pid %in% adiposemetabsubjects,]
  metamergemuscle <- metamergemuscle[metamergemuscle$pid %in% musclemetabsubjects,]
  metamergeblood <- metamergeblood[metamergeblood$pid %in% bloodmetabsubjects,]

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

  crosstissuemetab_lvsigtestdf[i,"Adipose_45min_Endurance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "EE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[1,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Adipose_45min_Resistance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "RE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[2,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Adipose_4hr_Endurance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "EE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[3,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Adipose_4hr_Resistance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "RE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[4,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Adipose_24hr_Endurance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "EE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[5,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Adipose_24hr_Resistance"] <- sign(median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "RE","Diff"])-median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"]))*(-log10(adiposetestout[6,"p.adj"][[1]]))

  crosstissuemetab_lvsigtestdf[i,"Muscle_15min_Endurance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "EE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[1,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Muscle_15min_Resistance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "RE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[2,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Muscle_3_5hr_Endurance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "EE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[3,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Muscle_3_5hr_Resistance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "RE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[4,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Muscle_24hr_Endurance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "EE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[5,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Muscle_24hr_Resistance"] <- sign(median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "RE","Diff"])-median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"]))*(-log10(muscletestout[6,"p.adj"][[1]]))

  crosstissuemetab_lvsigtestdf[i,"Blood_20minDuring_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[1,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_40minDuring_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[2,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_10min_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[3,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_10min_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[4,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_30min_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[5,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_30min_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[6,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_4hr_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[7,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_4hr_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[8,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_24hr_Endurance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "EE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[9,"p.adj"][[1]]))
  crosstissuemetab_lvsigtestdf[i,"Blood_24hr_Resistance"] <- sign(median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "RE","Diff"])-median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"]))*(-log10(bloodtestout[10,"p.adj"][[1]]))

}

crosstissuemetab_lvsigtestdfmeta <- data.frame(row.names = colnames(crosstissuemetab_lvsigtestdf),
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

crosstissuemetab_lvsigtestdfmeta$Tissue <- factor(crosstissuemetab_lvsigtestdfmeta$Tissue,levels = c("Adipose","Blood","Muscle"))
crosstissuemetab_lvsigtestdfmeta$Timepoint <- factor(crosstissuemetab_lvsigtestdfmeta$Timepoint,levels = c("D20M","D40M","P10M","P15-45M","P3.5/4H","P24H"))
crosstissuemetab_lvsigtestdfmeta$Modality <- factor(crosstissuemetab_lvsigtestdfmeta$Modality,levels = c("EE","RE"))

crosstissuemetab_lvsigtestdfmeta <- crosstissuemetab_lvsigtestdfmeta[order(crosstissuemetab_lvsigtestdfmeta$Tissue,crosstissuemetab_lvsigtestdfmeta$Modality,crosstissuemetab_lvsigtestdfmeta$Timepoint),]
crosstissuemetab_lvsigtestdf <- crosstissuemetab_lvsigtestdf[,rownames(crosstissuemetab_lvsigtestdfmeta)]


# Generating a heatmap detailing the signficant responses to exercise by tissue, mode and time pointfor each targeted LV
# We currently use this for the annotations to add to the bubble heatmap
png(file = "crosstissuemetab_toplvs_significanceheatmap_102925.png",width = 8,height = 6,units = "in",res = 600)
pheatmap(crosstissuemetab_lvsigtestdf[names(crosstissuemetablvmarkerlist),],cluster_cols = F,breaks = seq(-5,5,length.out = 101),color = colorpanel(101,"blue","white","red"),annotation_col = crosstissuemetab_lvsigtestdfmeta[,c("Timepoint","Modality","Tissue")],annotation_colors = precawg_ann_cols,show_colnames = F)
dev.off()

#####
# Similar to above, but now we generate a finalized form with bubbleheatmap
####

crosstissuemetab_lvsigtest_bhsize_df <- matrix(0L,nrow = 88,ncol = 22)
rownames(crosstissuemetab_lvsigtest_bhsize_df) <- paste("LV",c(1:88),sep = "")
colnames(crosstissuemetab_lvsigtest_bhsize_df) <- c("Adipose_45min_Endurance","Adipose_4hr_Endurance","Adipose_24hr_Endurance",
                                                    "Adipose_45min_Resistance","Adipose_4hr_Resistance","Adipose_24hr_Resistance",
                                                    "Blood_20minDuring_Endurance","Blood_40minDuring_Endurance","Blood_10min_Endurance",
                                                    "Blood_30min_Endurance","Blood_4hr_Endurance","Blood_24hr_Endurance",
                                                    "Blood_10min_Resistance","Blood_30min_Resistance","Blood_4hr_Resistance","Blood_24hr_Resistance",
                                                    "Muscle_15min_Endurance","Muscle_3_5hr_Endurance","Muscle_24hr_Endurance",
                                                    "Muscle_15min_Resistance","Muscle_3_5hr_Resistance","Muscle_24hr_Resistance")

crosstissuemetab_lvsigtest_bhcolor_df <- matrix(0L,nrow = 88,ncol = 22)
rownames(crosstissuemetab_lvsigtest_bhcolor_df) <- paste("LV",c(1:88),sep = "")
colnames(crosstissuemetab_lvsigtest_bhcolor_df) <- c("Adipose_45min_Endurance","Adipose_4hr_Endurance","Adipose_24hr_Endurance",
                                                     "Adipose_45min_Resistance","Adipose_4hr_Resistance","Adipose_24hr_Resistance",
                                                     "Blood_20minDuring_Endurance","Blood_40minDuring_Endurance","Blood_10min_Endurance",
                                                     "Blood_30min_Endurance","Blood_4hr_Endurance","Blood_24hr_Endurance",
                                                     "Blood_10min_Resistance","Blood_30min_Resistance","Blood_4hr_Resistance","Blood_24hr_Resistance",
                                                     "Muscle_15min_Endurance","Muscle_3_5hr_Endurance","Muscle_24hr_Endurance",
                                                     "Muscle_15min_Resistance","Muscle_3_5hr_Resistance","Muscle_24hr_Resistance")

for(i in 1:dim(crosstissuemetab_lvsigtestdf)[1]){

  metamerge <- trimcrosstissuemetabmetadata
  metamerge$B <- crosstissuemetabcombofin.plierResult.all.7$B[i,rownames(trimcrosstissuemetabmetadata)]
  metamerge$pid <- factor(metamerge$pid,levels = unique(metamerge$pid))
  metamerge$Diff <- 0

  metamergeadipose <- metamerge[metamerge$Tissue %in% "Adipose",]
  metamergemuscle <- metamerge[metamerge$Tissue %in% "Muscle",]
  metamergeblood <- metamerge[metamerge$Tissue %in% "Blood",]

  metamergeadipose <- metamergeadipose[metamergeadipose$pid %in% adiposemetabsubjects,]
  metamergemuscle <- metamergemuscle[metamergemuscle$pid %in% musclemetabsubjects,]
  metamergeblood <- metamergeblood[metamergeblood$pid %in% bloodmetabsubjects,]

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

  crosstissuemetab_lvsigtest_bhsize_df[i,"Adipose_45min_Endurance"] <- -log10(adiposetestout[1,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Adipose_45min_Resistance"] <- -log10(adiposetestout[2,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Adipose_4hr_Endurance"] <- -log10(adiposetestout[3,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Adipose_4hr_Resistance"] <- -log10(adiposetestout[4,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Adipose_24hr_Endurance"] <- -log10(adiposetestout[5,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Adipose_24hr_Resistance"] <- -log10(adiposetestout[6,"p.adj"][[1]])

  crosstissuemetab_lvsigtest_bhsize_df[i,"Muscle_15min_Endurance"] <- -log10(muscletestout[1,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Muscle_15min_Resistance"] <- -log10(muscletestout[2,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Muscle_3_5hr_Endurance"] <- -log10(muscletestout[3,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Muscle_3_5hr_Resistance"] <- -log10(muscletestout[4,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Muscle_24hr_Endurance"] <- -log10(muscletestout[5,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Muscle_24hr_Resistance"] <- -log10(muscletestout[6,"p.adj"][[1]])

  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_20minDuring_Endurance"] <- -log10(bloodtestout[1,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_40minDuring_Endurance"] <- -log10(bloodtestout[2,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_10min_Endurance"] <- -log10(bloodtestout[3,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_10min_Resistance"] <- -log10(bloodtestout[4,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_30min_Endurance"] <- -log10(bloodtestout[5,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_30min_Resistance"] <- -log10(bloodtestout[6,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_4hr_Endurance"] <- -log10(bloodtestout[7,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_4hr_Resistance"] <- -log10(bloodtestout[8,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_24hr_Endurance"] <- -log10(bloodtestout[9,"p.adj"][[1]])
  crosstissuemetab_lvsigtest_bhsize_df[i,"Blood_24hr_Resistance"] <- -log10(bloodtestout[10,"p.adj"][[1]])

  crosstissuemetab_lvsigtest_bhcolor_df[i,"Adipose_45min_Endurance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "EE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Adipose_45min_Resistance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "RE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P15-45M" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Adipose_4hr_Endurance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "EE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Adipose_4hr_Resistance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "RE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P3.5/4H" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Adipose_24hr_Endurance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "EE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Adipose_24hr_Resistance"] <- median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "RE","Diff"]) - median(metamergeadipose[metamergeadipose$Timepoint %in% "P24H" & metamergeadipose$Modality %in% "CON","Diff"])

  crosstissuemetab_lvsigtest_bhcolor_df[i,"Muscle_15min_Endurance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "EE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Muscle_15min_Resistance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "RE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P15-45M" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Muscle_3_5hr_Endurance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "EE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Muscle_3_5hr_Resistance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "RE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P3.5/4H" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Muscle_24hr_Endurance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "EE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Muscle_24hr_Resistance"] <- median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "RE","Diff"]) - median(metamergemuscle[metamergemuscle$Timepoint %in% "P24H" & metamergemuscle$Modality %in% "CON","Diff"])

  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_20minDuring_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "D20M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_40minDuring_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "D40M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_10min_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_10min_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P10M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_30min_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_30min_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P15-45M" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_4hr_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_4hr_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P3.5/4H" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_24hr_Endurance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "EE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"])
  crosstissuemetab_lvsigtest_bhcolor_df[i,"Blood_24hr_Resistance"] <- median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "RE","Diff"]) - median(metamergeblood[metamergeblood$Timepoint %in% "P24H" & metamergeblood$Modality %in% "CON","Diff"])
}

testout <- bubbleHeatmap(colorMat = crosstissuemetab_lvsigtest_bhcolor_df[levels(crosstissuemetabbmat_tissuecordf_forggplot$LV)[levels(crosstissuemetabbmat_tissuecordf_forggplot$LV) %in% names(crosstissuemetablvmarkerlist)],],sizeMat = crosstissuemetab_lvsigtest_bhsize_df[levels(crosstissuemetabbmat_tissuecordf_forggplot$LV)[levels(crosstissuemetabbmat_tissuecordf_forggplot$LV) %in% names(crosstissuemetablvmarkerlist)],],legendTitles = c("Significance","log2FC"),colorSeq = c("blue","white","red"),colorLim = c(-0.5,0.5),sizeLim = c(0,5))

# Figure 5F
png(file = "Figure5F_crosstissuemetab_toplvs_bubbleheatmap_102925.png",width = 10,height = 10,units = "in",res = 600)
grid.draw(testout)
dev.off()

#####
# Making the Metabolomics portion of the Supplemental Table S5
####

ourindexCol <- as.numeric(gsub("LV","",names(crosstissuemetablvmarkerlist)))

toplvZmat <- crosstissuemetabcombofin.plierResult.all.7$Z[Reduce(union,crosstissuemetablvmarkerlist),ourindexCol]
colnames(toplvZmat) <- paste("LV",ourindexCol,sep = "")
toplvUmat <- crosstissuemetabcombofin.plierResult.all.7$Uauc[,ourindexCol]
toplvUmat <- toplvUmat[apply(toplvUmat,1,max) > 0,]

MetabPLIER_supplementaltable <- rbind(toplvUmat,toplvZmat)
write.csv(MetabPLIER_supplementaltable,file = "Supplemental Table S5_Metab_102925.csv",row.names = T)

# Figure 5J
MotrpacHumanPreSuspension::plot_single_feature(feature = "CAR(8:0)",repo_local_dir = "~/GitHub/precovid-analyses/",epigen = F,selected_tissues = "muscle",selected_omes = "metab",output_file = "Figure5J_muscle_CAR8_0_metabpeak_plot_093025.png")
# Figure 5K
MotrpacHumanPreSuspension::plot_single_feature(feature = "CAR(8:0)",repo_local_dir = "~/GitHub/precovid-analyses/",epigen = F,selected_tissues = "blood",selected_omes = "metab",output_file = "Figure5K_blood_CAR8_0_metabpeak_plot_093025.png")

# In case you want to save progress
#save.image("precovid_metabplieranalysis_102925.RData")
