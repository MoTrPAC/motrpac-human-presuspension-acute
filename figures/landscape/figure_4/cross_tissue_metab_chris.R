library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacBicQC)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(gprofiler2)
library(ggpubr)
library(PLIER)
library(gridExtra)
library(biomaRt)
library(circlize)
library(bubbleHeatmap)
library(here)

# library(devtools)
# install_github("wgmao/PLIER")

setwd(file.path(here(), "figures", "landscape", "figure_5"))

# Loading the refmet metabolite classifications for use in PLIER's prior knowledge
refmetmetabolites <- Reduce(union,MOLECULAR_SIGNATURES$REFMET)
refmetmat <- matrix(0L,nrow = length(refmetmetabolites),ncol = length(names(MOLECULAR_SIGNATURES$REFMET)))
rownames(refmetmat) <- refmetmetabolites
colnames(refmetmat) <- names(MOLECULAR_SIGNATURES$REFMET)
for(i in 1:length(names(MOLECULAR_SIGNATURES$REFMET))){
  refmetmat[MOLECULAR_SIGNATURES$REFMET[i][[1]],i] <- 1
}

#in order to match what greg did for the rnaseq we don't want to directly
#combine all tissues then pivot longer, because we dont want to subset to
#all shared timepoints, its ok if some tissues only measure some tissues

combined_metab_list = list()
for(tissue in tissue_available_list()){
  metab_tissue = load_qc(selected_tissues = tissue,
                         selected_omes = "metab")
  combined_metab = combine_qc_matrixes(metab_tissue)
  rownames(combined_metab) = sub("\\..*", "", rownames(combined_metab))
  colnames(combined_metab) = paste(colnames(combined_metab), tissue, sep = "-")

  combined_metab_list[[tissue]] = combined_metab %>%
    t() %>%
    scale() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature_id") #and also scale per matrix here
}

combined_metab_joined = Reduce(function(x, y) inner_join(x, y, by = "feature_id"),
                               combined_metab_list) %>%
  tibble::column_to_rownames("feature_id")
#just subset the text before the "_" which is the pid + tp
pid_tp_df = data.frame(
  colname = colnames(combined_metab_joined),
  pid_tp = sub("-.*", "", colnames(combined_metab_joined))
)


crosstissue_metab_samplemetadata = load_pheno()$data %>%
  mutate(pid_tp = paste(pid, Timepoint, sep = "..")) %>%
  left_join(., pid_tp_df, by = "pid_tp") %>%
  filter(!is.na(colname)) %>%
  dplyr::select(pid, Timepoint, Sex, randomGroupCode, colname) %>%
  distinct(colname, .keep_all = TRUE) %>%
  dplyr::mutate(
    Tissue = dplyr::case_when(
      grepl("Adipose", colname, ignore.case = TRUE) ~ "Adipose",
      grepl("Blood", colname, ignore.case = TRUE) ~ "Blood",
      grepl("Muscle", colname, ignore.case = TRUE) ~ "Muscle",
      TRUE ~ NA_character_
    )
  ) %>%
  tibble::column_to_rownames("colname")

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

names(crosstissue_metab_samplemetadata)[4] <- "Modality"


# Converting the refmetmatrix into a prior knowledge format for running PLIER
refmetmatrix <- matrix(0L,nrow = length(rownames(combined_metab_joined)),ncol = dim(refmetmat)[2])
rownames(refmetmatrix) <- rownames(combined_metab_joined)
colnames(refmetmatrix) <- colnames(refmetmat)
for(i in 1:dim(refmetmatrix)[1]){
  ourmetab <- rownames(refmetmatrix)[i]
  if(ourmetab %in% rownames(refmetmat)){
    refmetmatrix[i,] <- refmetmat[ourmetab,]
  }
}

metab_finalpaths <- refmetmatrix

# Specifying the number of LVs we desire from PLIER. We typically use as input double the calculated number of PCs from the num.pc function
crosstissuemetabcombovalue <- min(num.pc(as.data.frame(combined_metab_joined), seed = 1), 50)
#chris: new implementation 72, so we'd get 144 LVs?
# previous dim of combined-metab_joined: 424 1473, 39 pcss -> new 455 1417, 72 pcs
# We cap at 50 (*2) LVs though, to keep things more interpretable.

# Running PLIER with default parameters
crosstissuemetabcombofin.plierResult.all.7=PLIER(as.matrix(combined_metab_joined),
                                                 metab_finalpaths,
                                                 k=2*crosstissuemetabcombovalue,trace=F, frac=0.7, scale=T, seed=1)

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
combined_metab_joined <- combined_metab_joined[,rownames(crosstissue_metab_samplemetadata)]


# Generating heatmaps of z-scored expression across all samples for top metabolites associated with each LV
# Heatmaps of select LVs of interest are part of Supplemental Figure S5
for(i in 1:dim(crosstissuemetabcombofin.plierResult.all.7$Z)[2]){

  top10rows <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,i]),][1:10,])
  datatoplot <- combined_metab_joined[top10rows,]

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
crosstissuemetablv15markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,15]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[11],])
crosstissuemetablv18markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,18]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[12],])
crosstissuemetablv21markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,21]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[17],])
crosstissuemetablv45markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,45]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[27],])
crosstissuemetablv86markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,86]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[38],])
crosstissuemetablv91markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,91]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[40],])
crosstissuemetablv97markers <- rownames(crosstissuemetabcombofin.plierResult.all.7$Z[order(-crosstissuemetabcombofin.plierResult.all.7$Z[,97]),][1:colSums(crosstissuemetabcombofin_Zmatscore > 2)[87],])

crosstissuemetablvmarkerlist <- list("LV15" = crosstissuemetablv15markers,
                                     "LV18" = crosstissuemetablv18markers,
                                     "LV21" = crosstissuemetablv21markers,
                                     "LV45" = crosstissuemetablv45markers,
                                     "LV86" = crosstissuemetablv86markers,
                                     "LV91" = crosstissuemetablv91markers,
                                     "LV97" = crosstissuemetablv97markers)

#####
# Now, we want to plot how similar/different the responses are in the three tissues
# We will get a score for relationship/correlation between each tissue
# First we need to divide the B matrix into three matrices for the tissues
# Then we need to rename the matrix columns to reflect pid, timepoint
####

crosstissuemetabbmat <- crosstissuemetabcombofin.plierResult.all.7$B
rownames(crosstissuemetabbmat) <- paste("LV",c(1:100),sep = "")
crosstissuemetabbmat_muscle <- crosstissuemetabbmat[,grep("muscle",colnames(crosstissuemetabbmat))]
colnames(crosstissuemetabbmat_muscle) = gsub("-muscle","", colnames(crosstissuemetabbmat_muscle))

crosstissuemetabbmat_adipose <- crosstissuemetabbmat[,grep("adipose",colnames(crosstissuemetabbmat))]
colnames(crosstissuemetabbmat_adipose) = gsub("-adipose","", colnames(crosstissuemetabbmat_adipose))

crosstissuemetabbmat_blood <- crosstissuemetabbmat[,grep("blood",colnames(crosstissuemetabbmat))]
colnames(crosstissuemetabbmat_blood) = gsub("-blood","", colnames(crosstissuemetabbmat_blood))

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


crosstissuemetablvmembershipbarplotdf <- data.frame("LV" = c(rep("LV15",length(crosstissuemetablvmarkerlist$LV15)),
                                                             rep("LV18",length(crosstissuemetablvmarkerlist$LV18)),
                                                             rep("LV21",length(crosstissuemetablvmarkerlist$LV21)),
                                                             rep("LV45",length(crosstissuemetablvmarkerlist$LV45)),
                                                             rep("LV86",length(crosstissuemetablvmarkerlist$LV86)),
                                                             rep("LV91",length(crosstissuemetablvmarkerlist$LV91)),
                                                             rep("LV97",length(crosstissuemetablvmarkerlist$LV97))))

crosstissuemetablvmembershipbarplotdf$LV <- factor(crosstissuemetablvmembershipbarplotdf$LV,levels = crosstissuemetabbmat_tissuecordf_forggplot$LV[order(crosstissuemetabbmat_tissuecordf$Blood_v_Muscle)])

# Figure 5G
png(file = "Figure5G_crosstissuemetab_lvtissuecorrelation_111725.png",width = 4,height = 4,units = "in",res = 600)
ggplot(crosstissuemetabbmat_tissuecordf_forggplot[crosstissuemetabbmat_tissuecordf_forggplot$LV %in% names(crosstissuemetablvmarkerlist),]) + geom_hline(yintercept = 0,size = 1,linetype = "dashed") + geom_point(aes(x=LV,y=Correlation,color=Comparison),size = 3) + theme_classic() + scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7")) + theme(axis.title.x=element_blank(),legend.position = "top",text = element_text(size=10),legend.key.spacing = unit(0,"cm")) + guides(col = guide_legend(ncol = 1))
dev.off()

# Figure 5H
png(file = "Figure5H_crosstissuemetab_lvmembershipbarplot_111725.png",width = 4,height = 3,units = "in",res = 600)
ggplot(crosstissuemetablvmembershipbarplotdf,aes(LV)) + geom_bar() + theme_classic() + ylab("Membership") + theme(axis.title.x=element_blank(),text = element_text(size=10))
dev.off()

crosstissuemetabcombofinUaucplottrim <- crosstissuemetabcombofin.plierResult.all.7$Uauc[,c(15,18,21,45,86,91,97)]
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
