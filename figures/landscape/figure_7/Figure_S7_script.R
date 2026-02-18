#comprehensive R script that generates all panels and tables for Figure S7

#load packages
library(MotrpacHumanPreSuspension)
library(MotrpacHumanPreSuspensionData)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(openxlsx)
library(tidyverse)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#### Figure S7A: Heatmap of TFEB ChIP-bound targets predicted by SC-ION ####

#this part get the table of TFEB ChIP-bound targets predicted by SC-ION
#list of TFEB targets in the SCION network
targets <- read.xlsx("TFEB_targets.xlsx")

#for all these sources, we will consider a peak if it is within 1kb of the TSS

#source #1 - Gambardella et al 2020
peaks <- read.xlsx("external_data/Gambardella-et-al-2020-table-s2.xlsx",sheet=2)

#search gene symbols. if there is a peak withink 1kb of the TSS, return Y
targets$Gambardella_et_al_2020 <- unlist(lapply(targets$SCION..target,function(x){
  my.peaks <- peaks %>% dplyr::filter(SYMBOL==x)
  if(sum(grepl("Promoter",my.peaks$annotation) & grepl("<=1",my.peaks$annotation)>0) ){
    return("Y")
  }else{
    return("N")
  }
}))

#source #2 - MSigDB-TFEB-TARGET-GENES
#the paper is Yevshin et al, 2018 where they reprocess multiple ChIP-seq datasets for many TFs
#this is already filtered for peaks within 1 kb of the TSS and is a list of gene symbols and gene IDs
#search gene symbols. if in the list, return Y
msigdb <- read_tsv("external_data/TFEB_TARGET_GENES.v2024.1.Hs.tsv",col_names=F)
peaks2 <- unlist(strsplit(msigdb$X2[18],",",1))
targets$MSigDb_TFEB_TARGET_GENES <- ifelse(targets$SCION..target%in%peaks2,"Y","N")

#source #3 - Settembre et al 2013
#ChIP qPCR where they show TFEB binds its own promoter in mice
#They also show TFEB binds PGC1A but PPARGC1A is not in the SCION network
#no table here since it's single target ChIP, so just create a column that indicates this
targets$Settembre_et_al_2013 = ifelse(targets$SCION..target=="TFEB","Y","N")

#source #4 - Doronzo et al, 2018
#ChIP-Seq on TFEB
#They only provided the BED file - we need to annotate it using annotatr package
#load the BED file and convert to GRanges
regions <- read_regions("external_data/GSM2354032_TFEB_e5_peaks.bed",genome="hg19")
#now annotate the regions
annotations <- build_annotations(genome="hg19",annotations="hg19_basicgenes")
regions_annotated <- annotate_regions(regions,annotations,ignore.strand=T,quiet=F)
regions_annotated_df <- data.frame(regions_annotated)

#for each symbol, check if there is a peak in promoter or 1to5kb upstream of the TSS.
targets$Doronzo_et_al_2018 <- unlist(lapply(targets$SCION..target,function(x){
  my.peaks <- regions_annotated_df %>% dplyr::filter(annot.symbol==x)
  if(sum(grepl("promoter",my.peaks$annot.type) | grepl("1to5kb",my.peaks$annot.type)>0) ){
    return("Y")
  }else{
    return("N")
  }
}))

#create a column that indicates any ChIP evidence
#this column is how we will filter the targets
targets$any <- ifelse(rowSums(targets=="Y")>0,"Y","N")

#save!
write_csv(targets,"TFEB_targets_with_peaks.csv")

#this part creates the heatmap
#get list of TFEB-bound targets
chip <- read_csv("TFEB_targets_with_peaks.csv")
chip_targets <- chip$SCION..target[chip$any=="Y"]

#get muscle DA
muscle_da <- MUSCLE_TRNSCRPT_DA
muscle_da <- muscle_da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")

#create two heatmaps of the RNAseq logFC
#first, use the the ChIP-bound genes (do not filter by gene set)
any_chip <- data.frame("gene_symbol"=chip_targets)
any_chip <- left_join(any_chip,muscle_da,by="gene_symbol") %>% dplyr::select("gene_symbol","contrast_short","logFC","p_value","adj_p_value")
any_chip_wide <- pivot_wider(any_chip, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
any_chip_wide <- any_chip_wide[,grepl("gene_symbol|delta-delta",colnames(any_chip_wide))]
rownames(any_chip_wide) <- any_chip_wide$gene_symbol

#separate out FC and adj p-value for the heatmap
fc <- any_chip_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- any_chip_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = any_chip_wide$gene_symbol

#set heatmap annotations
cdesc <- data.frame(group=c("ADUEndur","ADUEndur","ADUEndur","ADUResist","ADUResist","ADUResist"),
                    time=c("post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_15_30_45_min","post_3.5_4_hr","post_24_hr"))
ha <- HeatmapAnnotation(
  Modality = cdesc$group,
  Timepoint = cdesc$time,
  col = list(
    Modality = HUMAN_EXERCISE_GROUP_COLORS,
    Timepoint = HUMAN_ACUTE_TIMEPOINT_COLORS
  ))

#Heatmap dimensions
a <- unit(3, "mm")
width <- ncol(fc) * a
height <- nrow(fc) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(-1,0,1.5),
  colors = c("blue", "white", "red"))

#plot heatmap
ht <- Heatmap(as.matrix(fc),
              heatmap_legend_param = list(
                title = "logFC",
                legend_direction = "vertical",
                legend_width = unit(50, "mm")
              ),
              width = width,
              height = height,
              top_annotation = ha,
              #left_annotation = row_ha,
              cluster_columns = F,
              cluster_rows=T,
              show_column_names = F,
              show_row_names = F,
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "#555555"))
                
                #Point size
                if(pval[i, j] < 0.05) {
                  gb = textGrob("*")
                  gb_w = convertWidth(grobWidth(gb), "mm")
                  gb_h = convertHeight(grobHeight(gb), "mm")
                  grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                }
                
              },
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8)
)

#save to PDF
pdf("Fig_S7A_TFEB_ChIP_targets_heatmap.pdf",height=30,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#### Figure S7B - TFEB single feature plot ####
TFEB<-MotrpacHumanPreSuspension::plot_single_feature("TFEB",
                                                       selected_omes = "transcript-rna-seq")
ggsave(height=1.5,width=4.5,filename = "Fig_S7B_TFEB_multiplot_all_tps.pdf")

#### Figure S7C - HSP90AA1 single feature plot ####
HSP90AA1<-MotrpacHumanPreSuspension::plot_single_feature("HSP90AA1",
                                                     selected_omes = "transcript-rna-seq")
ggsave(height=1.5,width=4.5,filename = "Fig_S7C_HSP90AA1_multiplot_all_tps.pdf")

#### Figure S7D -PP2A RNA heatmap ####

#load DA and filter for PP2A subunits
da <- BLOOD_TRNSCRPT_DA
da <- da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
da <- left_join(da,HUMAN_FEATURE_TO_GENE,by="feature_id")

#all the annotated PP2A subunits
PP2a_subunits <- c("PPP2R1A","PPP2R1B","PPP2R2A","PPP2R2B","PPP2R2C","PPP2R2D","PPP2R3A","PPP2R3B","PPP2R3C","PPP2R4","PPP2R5A","PPP2R5B","PPP2R5C","PPP2R5D","PPP2R5E","PPP2CA","PPP2CB")
PP2A <- da %>% dplyr::filter(gene_symbol%in%PP2a_subunits)
PP2A_filter <- PP2A %>% dplyr::select("feature_id","gene_symbol","contrast_short","logFC","p_value","adj_p_value")
PP2A_wide <- pivot_wider(PP2A_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
PP2A_wide <- PP2A_wide[!grepl("PAR_Y",PP2A_wide$feature_id),grepl("feature_id|gene_symbol|delta-delta",colnames(PP2A_wide))]
rownames(PP2A_wide) <- PP2A_wide$gene_symbol

#repeat for muscle
da2 <- MUSCLE_TRNSCRPT_DA
da2 <- da2 %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
da2 <- left_join(da2,HUMAN_FEATURE_TO_GENE,by="feature_id")

PP2A2 <- da2 %>% dplyr::filter(gene_symbol%in%PP2a_subunits)
PP2A2_filter <- PP2A2 %>% dplyr::select("feature_id","gene_symbol","contrast_short","logFC","p_value","adj_p_value")
PP2A2_wide <- pivot_wider(PP2A2_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
PP2A2_wide <- PP2A2_wide[!grepl("PAR_Y",PP2A2_wide$feature_id),grepl("feature_id|gene_symbol|delta-delta",colnames(PP2A2_wide))]
rownames(PP2A2_wide) <- PP2A2_wide$gene_symbol

#combine blood and muscle into one heatmap for the figure
merge <- full_join(PP2A_wide,PP2A2_wide,by=c("feature_id","gene_symbol"),suffix=c(".blood",".muscle"))
rownames(merge) <- merge$gene_symbol
#remove the subunits not detected in both tissues
merge <- na.omit(merge)

#separate out FC and adj p-value for the heatmap
fc <- merge %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- merge %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = merge$gene_symbol

#set heatmap annotations
cdesc <- data.frame(tissue=c("blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","muscle","muscle","muscle","muscle","muscle","muscle"),
                    group=c("ADUEndur","ADUEndur","ADUEndur","ADUEndur","ADUEndur","ADUEndur","ADUResist","ADUResist","ADUResist","ADUResist","ADUEndur","ADUEndur","ADUEndur","ADUResist","ADUResist","ADUResist"),
                    time=c("during_20_min","during_40_min","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_15_30_45_min","post_3.5_4_hr","post_24_hr"))
ha <- HeatmapAnnotation(
  Tissue = cdesc$tissue,
  Modality = cdesc$group,
  Timepoint = cdesc$time,
  col = list(
    Tissue = HUMAN_TISSUE_COLORS,
    Modality = HUMAN_EXERCISE_GROUP_COLORS,
    Timepoint = HUMAN_ACUTE_TIMEPOINT_COLORS
  ))

#Heatmap dimensions
a <- unit(3, "mm")
width <- ncol(fc) * a
height <- nrow(fc) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(-1,0,1),
  colors = c("blue", "white", "red"))

#plot heatmap
ht <- Heatmap(fc,
              heatmap_legend_param = list(
                title = "logFC",
                legend_direction = "vertical",
                legend_width = unit(50, "mm")
              ),
              width = width,
              height = height,
              top_annotation = ha,
              #left_annotation = row_ha,
              cluster_columns = F,
              cluster_rows=T,
              show_column_names = F,
              show_row_names = T,
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "#555555"))
                
                #Point size
                if(pval[i, j] < 0.05) {
                  gb = textGrob("*")
                  gb_w = convertWidth(grobWidth(gb), "mm")
                  gb_h = convertHeight(grobHeight(gb), "mm")
                  grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                }
                
              },
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8)
)

#save to PDF
pdf("Fig_S7D_PP2A_RNA_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#### Figure S7E - key metabolite heatmap ####

#load DA 
da <- BLOOD_METAB_DA
da <- da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")

#metabolites of interest
metab_list <- c("Spermidine","N1-Acetylspermidine","Leucine","Isoleucine")
metab <- da %>% dplyr::filter(feature_id%in%metab_list)
metab <- metab[order(metab$feature_id),]
metab_filter <- metab %>% dplyr::select("feature_id","contrast_short","logFC","p_value","adj_p_value")
metab_wide <- pivot_wider(metab_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
metab_wide <- metab_wide[grepl("feature_id|feature_id|delta-delta",colnames(metab_wide))]
rownames(metab_wide) <- metab_wide$feature_id

#repeat for muscle
da2 <- MUSCLE_METAB_DA
da2 <- da2 %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")

metab2 <- da2 %>% dplyr::filter(feature_id%in%metab_list)
metab2 <- metab2[order(metab2$feature_id),]
metab2_filter <- metab2 %>% dplyr::select("feature_id","contrast_short","logFC","p_value","adj_p_value")
metab2_wide <- pivot_wider(metab2_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
metab2_wide <- metab2_wide[grepl("feature_id|delta-delta",colnames(metab2_wide))]
rownames(metab2_wide) <- metab2_wide$feature_id

#combine blood and muscle into one heatmap for the figure
merge <- full_join(metab_wide,metab2_wide,by=c("feature_id"),suffix=c(".blood",".muscle"))
rownames(merge) <- merge$feature_id
#remove the subunits not detected in both tissues
merge <- na.omit(merge)

#separate out FC and adj p-value for the heatmap
fc <- merge %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- merge %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = merge$feature_id

#set heatmap annotations
cdesc <- data.frame(tissue=c("blood","blood","blood","blood","blood","blood","blood","blood","blood","blood","muscle","muscle","muscle","muscle","muscle","muscle"),
                    group=c("ADUEndur","ADUEndur","ADUEndur","ADUEndur","ADUEndur","ADUEndur","ADUResist","ADUResist","ADUResist","ADUResist","ADUEndur","ADUEndur","ADUEndur","ADUResist","ADUResist","ADUResist"),
                    time=c("during_20_min","during_40_min","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_10_min","post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_15_30_45_min","post_3.5_4_hr","post_24_hr","post_15_30_45_min","post_3.5_4_hr","post_24_hr"))
ha <- HeatmapAnnotation(
  Tissue = cdesc$tissue,
  Modality = cdesc$group,
  Timepoint = cdesc$time,
  col = list(
    Tissue = HUMAN_TISSUE_COLORS,
    Modality = HUMAN_EXERCISE_GROUP_COLORS,
    Timepoint = HUMAN_ACUTE_TIMEPOINT_COLORS
  ))

#Heatmap dimensions
a <- unit(3, "mm")
width <- ncol(fc) * a
height <- nrow(fc) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(-1,0,1),
  colors = c("blue", "white", "red"))

#plot heatmap
ht <- Heatmap(fc,
              heatmap_legend_param = list(
                title = "logFC",
                legend_direction = "vertical",
                legend_width = unit(50, "mm")
              ),
              width = width,
              height = height,
              top_annotation = ha,
              #left_annotation = row_ha,
              cluster_columns = F,
              cluster_rows=F,
              show_column_names = F,
              show_row_names = T,
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "#555555"))
                
                #Point size
                if(pval[i, j] < 0.05) {
                  gb = textGrob("*")
                  gb_w = convertWidth(grobWidth(gb), "mm")
                  gb_h = convertHeight(grobHeight(gb), "mm")
                  grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                }
                
              },
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8)
)

#save to PDF
pdf("Fig_S7E_metab_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "left",
     legend_gap = unit(10, "mm"))
dev.off()

