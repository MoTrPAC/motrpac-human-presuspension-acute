#comprehensive R script that generates all panels and tables for Figure 7

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

#### Figure 7C - TFEB Phosphosite Heatmap ####

#load DA and filter for TFEB phosphosites
#note for public: this is equivalent to:
#`load_differential_analysis(selected_omes = "prot-ph", selected_tissues = "muscle", single_matrix = TRUE)`
#we recommend using `load_differential_analysis` unless you're very familiar with the data objects.

da <- MUSCLE_PROT_PH_DA
da <- da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
da <- left_join(da,HUMAN_FEATURE_TO_GENE,by="feature_id")

TFEB <- da %>% dplyr::filter(gene_symbol=="TFEB")
TFEB_filter <- TFEB %>% dplyr::select("feature_id","contrast_short","logFC","p_value","adj_p_value")
TFEB_wide <- pivot_wider(TFEB_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value"))
TFEB_wide <- TFEB_wide[,grepl("feature_id|delta-delta",colnames(TFEB_wide))]
sites <- unlist(strsplit(as.character(TFEB_wide$feature_id),"_"))
sites <- sites[seq(2,length(sites),by=2)]
rownames(TFEB_wide) <- paste("TFEB",sites,sep=".")

#separate out FC and adj p-value for the heatmap
fc <- TFEB_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- TFEB_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = TFEB_wide$feature_id

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
  breaks = c(min(fc),0,1),
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
pdf("Fig_7C_TFEB_phospho_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#### Figure 7D - TFEB ChIP target ORA ####

#This chunk of the code generates the TFEB ChIP targets table that is part of Datset S7
#list of all transcripts detected in muscle
muscle_trans <- MUSCLE_TRNSCRPT_DA %>% dplyr::select(feature_id) %>% unique()
#filter for just those with valid gene symbols
muscle_trans <- left_join(muscle_trans,HUMAN_FEATURE_TO_GENE,by="feature_id") %>%
  dplyr::select("feature_id":"ensembl_gene") %>%
  dplyr::filter(!is.na(gene_symbol))

#for all these sources, we will consider a peak if it is within 1kb of the TSS

#source #1 - Gambardella et al 2020
peaks <- read.xlsx("external_data/Gambardella-et-al-2020-table-s2.xlsx",sheet=2)

#search gene symbols. if there is a peak withink 1kb of the TSS, return Y
muscle_trans$Gambardella_et_al_2020 <- unlist(lapply(muscle_trans$gene_symbol,function(x){
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
muscle_trans$MSigDb_TFEB_TARGET_GENES <- ifelse(muscle_trans$gene_symbol%in%peaks2,"Y","N")

#source #3 - Settembre et al 2013
#ChIP qPCR where they show TFEB binds its own promoter in mice
#They also show TFEB binds PGC1A but PPARGC1A is not in the SCION network
#no table here since it's single target ChIP, so just create a column that indicates this
muscle_trans$Settembre_et_al_2013 = ifelse(muscle_trans$gene_symbol=="TFEB" | muscle_trans$gene_symbol=="PPARGC1A","Y","N")

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
muscle_trans$Doronzo_et_al_2018 <- unlist(lapply(muscle_trans$gene_symbol,function(x){
  my.peaks <- regions_annotated_df %>% dplyr::filter(annot.symbol==x)
  if(sum(grepl("promoter",my.peaks$annot.type) | grepl("1to5kb",my.peaks$annot.type)>0) ){
    return("Y")
  }else{
    return("N")
  }
}))

#create a column that indicates any ChIP evidence
#this column is how we will filter the targets
chip <- muscle_trans %>% dplyr::select("Gambardella_et_al_2020":"Doronzo_et_al_2018")
muscle_trans$any <- ifelse(rowSums(chip=="Y")>0,"Y","N")

#save!
write_csv(muscle_trans,"TFEB_ChIP_muscle.csv")

#This chunk of the code performs the ORA
#load muscle DA to obtain background
muscle_da <- MUSCLE_TRNSCRPT_DA
background <- data.frame("feature_id"=as.character(unique(muscle_da$feature_id)))
background <- right_join(HUMAN_FEATURE_TO_GENE,background,by="feature_id")

#filter muscle DA into the four lists of interest
muscle_da <- muscle_da %>% dplyr::filter(contrast_type=="exercise_with_controls" & contrast_category %in% c("EE-CON","RE-CON") & Timepoint=="post_3.5_4_hr")
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")
muscle_ee <- muscle_da %>% dplyr::filter(contrast_category=="EE-CON" & adj_p_value < 0.05)
muscle_re <- muscle_da %>% dplyr::filter(contrast_category=="RE-CON" & adj_p_value < 0.05)
muscle_ee_up <- muscle_ee %>% dplyr::filter(logFC>0)
muscle_ee_down <- muscle_ee %>% dplyr::filter(logFC<0)
muscle_re_up <- muscle_re %>% dplyr::filter(logFC>0)
muscle_re_down <- muscle_re %>% dplyr::filter(logFC<0)

#filter these lists for TFEB ChIP targets
chip <- read_csv("TFEB_ChIP_muscle.csv") %>% dplyr::filter(any=="Y")

input_list <- list("EE_3.5_4_hr_up" = as.character(muscle_ee_up$gene_symbol[muscle_ee_up$feature_id%in%chip$feature_id]),
                   "EE_3.5_4_hr_down" = as.character(muscle_ee_down$gene_symbol[muscle_ee_down$feature_id%in%chip$feature_id]),
                   "RE_3.5_4_hr_up" = as.character(muscle_re_up$gene_symbol[muscle_re_up$feature_id%in%chip$feature_id]),
                   "RE_3.5_4_hr_down" = as.character(muscle_re_down$gene_symbol[muscle_re_down$feature_id%in%chip$feature_id]))

#run on each gene set
ora.results.list <- lapply(input_list,function(x){
  enrich <- run_ORA(x,as.character(background$gene_symbol))
  return(enrich)
})
#save all results
names(ora.results.list) <- names(input_list)
write.xlsx(ora.results.list,"TFEB_ChIP_targets_DA_3.5_4_hr_ORA_oct2025.xlsx")

#create heatmap of key terms
terms <- read.xlsx("TFEB_ChIP_targets_DA_3.5_4_hr_ORA_oct2025_highlighted.xlsx")
#filter the results for just these terms
ora.results.filt <- lapply(ora.results.list,function(x){
  return(x[x$set%in%terms$set,])
})

#get data for heatmap
ora.results.heatmap <- lapply(names(ora.results.filt),function(x){
  table <- ora.results.filt[[x]] %>% dplyr::select("set_short","set_size_in_input","input_size","adj_p_value")
  table$input_ratio <- table$set_size_in_input/table$input_size
  colnames(table)[colnames(table)!="set_short"] <- paste(colnames(table)[colnames(table)!="set_short"],x,sep=".")
  return(table)
})
names(ora.results.heatmap) <- names(ora.results.filt)
ee_up <- ora.results.heatmap[[1]]
ee_dn <- ora.results.heatmap[[2]]
re_up <- ora.results.heatmap[[3]]
re_dn <- ora.results.heatmap[[4]]
data_heatmap <- full_join(ee_up,ee_dn,by="set_short")
data_heatmap <- full_join(data_heatmap,re_up,by="set_short")
data_heatmap <- full_join(data_heatmap,re_dn,by="set_short")
rownames(data_heatmap) <- data_heatmap$set_short

#separate out adjusted p-value and input_ratio for heatmap
#ratio_heatmap <- data_heatmap %>% dplyr::select(contains("input_ratio"))
pval_heatmap <- data_heatmap %>% dplyr::select(contains("adj_p_value"))

#convert to -log10(pval)
pval_heatmap <- -log10(pval_heatmap)
colnames(pval_heatmap) <- gsub("adj_p_value.","",colnames(pval_heatmap))

#Heatmap dimensions
a <- unit(5, "mm")
width <- ncol(pval_heatmap) * a
height <- nrow(pval_heatmap) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(0,max(pval_heatmap)),
  colors = c("white", "#57347d"))

#plot heatmap
ht <- Heatmap(as.matrix(pval_heatmap),
              heatmap_legend_param = list(
                title = "-log10(adjusted p-value)",
                legend_direction = "vertical",
                legend_width = unit(50, "mm"),
                at = c(0,round(max(pval_heatmap),2))
              ),
              width = width,
              height = height,
              #top_annotation = ha,
              #left_annotation = row_ha,
              cluster_columns = F,
              cluster_rows=T,
              show_column_names = T,
              column_names_side = "top",
              show_row_names = T,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              cell_fun = function(j, i, x, y, width, height, fill) {

                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "#555555",
                                    fill = ifelse(pval_heatmap[i,j] > -log10(0.05),"gray","white")))
                #scale pvalues so the circles don't get so big
                norm_per <- ifelse(pval_heatmap[i,j]>5,5,pval_heatmap[i,j])
                #make percentage 0 to 1 instead of 0 to 100
                #norm_per <- norm_per/100

                #Point size
                grid.circle(x = x,
                            y = y,
                            r = norm_per*0.02,
                            gp = gpar(fill = col_fun(pval_heatmap[i, j]),
                                      col = "black"))

              },
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8)
)

# lgd_list = list(
#   Legend(labels = c("<1","5","10","15","20"),
#          title = "% associated genes",
#          type = "points",
#          # ncol = 5,
#          pch = 16,
#          size = unit(1:5, 'mm'),
#          legend_gp = gpar(
#            col = rep("black",4)
#          ),
#          background = 'white'
#   )
# )

lgd_list = list(
  Legend(labels = c("< 0.05",">= 0.05"),
         title = "BH Adjusted p-value",
         type = "grid",
         # ncol = 5,
         #pch = 16,
         #size = unit(1:5, 'mm'),
         legend_gp = gpar(
           fill=c("gray","white")
         ),
         border = 'black'
  )
)

#save to PDF
pdf("Fig_7D_TFEB_ChIP_ORA_heatmap.pdf",height=16,width=10)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "left",
     annotation_legend_list = lgd_list,
     legend_gap = unit(10, "mm"))
dev.off()

#### Figure 7E - TFEB/RXRA/CHCHD3 shared targets Heatmap ####

#get da results
muscle_da <- MUSCLE_TRNSCRPT_DA %>% dplyr::filter(contrast_type=="exercise_with_controls" & contrast_category %in% c("EE-CON","RE-CON"))
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")

#list of genes for heatmap
genes <- data.frame("gene_symbol"=c("RXRA","CHAMP1","CEBPB","TFEB","NFYA","PLA2G4E","SPTLC1","CIRBP","TUBGCP6","C1QTNF6"))

#add the DA information for these genes
genes <- left_join(genes,muscle_da,by="gene_symbol") %>% dplyr::select("gene_symbol","contrast_short","logFC","p_value","adj_p_value")
genes_wide <- pivot_wider(genes, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value"))
rownames(genes_wide) <- genes_wide$gene_symbol

#separate out FC and adj p-value for the heatmap
fc <- genes_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- genes_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = genes_wide$gene_symbol

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
pdf("Fig_7E_TFEB_RXRA_CHCHD3_shared_targets_heatmap.pdf",height=6,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#### Figure 7F - Lactic acid single feature plot ####
lactic<-MotrpacHumanPreSuspension::plot_single_feature("Lactic acid",
                                                         selected_omes = "metab")
ggsave(height=1.5,width=4.5,filename = "Fig_7F_lactic_acid_multiplot_all_tps.pdf")
