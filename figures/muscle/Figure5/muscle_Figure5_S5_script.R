#comprehensive R script that generates panels for Figure 5/S5

#load packages
library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(clusterProfiler)
library(AnnotationDbi)
library(openxlsx)
library(purrr)
library(glue)
library(cmapR)
library(UpSetR)

# #code to prepare phospho data for PTM-SEA and reformat output files
# #load the muscle phospho DA 
# #filter for only the needed contrasts
# muscle_da <- MUSCLE_PROT_PH_DA %>% dplyr::filter(contrast_type=="exercise_with_controls" & contrast_category %in% c("EE-CON","RE-CON"))
# 
# #obtain list of confident sites from the phospho metadata
# conf_sites <- MUSCLE_PROT_PH_QC$feature_metadata %>% dplyr::filter(confident_site)
# muscle_da <- muscle_da %>% dplyr::filter(feature_id %in% conf_sites$id)
# 
# #need to make this a list format to work with run_PTMSEA
# DA_list <- list("muscle.prot-ph"=muscle_da)
# 
# run_PTMSEA(DA_list = DA_list,
#            repo_local_dir = "/Users/nclark/Desktop/MotrpacHumanPreSuspension",
#            selected_tissues = "muscle",
#            min_size = 5L)
# 
# #update the signature GCT files to incorporate metadata for the flanking sequences
# files <- list.files(path="muscle-signature_gct",pattern=".gct",full.names=T)
# for(i in files){
#   #read in the file
#   gct <- parse_gctx(i)
#   #obtain the flanking sequences
#   rid <- gct@rid
#   #for each flanking sequence, get all of the ids and the corresponding gene symbols
#   #if multiple, collapse using |
#   rdesc <- lapply(rid, function(j){
#     meta <- conf_sites[grepl(j,conf_sites$flanking_sequence),] %>% dplyr::select("ptm_id","gene_symbol")
#     meta2 <- data.frame("ptm_id"=paste(meta$ptm_id,collapse="|"),
#                         "gene_symbol"=paste(meta$gene_symbol,collapse="|"))
#     return(meta2)
#   })
#   rdesc <- Reduce(rbind,rdesc)
#   rdesc$id <- rid
#   #add to current rdesc and re-save the gct file
#   gct@rdesc <- left_join(gct@rdesc,rdesc,by="id")
#   gct@rid <- rid
#   fn <- gsub("muscle-signature_gct","muscle-signature_gct_fixed",i)
#   write_gct(gct,fn,appenddim=F)
# }

#Figure 5A and S5A, bubble heatmaps and upset plot
#PTMSEA results
ptmsea <- parse_gctx("ptm-sea-results-combined.gct")
rdesc <- ptmsea@rdesc
mat <- ptmsea@mat

#create upset plot
upset <- list("EE post 15 min"=rownames(rdesc)[rdesc$fdr.pvalue.muscle.z.std_Endur.post_15_30_45_min...Control.post_15_30_45_min..delta.delta.<=0.05],
              "EE post 3.5 hr"=rownames(rdesc)[rdesc$fdr.pvalue.muscle.z.std_Endur.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.<=0.05],
              "EE post 24 hr"=rownames(rdesc)[rdesc$fdr.pvalue.muscle.z.std_Endur.post_24_hr...Control.post_24_hr..delta.delta.<=0.05],
              "RE post 15 min"=rownames(rdesc)[rdesc$fdr.pvalue.muscle.z.std_Resist.post_15_30_45_min...Control.post_15_30_45_min..delta.delta.<=0.05],
              "RE post 3.5 hr"=rownames(rdesc)[rdesc$fdr.pvalue.muscle.z.std_Resist.post_3.5_4_hr...Control.post_3.5_4_hr..delta.delta.<=0.05],
              "RE post 24 hr"=rownames(rdesc)[rdesc$fdr.pvalue.muscle.z.std_Resist.post_24_hr...Control.post_24_hr..delta.delta.<=0.05])

#create an overlap table to save
binary_matrix <- data.frame(
  row.names = unique(unlist(upset)),
  EE_15min = as.numeric(unique(unlist(upset)) %in% upset$`EE post 15 min`),
  EE_3.5hr = as.numeric(unique(unlist(upset)) %in% upset$`EE post 3.5 hr`),
  EE_24hr = as.numeric(unique(unlist(upset)) %in% upset$`EE post 24 hr`),
  RE_15min = as.numeric(unique(unlist(upset)) %in% upset$`RE post 15 min`),
  RE_3.5hr = as.numeric(unique(unlist(upset)) %in% upset$`RE post 3.5 hr`),
  RE_24hr = as.numeric(unique(unlist(upset)) %in% upset$`RE post 24 hr`)
)
write.csv(binary_matrix,"PTMSEA-overlap.csv")

#create upset plot
pdf("FigS5A_PTMSEA-muscle-upset-plot.pdf",height=6,width=10)
upset(fromList(upset),nsets=6,nintersects=NA,keep.order=T,order.by="freq",decreasing=T,text.scale=2)
dev.off()

#for a list of terms, produce a leading edge heatmap of the genes in those terms
ee_terms <- c(
  "KINASE-PSP_Src/SRC",
  "KINASE-PSP_PKG2/PRKG2",
  "KINASE-PSP_Chk1/CHEK1",
  "KINASE-PSP_BRAF",
  "KINASE-PSP_AMPKA1/PRKAA1",
  "KINASE-iKiP_TSSK2",
  "KINASE-iKiP_TAOK2",
  "KINASE-iKiP_SIK2",
  "KINASE-iKiP_SIK1",
  "KINASE-iKiP_SGK3",
  "KINASE-iKiP_RPS6KB1",
  "KINASE-iKiP_RPS6KA6",
  "KINASE-iKiP_MAST2",
  "KINASE-iKiP_EEF2K",
  "KINASE-iKiP_DCLK2",
  "KINASE-iKiP_CHUK.IKKA",
  "KINASE-iKiP_CGK2",
  "KINASE-iKiP_CDK9-CCNK",
  "KINASE-iKiP_CAMK2G",
  "KINASE-iKiP_CAMK2D",
  "KINASE-iKiP_CAMK2B",
  "KINASE-iKiP_CAMK2A",
  "KINASE-iKiP_AURKB",
  "KINASE-iKiP_AKT1")

#separate out z-score and adj p-value for the heatmap
z <- mat[row.names(mat)%in%ee_terms,]
pval <- rdesc %>% dplyr::select(contains("fdr")) %>% data.frame()
pval <- pval[row.names(pval)%in%ee_terms,]

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
a <- unit(5, "mm")
width <- ncol(z) * a
height <- nrow(z) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(min(z),0,max(z)),
  colors = c("#4266F6","white", "#953C2A"))

#plot heatmap
ht <- Heatmap(as.matrix(z),
              heatmap_legend_param = list(
                title = "NES",
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
              column_names_side = "top",
              show_row_names = T,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              column_split = cdesc$group,
              column_gap = unit(2,"mm"),
              border_gp = gpar(col = "black"),
              column_title=NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "gray",
                                    fill = ifelse(pval[i,j] <= 0.05,"gray","white")))
                #scale pvalues so the circles don't get so big
                norm_per <- ifelse(z[i,j]>10,10,z[i,j])
                #make percentage 0 to 1 instead of 0 to 100
                #norm_per <- norm_per/100
                
                #Point size
                grid.circle(x = x,
                            y = y,
                            r = norm_per*0.025,
                            gp = gpar(fill = col_fun(z[i, j]),
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
pdf("FigS5A_PTMSEA_EE_only_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_list = lgd_list,
     legend_gap = unit(10, "mm"))
dev.off()

re_terms <- c(
  "KINASE-PSP_Ret/RET",
  "KINASE-PSP_MARK2",
  "KINASE-PSP_MARK1",
  "KINASE-PSP_ERK1/MAPK3",
  "KINASE-PSP_ATM",
  "KINASE-iKiP_TTBK1",
  "KINASE-iKiP_TGFBR1",
  "KINASE-iKiP_SPHK2",
  "KINASE-iKiP_RIPK5",
  "KINASE-iKiP_MAPK11",
  "KINASE-iKiP_GRK7",
  "KINASE-iKiP_GRK2",
  "KINASE-iKiP_BMPR1B",
  "KINASE-iKiP_ALK",
  "KINASE-iKiP_AKT2",
  "KINASE-PSP_P38A/MAPK14",
  "KINASE-iKiP_STK11",
  "KINASE-iKiP_MAPK9.JNK2",
  "KINASE-iKiP_MAPK8.JNK1",
  "KINASE-iKiP_MAPK14",
  "KINASE-iKiP_MAPK13",
  "KINASE-iKiP_MAPK10.JNK3",
  "KINASE-iKiP_GRK5",
  "KINASE-iKiP_EIF2AK1")

#separate out z-score and adj p-value for the heatmap
z <- mat[row.names(mat)%in%re_terms,]
pval <- rdesc %>% dplyr::select(contains("fdr")) %>% data.frame()
pval <- pval[row.names(pval)%in%re_terms,]

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
a <- unit(5, "mm")
width <- ncol(z) * a
height <- nrow(z) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(min(z),0,max(z)),
  colors = c("#4266F6","white", "#953C2A"))

#plot heatmap
ht <- Heatmap(as.matrix(z),
              heatmap_legend_param = list(
                title = "NES",
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
              column_names_side = "top",
              show_row_names = T,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              column_split = cdesc$group,
              column_gap = unit(2,"mm"),
              border_gp = gpar(col = "black"),
              column_title=NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "gray",
                                    fill = ifelse(pval[i,j] <= 0.05,"gray","white")))
                #scale pvalues so the circles don't get so big
                norm_per <- ifelse(z[i,j]>10,10,z[i,j])
                #make percentage 0 to 1 instead of 0 to 100
                #norm_per <- norm_per/100
                
                #Point size
                grid.circle(x = x,
                            y = y,
                            r = norm_per*0.02,
                            gp = gpar(fill = col_fun(z[i, j]),
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
pdf("FigS5A_PTMSEA_RE_only_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_list = lgd_list,
     legend_gap = unit(10, "mm"))
dev.off()

fig4a <- c(
  "PERT-PSP_EGF",
  "PATH-BI_ISCHEMIA",
  "KINASE-PSP_MAPKAPK2",
  "KINASE-iKiP_MAPKAPK3",
  "PERT-PSP_TNF",
  "PERT-PSP_IL_1B",
  "PERT-PSP_IGF_1",
  "PATH-NP_IL11_PATHWAY",
  "KINASE-PSP_P38B/MAPK11",
  "PERT-PSP_THROMBIN",
  "PERT-PSP_IL_2",
  "PATH-NP_IL6_PATHWAY",
  "PERT-PSP_EXERCISE",
  "KINASE-PSP_P38A/MAPK14",
  "KINASE-PSP_ERK1/MAPK3",
  "KINASE-iKiP_CAMK2G",
  "KINASE-iKiP_CAMK2D",
  "KINASE-iKiP_CAMK2A",
  "KINASE-PSP_AMPKA1/PRKAA1",
  "KINASE-iKiP_CAMK2B",
  "KINASE-iKiP_AKT2",
  "KINASE-iKiP_AKT1",
  "KINASE-iKiP_AMPKA2",
  "KINASE-iKiP_AMPKA1",
  "KINASE-iKiP_MAPKAPK5",
  "KINASE-iKiP_MAPKAPK2",
  "KINASE-iKiP_MAPK3.ERK1",
  "KINASE-iKiP_MAPK1.ERK2",
  "KINASE-iKiP_MAPK9.JNK2",
  "KINASE-iKiP_MAPK8.JNK1",
  "KINASE-iKiP_MAPK10.JNK3",
  "KINASE-iKiP_MAPK13",
  "KINASE-iKiP_MAPK12",
  "KINASE-iKiP_MAPK14",
  "KINASE-iKiP_MAPK11",
  "KINASE-PSP_mTOR/MTOR",
  "KINASE-iKiP_TGFBR1",
  "KINASE-iKiP_GRK7",
  "KINASE-iKiP_GRK2",
  "KINASE-iKiP_HIPK3",
  "KINASE-iKiP_HIPK2",
  "KINASE-iKiP_GRK5",
  "KINASE-iKiP_EIF2AK1",
  "KINASE-iKiP_CSNK1G2",
  "KINASE-iKiP_CSNK1G3",
  "KINASE-iKiP_CSNK2A2.CK2A2",
  "KINASE-iKiP_CSNK2A1.CK2A1")

#separate out z-score and adj p-value for the heatmap
z <- mat[row.names(mat)%in%fig4a,]
pval <- rdesc %>% dplyr::select(contains("fdr")) %>% data.frame()
pval <- pval[row.names(pval)%in%fig4a,]

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
a <- unit(5, "mm")
width <- ncol(z) * a
height <- nrow(z) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(min(z),0,max(z)),
  colors = c("#4266F6","white", "#953C2A"))

#plot heatmap
ht <- Heatmap(as.matrix(z),
              heatmap_legend_param = list(
                title = "NES",
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
              column_names_side = "top",
              show_row_names = T,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              column_split = cdesc$group,
              column_gap = unit(2,"mm"),
              border_gp = gpar(col = "black"),
              column_title=NULL,
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #Heatmap grid
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "gray",
                                    fill = ifelse(pval[i,j] <= 0.05,"gray","white")))
                #scale pvalues so the circles don't get so big
                norm_per <- ifelse(z[i,j]>9,9,z[i,j])
                #make percentage 0 to 1 instead of 0 to 100
                #norm_per <- norm_per/100
                
                #Point size
                grid.circle(x = x,
                            y = y,
                            r = norm_per*0.015,
                            gp = gpar(fill = col_fun(z[i, j]),
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
pdf("Fig5A_PTMSEA_heatmap.pdf",height=12,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_list = lgd_list,
     legend_gap = unit(10, "mm"))
dev.off()

#Figure 5B, HSPB1 heatmap
#load DA for all three omes, merge, then add gene symbols
da <- rbind(MUSCLE_TRNSCRPT_DA,MUSCLE_PROT_PR_DA,MUSCLE_PROT_PH_DA)
da <- da %>% dplyr::filter(contrast_type=="exercise_with_controls")
da <- left_join(da,HUMAN_FEATURE_TO_GENE,by="feature_id")

HSPB1 <- da %>% dplyr::filter(gene_symbol=="HSPB1")
HSPB1_filter <- HSPB1 %>% dplyr::select("feature_id","gene_symbol","contrast_short","logFC","p_value","adj_p_value")
HSPB1_wide <- pivot_wider(HSPB1_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) %>% data.frame()
rownames(HSPB1_wide) <- HSPB1_wide$feature_id

#separate out FC and adj p-value for the heatmap
fc <- HSPB1_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- HSPB1_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = paste(HSPB1_wide$gene_symbol,gsub(".*_", "", HSPB1_wide$feature_id),sep="_")

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
  breaks = c(-1,0,2),
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
                if(pval[i, j] <= 0.05) {
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
pdf("Fig5B_HSPB1_single_feature_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#Figure 5D, HIPK3 signature heatmap
#create heatmap for HIPK3 signature
sig <- parse_gctx("KINASE.iKiP_HIPK3_n6x22.gct")
mat <- sig@mat
rdesc <- sig@rdesc

#separate out NES and adj p-value for the heatmap
nes <- mat
pval <- rdesc %>% dplyr::select(contains("adj_p_value")) %>% data.frame()
rownames(nes) = rownames(pval) = paste(rdesc$gene_symbol,gsub(".*_", "", rdesc$feature_id),sep="_")

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
width <- ncol(nes) * a
height <- nrow(nes) * a

#color scale
col_fun <- circlize::colorRamp2(
  breaks = c(min(nes),0,max(nes)),
  colors = c("blue", "white", "red"))

#plot heatmap
ht <- Heatmap(as.matrix(nes),
              heatmap_legend_param = list(
                title = "NES",
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
                if(pval[i, j] <= 0.05) {
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
pdf("Fig5D_KINASE.iKiP_HIPK3_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#Figure 5E, Figure S5D and S5E, ORA on HIPK3 nearest neighbor analysis results
#nearest neighbor sites
sites <- read_csv("HIPK3_Y359_NearestNeighbor_Pearson_0p8.csv")

#get the DA information so we can filter for DA sites only before performing ORA
muscle_da <- MUSCLE_PROT_PH_DA %>% dplyr::filter(contrast_type=="exercise_with_controls")
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")
background <- muscle_da

#filter for only the significant sites for the nearest neighbor analysis
muscle_sig <- muscle_da %>% dplyr::filter(adj_p_value <= 0.05)
sites_da <- sites[sites$feature_id%in%muscle_sig$feature_id,]
write_csv(sites_da,"HIPK3_Y359_NearestNeighbor_Pearson_0p8_DA.csv")

#run ORA
input_list <- list("HIPK3"=sites_da)
ora.results.list <- lapply(names(input_list),function(x){
  message(glue("Performing enrichment test for gene list {x}"))
  gene.sig <- input_list[[x]] %>% dplyr::select(gene_symbol)  #get features in this cluster
  gene.sig <- gene.sig[gene.sig$gene_symbol%in%background$gene_symbol,] 
  #C2, C5, MITOCARTA separately
  C2 <- run_ORA(gene.sig$gene_symbol,as.character(background$gene_symbol),c("REACTOME","WP","PID","BIOCARTA","KEGG_MEDICUS"))
  C5 <- run_ORA(gene.sig$gene_symbol,as.character(background$gene_symbol),c("GOBP","GOMF","GOCC"))
  mito <- run_ORA(gene.sig$gene_symbol,as.character(background$gene_symbol),"MITOCARTA")
  return(list("C2"=C2,"C5"=C5,"MITOCARTA"=mito))
})

#save all results
ora.results <- ora.results.list[[1]]
write.xlsx(ora.results,"ORA_HIPK3_nearest_neighbor_0p8_DA.xlsx")

#terms to use for leading edge analysis and bar plot
sig_terms <- c("WP_STRIATED_MUSCLE_CONTRACTION_PATHWAY",
               "GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
               "GOBP_SARCOMERE_ORGANIZATION",
               "GOBP_DETECTION_OF_MUSCLE_STRETCH",
               "GOMF_MUSCLE_ALPHA_ACTININ_BINDING"
)

#bar plot of the significant terms
ora_sig <- rbind(ora.results$C2[ora.results$C2$set_short%in%sig_terms,],
                 ora.results$C5[ora.results$C5$set_short%in%sig_terms,])
go_data <- data.frame(term=ora_sig$set_short,
                      description=ora_sig$set_short,
                      p_value=ora_sig$p_value)

# Calculate -log10(p-value) for plotting
go_data$neg_log10_p <- -log10(go_data$p_value)

# Create the horizontal bar chart
pdf(glue("Fig5E_HIPK3-ORA-bar-chart.pdf"),height=3,width=10)
ggplot(go_data, aes(x = neg_log10_p, y = reorder(description, neg_log10_p))) +
  geom_col(aes(fill = neg_log10_p),color = "black", linewidth=0.5) +
  scale_fill_gradient(low = "#f0e6ff", high = "#57347d", 
                      name = "-log10(p-value)") +
  labs(
    x = "-log10(p-value)",
    y = "Term"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
dev.off()

sig_terms <- c("GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE",
               "GOBP_SARCOMERE_ORGANIZATION"
)

for(i in sig_terms){
  #get genes in this term
  if(grepl("_",i)){
    db <- str_split_i(i,"_",1)
    sig_list <- MOLECULAR_SIGNATURES[[db]][[i]]
  }else{
    sig_list <- MOLECULAR_SIGNATURES[["MITOCARTA"]][[paste("MITOCARTA",i,sep="_")]]
  }
  my_sites <- sites_da[sites_da$gene_symbol%in%sig_list,]
  
  #separate out FC and adj p-value for the heatmap
  fc <- my_sites %>% dplyr::select(contains("logFC")) %>% data.frame()
  pval <- my_sites %>% dplyr::select(contains("adj_p_value")) %>% data.frame()
  rownames(fc) = rownames(pval) = paste(my_sites$gene_symbol,gsub(".*_", "", my_sites$feature_id),sep="_")
  
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
    breaks = c(-1,0,1),
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
                  if(pval[i, j] <= 0.05) {
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
  if(i=="GOMF_STRUCTURAL_CONSTITUENT_OF_MUSCLE"){
    pdf(glue("FigS5E_{i}.pdf"),height=8,width=8)
  }else{
    pdf(glue("FigS5D_{i}.pdf"),height=8,width=8)
  }
  draw(ht,
       merge_legends = TRUE,
       heatmap_legend_side = "right",
       legend_gap = unit(10, "mm"))
  dev.off()
  
}

#Figure S5F, LMOD2 heatmap

#filter for only significant LMOD2 sites
LMOD2 <- da %>% dplyr::filter(gene_symbol=="LMOD2")
LMOD2_filter <- LMOD2 %>% dplyr::select("feature_id","gene_symbol","contrast_short","logFC","p_value","adj_p_value")
LMOD2_wide <- pivot_wider(LMOD2_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) %>% data.frame()
pval <- LMOD2_wide %>% dplyr::select(contains("adj_p_value"))  %>% as.matrix()
LMOD2_wide <- LMOD2_wide[rowSums(pval <= 0.05)>0,]
rownames(LMOD2_wide) <- LMOD2_wide$feature_id

#separate out FC and adj p-value for the heatmap
fc <- LMOD2_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- LMOD2_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = paste(LMOD2_wide$gene_symbol,gsub(".*_", "", LMOD2_wide$feature_id),sep="_")

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
  breaks = c(-2,0,2),
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
                if(pval[i, j] <= 0.05) {
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
pdf("FigS5F_LMOD2_single_feature_heatmap_sig_only.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

