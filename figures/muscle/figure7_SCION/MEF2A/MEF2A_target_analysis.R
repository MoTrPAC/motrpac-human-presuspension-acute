#MEF2A target analysis

library(MotrpacHumanPreSuspension)
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

#load MEF2A targets
mef2a <- read_csv("MEF2A_targets.csv") %>% dplyr::select("name") %>% dplyr::filter(name!="MEF2A") %>% dplyr::rename(GeneSym=name)

#ENCODE resource
encode <- read_csv("MEF2A_encode.csv") %>% dplyr::rename("ENCODE_2015"="MEF2A")
mef2a <- left_join(mef2a,encode,by="GeneSym")

#mouse ChIP
mouse_chip <- read.xlsx("41467_2019_12812_MOESM4_ESM.xlsx") %>% dplyr::filter(TF=="Mef2a")
#need to convert mouse to human gene symbols
mapIt <- function(x) {
  require("Orthology.eg.db", character.only = TRUE)
  require("org.Mm.eg.db", character.only = TRUE)
  require("org.Hs.eg.db", character.only = TRUE)
  mouse <- mapIds(org.Mm.eg.db, x, "ENTREZID", "REFSEQ")
  mapped <- select(Orthology.eg.db, mouse,"Homo.sapiens","Mus.musculus")
  human <- mapIds(org.Hs.eg.db, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  human <- do.call(c, lapply(human, function(x) if(is.null(x)) return(NA) else return(x)))
  cbind(x, mapped, human)
}
symbols <- mapIt(mouse_chip$Nearest.Refseq) %>% dplyr::rename(GeneSym=human,Nearest.Refseq=x)
mouse_chip <- left_join(mouse_chip,symbols,by="Nearest.Refseq",relationship="many-to-many")
mouse_chip <- mouse_chip[!is.na(mouse_chip$GeneSym),]

#for each gene, check if peak within 1kb of TSS. if yes, mark 1.
mef2a$Akerberg_et_al_2019 <- unlist(lapply(mef2a$GeneSym,function(x){
  peaks <- mouse_chip[mouse_chip$GeneSym==x,]
  ifelse(sum(peaks$Distance.to.TSS<0 & peaks$Distance.to.TSS>=-1000)>0,1,0)
}))
mef2a$any <- ifelse(rowSums(mef2a[,-1],na.rm=T)>0,1,0)

#save!
write_csv(mef2a,"MEF2A_targets_with_peaks.csv")

#perform ORA
#load list of TFEB targets in each network (EE or RE)
ee <- read_csv("MEF2A_direct_targets_EE.csv")
re <- read_csv("MEF2A_direct_targets_RE.csv")
either <- c(ee$name,re$name) #names in either EE OR RE
both <- either[duplicated(either)] #names in both EE and RE

#get background list which is all detected transcripts that went into the DA
muscle_da <- MUSCLE_TRNSCRPT_DA
muscle_da <- muscle_da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")
background <- as.character(unique(muscle_da$gene_symbol))

#input is all chip-bound targets of MEF2A
chip_targets <- mef2a$GeneSym[mef2a$any==1]

#filter for entries with gene symbols that are ChIP-bound
ee_targets <- ee$name[ee$name%in%background & ee$name%in%chip_targets]
re_targets <- re$name[re$name%in%background & re$name%in%chip_targets]
either <- either[either%in%background & either%in%chip_targets]
both <- both[both%in%background & both%in%chip_targets]

#perform ORA on each list
targets <- list("EE"=ee_targets,"RE"=re_targets, "Either"=either,"Both"=both)
ora_results <- lapply(targets,function(x){
  return(run_ORA(x,background))
})

#for paper make a joint heatmap with interesting terms from C2 and C5
re_C5 <- ora_results$RE[ora_results$RE$collection=="C5",]
ee_C5 <- ora_results$EE[ora_results$EE$collection=="C5",]
either_C5 <- ora_results$Either[ora_results$Either$collection=="C5",]
re_C2 <- ora_results$RE[ora_results$RE$collection=="C2",]
ee_C2 <- ora_results$EE[ora_results$EE$collection=="C2",]
either_C2 <- ora_results$Either[ora_results$Either$collection=="C2",]
ora.results.filt <- list("RE"=rbind(re_C5,re_C2),
                         "EE"=rbind(ee_C5,ee_C2),
                         "Either"=rbind(either_C5,either_C2))

#use just either as that encompasses everything we want, and take only the pathways with at least 10 genes to filter evne further
#for C5, lots of terms, focus on GOBP
#additionally, remove negative/positive regulation terms as they are more broad and skew the p-values
sig_terms_either_C5 <- either_C5$set_short[either_C5$adj_p_value<0.05 & either_C5$set_size_in_input >= 10 & grepl("GOBP",either_C5$set_short) & !grepl("POSITIVE_REGULATION_|NEGATIVE_REGULATION_",either_C5$set_short)]
sig_terms_either_C2 <- either_C2$set_short[either_C2$adj_p_value<0.05 & either_C2$set_size_in_input >= 10]
sig_terms <- c(sig_terms_either_C5,sig_terms_either_C2)

#filter results to just these terms and calculate the input ratio for the heatmap
ora.results.heatmap <- lapply(names(ora.results.filt),function(x){
  table <- ora.results.filt[[x]] %>% dplyr::filter(set_short %in% sig_terms) %>% dplyr::select("set_short","set_size_in_input","input_size","adj_p_value")
  table$input_ratio <- table$set_size_in_input/table$input_size
  colnames(table)[colnames(table)!="set_short"] <- paste(colnames(table)[colnames(table)!="set_short"],x,sep=".")
  return(table)
})
names(ora.results.heatmap) <- names(ora.results.filt)
ee <- ora.results.heatmap[[1]]
re <- ora.results.heatmap[[2]]
either <- ora.results.heatmap[[3]]
data_heatmap <- full_join(re,ee,by="set_short")
data_heatmap <- full_join(data_heatmap,either,by="set_short")
rownames(data_heatmap) <- data_heatmap$set_short

#separate out adjusted p-value and input_ratio for heatmap
#ratio_heatmap <- data_heatmap %>% dplyr::select(contains("input_ratio"))
pval_heatmap <- data_heatmap %>% dplyr::select(contains("adj_p_value"))

#convert to -log10(pval)
pval_heatmap <- -log10(pval_heatmap)

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
                            r = norm_per*0.03,
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
pdf("MEF2A_ChIP_ORA_heatmap.pdf",height=16,width=10)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "left",
     annotation_legend_list = lgd_list,
     legend_gap = unit(10, "mm"))
dev.off()

# #lots of cool stuff in C2 and C5, EE only, RE only, and either
# #let's make a heatmap for both C2 and C5, starting with the significant terms in either
# 
# #C2 heatmap
# re_C2 <- ora_results$RE$C2
# ee_C2 <- ora_results$EE$C2
# either_C2 <- ora_results$Either$C2
# sig_terms <- either_C2$set_short[either_C2$adj_p_value<0.05]
# 
# #get data for heatmap
# re_heatmap <- re_C2[re_C2$set_short%in%sig_terms,c("set_short","set_size_in_input","input_size","adj_p_value")]
# re_heatmap$input_ratio <- re_heatmap$set_size_in_input/re_heatmap$input_size
# 
# ee_heatmap <- ee_C2[ee_C2$set_short%in%sig_terms,c("set_short","set_size_in_input","input_size","adj_p_value")]
# ee_heatmap$input_ratio <- ee_heatmap$set_size_in_input/ee_heatmap$input_size
# 
# either_heatmap <- either_C2[either_C2$set_short%in%sig_terms,c("set_short","set_size_in_input","input_size","adj_p_value")]
# either_heatmap$input_ratio <- either_heatmap$set_size_in_input/either_heatmap$input_size
# 
# data_heatmap <- list(ee_heatmap,re_heatmap,either_heatmap) %>% purrr::reduce(full_join, by="set_short")
# rownames(data_heatmap) <- data_heatmap$set_short
# 
# #filter for sets with at least 5 genes
# data_heatmap <- data_heatmap[data_heatmap$set_size_in_input>=5,]
# 
# #separate out adjusted p-value and input_ratio for heatmap
# ratio_heatmap <- data_heatmap %>% dplyr::select("input_ratio.x","input_ratio.y","input_ratio")
# pval_heatmap <- data_heatmap %>% dplyr::select("adj_p_value.x","adj_p_value.y","adj_p_value")
# 
# #Heatmap dimensions
# a <- unit(5, "mm")
# width <- ncol(pval_heatmap) * a
# height <- nrow(pval_heatmap) * a
# 
# #color scale
# col_fun <- circlize::colorRamp2(
#   breaks = c(0,0.05,1),
#   colors = c("#b2182b", "#d6604d", "white"))
# 
# #plot heatmap
# ht <- Heatmap(as.matrix(pval_heatmap),
#               heatmap_legend_param = list(
#                 title = "adjusted p-value",
#                 legend_direction = "vertical",
#                 legend_width = unit(50, "mm"),
#                 at = c(0,0.05,1)
#               ),
#               width = width,
#               height = height,
#               #top_annotation = ha,
#               #left_annotation = row_ha,
#               cluster_columns = F,
#               cluster_rows=T,
#               show_column_names = T,
#               show_row_names = T,
#               column_names_side="top",
#               col = col_fun,
#               rect_gp = gpar(type = "none"),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 
#                 #Heatmap grid
#                 grid.rect(x = x, y = y, width, height,
#                           gp = gpar(col = "#555555", fill = "white"))
#                 
#                 #Point size
#                 grid.circle(x = x,
#                             y = y,
#                             r = ratio_heatmap[i,j] * a * 3,
#                             gp = gpar(fill = col_fun(pval_heatmap[i, j]),
#                                       col = "#555555"))
#                 
#               },
#               column_title_rot = 90,
#               row_title_rot = 0,
#               column_names_gp = grid::gpar(fontsize = 8),
#               row_names_gp = grid::gpar(fontsize = 8)
# )
# 
# lgd_list = list(
#   Legend(labels = c("< 5","10","15"),
#          title = "% associated genes",
#          type = "points",
#          # ncol = 5,
#          pch = 16,
#          size = unit(c(1,3,5), 'mm'),
#          legend_gp = gpar(
#            col = rep("black",3)
#          ),
#          background = 'white'
#   )
# )
# 
# #save to PDF
# pdf("MEF2A_ORA_C2_heatmap.pdf",height=8,width=12)
# draw(ht,
#      merge_legends = TRUE,
#      heatmap_legend_side = "left",
#      annotation_legend_list = lgd_list,
#      legend_gap = unit(10, "mm"))
# dev.off()
# 
# #C5 heatmap
# re_C5 <- ora_results$RE$C5
# ee_C5 <- ora_results$EE$C5
# either_C5 <- ora_results$Either$C5
# sig_terms <- either_C5$set_short[either_C5$adj_p_value<0.05]
# 
# #get data for heatmap
# re_heatmap <- re_C5[re_C5$set_short%in%sig_terms,c("set_short","set_size_in_input","input_size","adj_p_value")]
# re_heatmap$input_ratio <- re_heatmap$set_size_in_input/re_heatmap$input_size
# 
# ee_heatmap <- ee_C5[ee_C5$set_short%in%sig_terms,c("set_short","set_size_in_input","input_size","adj_p_value")]
# ee_heatmap$input_ratio <- ee_heatmap$set_size_in_input/ee_heatmap$input_size
# 
# either_heatmap <- either_C5[either_C5$set_short%in%sig_terms,c("set_short","set_size_in_input","input_size","adj_p_value")]
# either_heatmap$input_ratio <- either_heatmap$set_size_in_input/either_heatmap$input_size
# 
# data_heatmap <- list(ee_heatmap,re_heatmap,either_heatmap) %>% purrr::reduce(full_join, by="set_short")
# rownames(data_heatmap) <- data_heatmap$set_short
# 
# #filter for sets with at least 5 genes
# data_heatmap <- data_heatmap[data_heatmap$set_size_in_input>=5,]
# 
# #separate out adjusted p-value and input_ratio for heatmap
# ratio_heatmap <- data_heatmap %>% dplyr::select("input_ratio.x","input_ratio.y","input_ratio")
# pval_heatmap <- data_heatmap %>% dplyr::select("adj_p_value.x","adj_p_value.y","adj_p_value")
# 
# #Heatmap dimensions
# a <- unit(3, "mm")
# width <- ncol(pval_heatmap) * a
# height <- nrow(pval_heatmap) * a
# 
# #color scale
# col_fun <- circlize::colorRamp2(
#   breaks = c(0,0.05,1),
#   colors = c("#b2182b", "#d6604d", "white"))
# 
# #plot heatmap
# ht <- Heatmap(as.matrix(pval_heatmap),
#               heatmap_legend_param = list(
#                 title = "adjusted p-value",
#                 legend_direction = "vertical",
#                 legend_width = unit(50, "mm"),
#                 at = c(0,0.05,1)
#               ),
#               width = width,
#               height = height,
#               #top_annotation = ha,
#               #left_annotation = row_ha,
#               cluster_columns = F,
#               cluster_rows=T,
#               show_column_names = T,
#               show_row_names = T,
#               column_names_side="top",
#               col = col_fun,
#               rect_gp = gpar(type = "none"),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                 
#                 #Heatmap grid
#                 grid.rect(x = x, y = y, width, height,
#                           gp = gpar(col = "#555555", fill = "white"))
#                 
#                 #Point size
#                 grid.circle(x = x,
#                             y = y,
#                             r = ratio_heatmap[i,j] * a,
#                             gp = gpar(fill = col_fun(pval_heatmap[i, j]),
#                                       col = "#555555"))
#                 
#               },
#               column_title_rot = 90,
#               row_title_rot = 0,
#               column_names_gp = grid::gpar(fontsize = 8),
#               row_names_gp = grid::gpar(fontsize = 8)
# )
# 
# lgd_list = list(
#   Legend(labels = c("< 5","10","20","30"),
#          title = "% associated genes",
#          type = "points",
#          # ncol = 5,
#          pch = 16,
#          size = unit(1:4, 'mm'),
#          legend_gp = gpar(
#            col = rep("black",3)
#          ),
#          background = 'white'
#   )
# )
# 
# #save to PDF
# pdf("MEF2A_ORA_C5_heatmap.pdf",height=12,width=12)
# draw(ht,
#      merge_legends = TRUE,
#      heatmap_legend_side = "left",
#      annotation_legend_list = lgd_list,
#      legend_gap = unit(10, "mm"))
# dev.off()

#heatmap of the ChIP-bound genes
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
  breaks = c(-2,0,2),
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
pdf("MEF2A_ChIP_targets_heatmap.pdf",height=30,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#for a list of terms, produce a leading edge heatmap of the genes in those terms
sig_terms <- c("WP_VEGFA_VEGFR2_SIGNALING",
               "REACTOME_DISEASES_OF_SIGNAL_TRANSDUCTION_BY_GROWTH_FACTOR_RECEPTORS_AND_SECOND_MESSENGERS",
               #"GOBP_NEGATIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II",
               #"GOBP_NEGATIVE_REGULATION_OF_RNA_BIOSYNTHETIC_PROCESS",
               #"GOBP_POSITIVE_REGULATION_OF_RNA_METABOLIC_PROCESS",
               #"GOBP_POSITIVE_REGULATION_OF_TRANSCRIPTION_BY_RNA_POLYMERASE_II",
               #"PID_IL6_7_PATHWAY",
               #"PID_HIF1_TFPATHWAY",
               #"WP_WNT_SIGNALING",
               #"WP_TGF_BETA_SIGNALING_PATHWAY",
               "GOBP_REGULATION_OF_AUTOPHAGY")

for(i in sig_terms){
  #get genes in this term
  db <- str_split_i(i,"_",1)
  sig_list <- MOLECULAR_SIGNATURES[[db]][[i]]
  chip_sig_wide <- any_chip_wide[any_chip_wide$gene_symbol%in%sig_list,]
  
  #separate out FC and adj p-value for the heatmap
  fc <- chip_sig_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
  pval <- chip_sig_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
  rownames(fc) = rownames(pval) = chip_sig_wide$gene_symbol
  
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
  pdf(glue("{i}.pdf"),height=8,width=8)
  draw(ht,
       merge_legends = TRUE,
       heatmap_legend_side = "right",
       legend_gap = unit(10, "mm"))
  dev.off()
  
}
