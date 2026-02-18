#nfic target analysis

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

#load NFIC targets
nfic <- read_csv("NFIC_targets.csv") %>% dplyr::select("name") %>% dplyr::filter(name!="NFIC") %>% dplyr::rename(GeneSym=name)

#ENCODE resource
encode <- read_csv("NFIC-encode.csv") %>% dplyr::rename("ENCODE_2015"="NFIC")
nfic <- left_join(nfic,encode,by="GeneSym")

#Pinar et al 2025
pinar <- read_csv("NFIC_pinar_et_al_2025.csv") 
nfic$Pinar_2025 <- ifelse(nfic$GeneSym%in%pinar$NFIC.targets,1,NA)
nfic$any <- ifelse(rowSums(nfic[,-1],na.rm=T)>0,1,0)

#save!
write_csv(nfic,"NFIC_targets_with_peaks.csv")

#perform ORA
#load list of NFIC targets in each network (EE or RE)
ee <- read_csv("NFIC_targets_EE.csv")
re <- read_csv("NFIC_targets_RE.csv")
either <- c(ee$Target,re$Target) #names in either EE OR RE
both <- either[duplicated(either)] #names in both EE and RE

#get background list which is all detected transcripts that went into the DA
muscle_da <- MUSCLE_TRNSCRPT_DA
muscle_da <- muscle_da %>% dplyr::filter(contrast_type=="exercise_with_controls")
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")
muscle_da$contrast_short = as.character(muscle_da$contrast_short)
background <- as.character(unique(muscle_da$gene_symbol))

#input is all chip-bound targets of nfic
chip_targets <- nfic$GeneSym[nfic$any==1]

#filter for entries with gene symbols that are ChIP-bound
ee_targets <- ee$Target[ee$Target%in%background & ee$Target%in%chip_targets]
re_targets <- re$Target[re$Target%in%background & re$Target%in%chip_targets]
either <- either[either%in%background & either%in%chip_targets]
both <- both[both%in%background & both%in%chip_targets]

#perform ORA on each list
targets <- list("EE"=ee_targets,"RE"=re_targets, "Either"=either,"Both"=both)
ora_results <- lapply(targets,function(x){
  return(run_ORA(x,background))
})
write.xlsx(ora_results,"NFIC-target-ORA-results.xlsx")

#for paper make a joint heatmap with interesting terms from C2, C5, and MITOCARTA
re_C5 <- ora_results$RE[ora_results$RE$collection=="C5",]
ee_C5 <- ora_results$EE[ora_results$EE$collection=="C5",]
either_C5 <- ora_results$Either[ora_results$Either$collection=="C5",]
re_C2 <- ora_results$RE[ora_results$RE$collection=="C2",]
ee_C2 <- ora_results$EE[ora_results$EE$collection=="C2",]
either_C2 <- ora_results$Either[ora_results$Either$collection=="C2",]
re_mito <- ora_results$RE[ora_results$RE$collection=="MITOCARTA",]
ee_mito <- ora_results$EE[ora_results$EE$collection=="MITOCARTA",]
both_mito <- ora_results$Both[ora_results$Both$collection=="MITOCARTA",]
ora.results.filt <- list("RE"=rbind(re_C5,re_C2,re_mito),
                         "EE"=rbind(ee_C5,ee_C2,ee_mito),
                         "Either"=rbind(either_C5,either_C2,either_mito))

#use just either as that encompasses everything we want, and take only the pathways with at least 30 genes to filter evne further
#for C5, lots of terms, focus on GOBP
#use more stringent cutoff of 0.01
sig_terms_either_C5 <- either_C5$set_short[either_C5$adj_p_value<0.01 & either_C5$set_size_in_input >= 30 & grepl("GOBP",either_C5$set_short)]
sig_terms_either_C2 <- either_C2$set_short[either_C2$adj_p_value<0.01 & either_C2$set_size_in_input >= 30 & grepl("REACTOME",either_C2$set_short)]
sig_terms_either_mito <- either_mito$set_short[either_mito$adj_p_value<0.01 & either_mito$set_size_in_input >= 30]
sig_terms <- c(sig_terms_either_C5,sig_terms_either_C2,sig_terms_either_mito)

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
pdf("NFIC_ChIP_ORA_heatmap.pdf",height=16,width=10)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "left",
     annotation_legend_list = lgd_list,
     legend_gap = unit(10, "mm"))
dev.off()

#ChIP-bound genes
any_chip <- data.frame("gene_symbol"=chip_targets)
any_chip <- left_join(any_chip,muscle_da,by="gene_symbol") %>% dplyr::select("gene_symbol","ensembl_gene","contrast_short","logFC","p_value","adj_p_value")
#there are two duplicated gene symbols, MATR3 and MKKS. for these, use the first reported gene ID as this is the canonical variant
dupes <- any_chip[duplicated(any_chip[,c("gene_symbol","contrast_short")]),]
any_chip <- any_chip %>% dplyr::filter(!ensembl_gene %in% dupes$ensembl_gene) %>% dplyr::mutate(ensembl_gene=NULL)
any_chip_wide <- pivot_wider(any_chip, id_cols=gene_symbol, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
any_chip_wide <- any_chip_wide[,grepl("gene_symbol|delta-delta",colnames(any_chip_wide))] %>% data.frame()
rownames(any_chip_wide) <- any_chip_wide$gene_symbol

#for a list of terms, produce a leading edge heatmap of the genes in those terms
sig_terms <- c("GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
               "REACTOME_MRNA_SPLICING",
               "Mitochondrial ribosome",
               "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
               "GOBP_AEROBIC_RESPIRATION",
               "Metals and cofactors",
               "GOBP_CELLULAR_RESPIRATION",
               "OXPHOS",
               "GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MEMBRANE",
               "GOBP_MITOCHONDRIAL_GENE_EXPRESSION",
               "REACTOME_TRANSLATION")

for(i in sig_terms){
  #get genes in this term
  if(grepl("_",i)){
    db <- str_split_i(i,"_",1)
    sig_list <- MOLECULAR_SIGNATURES[[db]][[i]]
  }else{
    sig_list <- MOLECULAR_SIGNATURES[["MITOCARTA"]][[paste("MITOCARTA",i,sep="_")]]
  }
  chip_sig_wide <- any_chip_wide[any_chip_wide$gene_symbol%in%sig_list,]
  
  #separate out FC and adj p-value for the heatmap
  fc <- chip_sig_wide %>% dplyr::select(contains("logFC")) %>% data.frame()
  pval <- chip_sig_wide %>% dplyr::select(contains("adj_p_value")) %>% data.frame()
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

#remake the heatmap for the metals and cofactors with just a subset of the genes for the main figure
#select only the genes that are up at 24hr in EE
sig_list = c("CYCS","ALAS1","SLC25A20","COQ3","COX10","NTHL1","SIRT5","HSPA9","FDX1","CISD1","NUDT8","BOLA3","ETFDH","GRPEL1","SLC25A3","NDUFAB1","CISD3","SDHB","CYC1","FAM210B","NDUFS7","CYB5B","NDUFV1","NDUFS1","COQ9","NFS1","GLRX5","SLC25A42")
chip_sig_wide <- any_chip_wide[any_chip_wide$gene_symbol%in%sig_list,]

#separate out FC and adj p-value for the heatmap
fc <- chip_sig_wide %>% dplyr::select(contains("logFC")) %>% data.frame()
pval <- chip_sig_wide %>% dplyr::select(contains("adj_p_value")) %>% data.frame()
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
pdf(glue("Metals and cofactors filtered.pdf"),height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()
