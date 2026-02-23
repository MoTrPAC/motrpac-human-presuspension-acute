#comprehensive R script that generates most panels for Figure 7/S7
#Greg will have his own script(s) which generate the rest of the panels

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

#Figure 7B, FOXO1/FOXO3 phosphosite heatmap

#load DA and filter for FOXO phosphosites
da <- MUSCLE_PROT_PH_DA
da <- da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
da <- left_join(da,HUMAN_FEATURE_TO_GENE,by="feature_id")

FOXO <- da %>% dplyr::filter(gene_symbol%in%c("FOXO1","FOXO3"))
FOXO_filter <- FOXO %>% dplyr::select("feature_id","gene_symbol","contrast_short","logFC","p_value","adj_p_value")
FOXO_wide <- pivot_wider(FOXO_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
FOXO_wide <- FOXO_wide[,grepl("feature_id|gene_symbol|delta-delta",colnames(FOXO_wide))]
sites <- unlist(strsplit(as.character(FOXO_wide$feature_id),"_"))
rownames(FOXO_wide) <- paste(FOXO_wide$gene_symbol,sites[seq(2,length(sites),by=2)],sep=".")

#separate out FC and adj p-value for the heatmap
fc <- FOXO_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- FOXO_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = rownames(FOXO_wide)

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
  breaks = c(-0.5,0,0.5),
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
              row_split=ifelse(grepl("FOXO1",rownames(fc)),"FOXO1","FOXO3"),
              cluster_columns = F,
              cluster_rows=T,
              show_column_names = F,
              show_row_names = T,
              col = col_fun,
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8),
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #grid border
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "#555555"))
                
                #Point size
                if(pval[i, j] < 0.05) {
                  gb = textGrob("*")
                  gb_w = convertWidth(grobWidth(gb), "mm")
                  gb_h = convertHeight(grobHeight(gb), "mm")
                  grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                }
              }
)

#save to PDF
pdf("Fig7B_FOXO_phospho_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#Figure 7C, ZEB1 phosphosite heatmap
#load DA and filter for ZEB1 phosphosites
da <- MUSCLE_PROT_PH_DA
da <- da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
da <- left_join(da,HUMAN_FEATURE_TO_GENE,by="feature_id")

ZEB1 <- da %>% dplyr::filter(gene_symbol=="ZEB1")
ZEB1_filter <- ZEB1 %>% dplyr::select("feature_id","gene_symbol","contrast_short","logFC","p_value","adj_p_value")
ZEB1_wide <- pivot_wider(ZEB1_filter, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
ZEB1_wide <- ZEB1_wide[,grepl("feature_id|gene_symbol|delta-delta",colnames(ZEB1_wide))]
sites <- unlist(strsplit(as.character(ZEB1_wide$feature_id),"_"))
rownames(ZEB1_wide) <- paste(ZEB1_wide$gene_symbol,sites[seq(2,length(sites),by=2)],sep=".")

#separate out FC and adj p-value for the heatmap
fc <- ZEB1_wide %>% dplyr::select(contains("logFC")) %>% as.matrix()
pval <- ZEB1_wide %>% dplyr::select(contains("adj_p_value")) %>% as.matrix()
rownames(fc) = rownames(pval) = rownames(ZEB1_wide)

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
  breaks = c(-0.5,0,0.5),
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
              column_title_rot = 90,
              row_title_rot = 0,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8),
              cell_fun = function(j, i, x, y, width, height, fill) {
                
                #grid border
                grid.rect(x = x, y = y, width, height,
                          gp = gpar(col = "#555555"))
                
                #Point size
                if(pval[i, j] < 0.05) {
                  gb = textGrob("*")
                  gb_w = convertWidth(grobWidth(gb), "mm")
                  gb_h = convertHeight(grobHeight(gb), "mm")
                  grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
                }
              }
)

#save to PDF
pdf("Fig7C_ZEB1_phospho_heatmap.pdf",height=8,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()

#Figure S7A: ORA on ZEB1 targets

#load ZEB1 targets
ZEB1 <- read_csv("ZEB1-targets.csv") %>% dplyr::select("Target") %>% unique() %>% dplyr::rename(GeneSym=Target)

#ENCODE resource
encode <- read_csv("external_data/ZEB1-encode.csv") %>% dplyr::rename("ENCODE_2015"="ZEB1")
ZEB1 <- left_join(ZEB1,encode,by="GeneSym")

#Rosmaninho et al 2018
ros_chip <- read.xlsx("external_data/embj201797115-sup-0004-tableev2.xlsx")
#for each gene, check if peak within 1kb of TSS. if yes, mark 1.
ZEB1$Rosmainho_et_al_2018 <- unlist(lapply(ZEB1$GeneSym,function(x){
  peaks <- ros_chip[ros_chip$Gene.Symbol==x,]
  #if not in dataset, skip
  if(dim(peaks)[1]==0){
    return(0)
  }
  # Extract numbers (handles both + and - signs) and check if abs value <= 1000
  numbers <- as.numeric(gsub(".*\\(([+-]\\d+)\\).*", "\\1", peaks$`Associated.Peak(s)`))
  ifelse(any(abs(numbers) <= 1000),1,0)
})) #this really doesn't add any more targets, so we can scrap it in the final version
ZEB1$any <- ifelse(rowSums(ZEB1[,-1],na.rm=T)>0,1,0)

#save!
write_csv(ZEB1,"ZEB1_targets_with_peaks.csv")

#perform ORA
#load list of TFEB targets in each network (EE or RE)
ee <- read_csv("ZEB1-targets-EE.csv")
re <- read_csv("ZEB1-targets-RE.csv")
either <- c(ee$Target,re$Target) #names in either EE OR RE
both <- either[duplicated(either)] #names in both EE and RE

#get background list which is all detected transcripts that went into the DA
muscle_da <- MUSCLE_TRNSCRPT_DA
muscle_da <- muscle_da %>% dplyr::filter(contrast_category=="EE-CON" | contrast_category=="RE-CON")
muscle_da <- right_join(HUMAN_FEATURE_TO_GENE,muscle_da,by="feature_id")
background <- as.character(unique(muscle_da$gene_symbol))

#input is all chip-bound targets of ZEB1
chip_targets <- ZEB1$GeneSym[ZEB1$any==1]

#filter for entries with gene symbols that are ChIP-bound
ee_targets <- ee$Target[ee$Target%in%background & ee$Target%in%chip_targets]
re_targets <- re$Target[re$Target%in%background & re$Target%in%chip_targets]
either <- either[either%in%background & either%in%chip_targets]
#both <- both[both%in%background & both%in%chip_targets]

#perform ORA on each list
#don't do both as there are too few (3)
targets <- list("EE"=ee_targets,"RE"=re_targets, "Either"=either)
ora_results <- lapply(targets,function(x){
  return(run_ORA(x,background))
})
write.xlsx(ora_results,"ZEB1_ORA_results.xlsx")

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

#use just either as that encompasses everything we want, and take only the pathways with at least 5 genes to filter evne further
#additionally, remove negative/positive regulation terms as they are more broad and skew the p-values
sig_terms_either_C5 <- either_C5$set_short[either_C5$adj_p_value<0.05 & either_C5$set_size_in_input >= 5 & !grepl("POSITIVE_REGULATION_|NEGATIVE_REGULATION_",either_C5$set_short)]
sig_terms_either_C2 <- either_C2$set_short[either_C2$adj_p_value<0.05 & either_C2$set_size_in_input >= 5]
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

# Create the horizontal bar chart
either$adj_p_value.Either= -log10(either$adj_p_value.Either)
pdf(glue("FigS7A_ZEB1-ChIP-ORA-bar-chart.pdf"),height=6,width=10)
ggplot(either, aes(x = adj_p_value.Either, y = reorder(set_short, adj_p_value.Either))) +
  geom_col(aes(fill = adj_p_value.Either),color = "black", linewidth=0.5) +
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

#for a list of terms, produce a leading edge heatmap of the genes in those terms
any_chip <- data.frame("gene_symbol"=chip_targets)
any_chip <- left_join(any_chip,muscle_da,by="gene_symbol") %>% dplyr::select("gene_symbol","contrast_short","logFC","p_value","adj_p_value")
any_chip_wide <- pivot_wider(any_chip, names_from=contrast_short, values_from=c("logFC","p_value","adj_p_value")) 
any_chip_wide <- any_chip_wide[,grepl("gene_symbol|delta-delta",colnames(any_chip_wide))]
rownames(any_chip_wide) <- any_chip_wide$gene_symbol

sig_terms <- c("GOBP_CELLULAR_RESPONSE_TO_STRESS"
)

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
  pdf(glue("FigS7A_{i}.pdf"),height=8,width=8)
  draw(ht,
       merge_legends = TRUE,
       heatmap_legend_side = "right",
       legend_gap = unit(10, "mm"))
  dev.off()
  
}

#Figure S7C, heatmap of FOXO1 targets

#load FOXO1 targets
FOXO1 <- read_csv("FOXO1-targets.csv") %>% dplyr::select("Target") %>% unique() %>% dplyr::rename(GeneSym=Target)

#santini et al 2024
#mouse ChIP, must map to human
mouse_chip <- read.xlsx("external_data/santini_et_al_2024.xlsx")
#need to convert mouse to human gene symbols
mapIt <- function(x) {
  require("Orthology.eg.db", character.only = TRUE)
  require("org.Mm.eg.db", character.only = TRUE)
  require("org.Hs.eg.db", character.only = TRUE)
  mouse <- mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL")
  mapped <- AnnotationDbi::select(Orthology.eg.db, mouse,"Homo.sapiens","Mus.musculus")
  human <- mapIds(org.Hs.eg.db, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  human <- do.call(c, lapply(human, function(x) if(is.null(x)) return(NA) else return(x)))
  cbind(x, mapped, human)
}
symbols <- mapIt(mouse_chip$geneId) %>% dplyr::rename(GeneSym=human,geneId=x)
mouse_chip <- left_join(mouse_chip,symbols,by="geneId",relationship="many-to-many")
mouse_chip <- mouse_chip[!is.na(mouse_chip$GeneSym),]

#for each gene, check if peak. if yes, mark 1.
FOXO1$Santini_et_al_2024 <- unlist(lapply(FOXO1$GeneSym,function(x){
  peaks <- mouse_chip[mouse_chip$GeneSym==x,]
  if(dim(peaks)[1]==0){
    return(0)
  }
  ifelse(peaks$is_FoxO1_target,1,0)
}))

#peaks from Greg already filtered to within 1kb of the TSS
atac_peaks <- read_csv("foxo1bindingpeaks_proxprom.csv") %>% dplyr::rename(GeneSym=gene_symbol)

FOXO1$ATAC <- unlist(lapply(FOXO1$GeneSym,function(x){
  peaks <- atac_peaks[atac_peaks$GeneSym==x,] %>% na.omit()
  ifelse(dim(peaks)[1]==0,0,1)
}))

#integrate
FOXO1$any <- ifelse(rowSums(FOXO1[,-1],na.rm=T)>0,1,0)

#save!
write_csv(FOXO1,"FOXO1_targets_with_peaks.csv")

#create heatmap of ChIP targets
chip_targets <- FOXO1$GeneSym[FOXO1$any==1]

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

#filter for FDR < 0.01
pval <- pval[rowSums(pval<=0.01)>0,]
fc <- fc[rownames(fc)%in%rownames(pval),]

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
  breaks = c(min(fc),0,2),
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
pdf("FigS7C_FOXO1_ChIP_targets_heatmap_FDR_0.01.pdf",height=15,width=8)
draw(ht,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     legend_gap = unit(10, "mm"))
dev.off()