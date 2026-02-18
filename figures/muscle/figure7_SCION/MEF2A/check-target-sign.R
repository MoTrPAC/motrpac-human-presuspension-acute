#check if target FC sign matches MEF2A FC sign
library(tidyverse)
library(MotrpacHumanPreSuspension)

#get edge table which has site information
#filter for edges only from SCION that come from MEF2A
table <- read_csv("MEF2A-edge-table.csv") %>% dplyr::filter(interaction=="regulates" & grepl("MEF2A ", name)) 

#load the DA results for phospho and transcript
phos <- MUSCLE_PROT_PH_DA %>% dplyr::mutate(contrast_short=as.character(contrast_short), contrast_type= as.character(contrast_type)) %>% dplyr::filter(contrast_type=="exercise_with_controls" & grepl("delta\\-delta",contrast_short))
trans <- MUSCLE_TRNSCRPT_DA %>% dplyr::mutate(contrast_short=as.character(contrast_short), contrast_type= as.character(contrast_type)) %>% dplyr::filter(contrast_type=="exercise_with_controls" & grepl("delta\\-delta",contrast_short))

#get the info for MEF2A
phos_mef2a <- phos %>% dplyr::filter(feature_id %in% table$regulator_id & adj_p_value < 0.05)

#get the info for all the targets
trans_targets <- trans %>% dplyr::filter(feature_id %in% table$target_id & adj_p_value < 0.05)

#for each target, get the sign of FC for that modality and check if it matches the site
for(i in 1:dim(table)[1]){
  site <- table$regulator_id[i]
  target <- table$target_id[i]
  site_sign <- sign(phos_mef2a$logFC)[phos_mef2a$feature_id==site]
  target_sign <- sign(trans_targets$logFC[trans_targets$feature_id==target])
  table$MEF2A_sign[i] <- site_sign
  table$target_sign[i] <- paste(target_sign,collapse="|")
  #if multiple target FCs, take the consensus. if it is split 50/50, it's unknown.
  counts <- table(target_sign)
  max_count <- max(counts)
  majority_values <- names(counts[counts == max_count])
  if (length(majority_values) == 1 && max_count > length(target_sign) / 2) {
    target_sign_consensus <- as.numeric(majority_values)
  } else {
    target_sign_consensus <- NA  # 
  }
  table$target_sign_consensus[i] <- target_sign_consensus
  table$directionality[i] <- ifelse(is.na(target_sign_consensus),NA,
                                    ifelse(target_sign_consensus==site_sign,1,-1))
}
write_csv(table,"MEF2A-edge-table-with-directionality.csv")
