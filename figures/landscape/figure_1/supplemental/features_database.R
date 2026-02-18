#the goal of this document is to track every set of features that is used
#in (almost) any aspect of the analysis.
#this corresponds to table S2 in the supplementary files

#the sets I'm choosing to include:
#Named metabolites in any aspect of the qc_norm. If a metabolite is measured on
#multiple platforms, only the platform with the lowest CV is chosen, because
#that is how it is most commonly presented in our dataset.

library(MotrpacHumanPreSuspensionAnalysis)
library(dplyr)
library(tidyr)
config = jsonlite::fromJSON("~/config.json")
repo_local_dir = file.path(config$precovid_repo_path, "data", "tmp")

all_dataset = load_qc(epigen = TRUE,
                      repo_local_dir = repo_local_dir,
                      remove_redundant_metab = TRUE)

#I previously was using combine_qc_matrixes. this doesn't seem to work, i believe
#because of missingness and subsetting to shared participants.
all_tissue_data = data.frame()
for(specific_tissue in tissue_available_list()){
  for(assay in names(all_dataset[[specific_tissue]])){
    qc_norm = all_dataset[[specific_tissue]][[assay]][["qc_norm"]]
    if(nrow(qc_norm) == 0) next
    tissue_data = data.frame(feature_id = rownames(qc_norm),
                             assay = assay,
                             tissue = specific_tissue)
    all_tissue_data = rbind(all_tissue_data, tissue_data)
  }
}

saveRDS(all_tissue_data, file.path(repo_local_dir, "figures", "table_S2.RDS"))
# all_tissue_data = readRDS(file.path(repo_local_dir, "figures", "table_S2.RDS"))


#if you want the total number of features measured on each tissue/ome
# all_tissue_data = readRDS(file.path(repo_local_dir, "figures", "table_S2.RDS"))
counts_per_tissue_assay = all_tissue_data %>%
  group_by(tissue, assay) %>%
  count() %>%
  tidyr::pivot_wider(names_from = "tissue",
                     values_from = "n") %>%
  arrange(assay)
write.csv(counts_per_tissue_assay,
          file = file.path(repo_local_dir, "figures", "table_S2.csv"),
          row.names = FALSE)
