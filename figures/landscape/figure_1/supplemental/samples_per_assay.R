#The goal of this document is to identify the numbers of samples present in every
#assay at every timepoint. This corresponds to supplementary table S1

#the sum total number of samples relevant to Pre-COVID is equal to the number of
#of samples that pass QC + the number of samples that don't. So I load qc and subset
# using the OUTLIERS object, which represent the samples that don't pass QC.

library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)

# Note:

# This figure(s) requires access to `MotrpacHumanPreSuspensionData` to fully recreate.

# config = jsonlite::fromJSON("~/config.json")
# repo_local_dir = file.path(config$precovid_repo_path, "data", "tmp")

all_dataset = load_qc(epigen = TRUE,
                      repo_local_dir = "~/Downloads/",
                      remove_redundant_metab = TRUE)

split_tissues = unlist(all_dataset, recursive = FALSE)
split_tissue_assay = unlist(split_tissues, recursive = FALSE)
qc_norm_only = split_tissue_assay[grep("qc_norm", names(split_tissue_assay))]

#makes a big list of tissues, assays
samples_per_plat = data.frame(
  list_name = names(qc_norm_only),
  vialLabels = sapply(qc_norm_only, function(x) paste(colnames(x), collapse = ";"))
) %>%
  mutate(
    tissue = sapply(strsplit(list_name, "\\."), `[`, 1),
    assay = sapply(strsplit(list_name, "\\."), `[`, 2)
  ) %>%
  select(-list_name)

#then uses the pheno file to annotate group and timepoint info and split the vector
qc_norm_counts = lapply(seq_len(nrow(samples_per_plat)), function(row) {
  vials_vec = strsplit(samples_per_plat[row, "vialLabels"], ";")[[1]] %>%
    base::trimws()

  pheno$data %>%
    filter(vialLabel %in% vials_vec) %>%
    group_by(randomGroupCode, Timepoint) %>%
    summarise(
      n = n(),
      vialLabels = paste(vialLabel, collapse = ","),
      .groups = "drop"
    ) %>%
    mutate(
      tissue = samples_per_plat$tissue[row],
      assay = samples_per_plat$assay[row]
    )
}) %>%
  bind_rows() %>%
  select(tissue, assay, randomGroupCode, Timepoint, n, vialLabels)

full_counts_by_tp = qc_norm_counts %>%
  select(-vialLabels) %>%
  pivot_wider(names_from = c("tissue", "Timepoint"),
              values_from = "n") %>%
  arrange(assay)

saveRDS(full_counts_by_tp, file.path(repo_local_dir, "figures", "table_S1_participants_per_tp.RDS"))
write.csv(full_counts_by_tp, file.path(repo_local_dir, "figures", "table_S1_participants_per_tp.csv"),
          row.names = FALSE)


#then we append outlier info to append some extra info
outliers_append = OUTLIERS %>%
  select(-reason) %>%
  group_by(tissue, ome) %>%
  count() %>%
  pivot_wider(names_from = "tissue",
              values_from = "n") %>%
  arrange(ome)


saveRDS(outliers_append, file.path(repo_local_dir, "figures", "table_S1_outliers.RDS"))
write.csv(outliers_append, file.path(repo_local_dir, "figures", "table_S1_outliers.csv"),
          row.names = FALSE)
