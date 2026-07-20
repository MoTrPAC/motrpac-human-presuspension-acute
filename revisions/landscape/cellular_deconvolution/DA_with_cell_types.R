library(here)
library(MotrpacHumanPreSuspensionAnalysis)
library(MotrpacHumanPreSuspensionData)
library(BiocParallel)
library(dplyr)
library(variancePartition)
library(reformulas)

# First which cell types are we using + how are we handling them? Things I'd like to avoid:

# Any cells with observations with 0s ->
# The low abundance and 0 meausrements makes me think that the values are not going
# to be estimated correctly for some of the participants. Even if we scale it,
# the model distribution would be affected by inclusion of 0s.
# so here I include

config = jsonlite::fromJSON("~/config.json")
gitdir = config$gitdir

analysis_data_raw = file.path(gitdir,
                              "MotrpacHumanPreSuspensionAnalysis", "data-raw")
differential_analysis_scripts = list.files(file.path(analysis_data_raw, "generate_differential_analysis"),
                                           full.names = TRUE)
#dont include the clinical/sex DA because that was implemented slightly differently.
differential_analysis_scripts = c(differential_analysis_scripts, file.path(analysis_data_raw, "gsutil_path_parsing.R"))
differential_analysis_scripts = differential_analysis_scripts[!grepl("clini|sex_", differential_analysis_scripts)]

lapply(differential_analysis_scripts, source)

# we start with just blood transcriptomics with the extra stuff
#--------------------------------

desired_ome = 'transcript-rna-seq'; tissue = 'blood'
local_path = file.path(gitdir, "precovid-analyses", "data/tmp/") #just for output
# counts_data_path = .find_path_name(desired_ome = desired_ome, tissue = tissue, data_type = "rsem-genes-count")
counts_data_path = "gs://motrpac-data-hub/quant-id/human-precovid/v1.0/transcriptomics/t04-blood-rna/transcript-rna-seq/motrpac_human-precovid_t04-blood-rna_transcript-rna-seq_rsem-genes-count_v1.0.txt"
raw_counts_input = MotrpacBicQC::dl_read_gcp(counts_data_path, sep = '\t', tmpdir = local_path)
parsed_qc_norm = load_qc(selected_omes = desired_ome,
                         selected_tissues = tissue,
                         load_acute_only = TRUE)
#-> so here we want to make sure the genes we chose and the participants we chose are the same as the ones that pass through expression logcpm cutoffs
metadata = parsed_qc_norm[[tissue]][[desired_ome]][['sample_metadata']]

# sapply(cell_type_metadata, min) -> use to check if there's 0s
cell_type_cols = c("T.cells", "B.cells", "NK.cells", "Monocytes", "Neutrophils")

cell_type_metadata = read.csv(file.path(here(),
                                        "revisions",
                                        "landscape",
                                        "cellular_deconvolution",
                                        "cibersort_with_metadata.csv")) %>%
  dplyr::select(vialLabel, dplyr::all_of(cell_type_cols)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(cell_type_cols), scale))

full_metadata = metadata %>%
  right_join(cell_type_metadata, by = "vialLabel")
rownames(full_metadata) = full_metadata$vialLabel

raw_counts_input = raw_counts_input %>%
  dplyr::filter(gene_id %in% rownames(parsed_qc_norm[[tissue]][[desired_ome]][['qc_norm']]))  %>% #filter to select genes
  tibble::column_to_rownames("gene_id") %>% #set to rownames
  dplyr::select(as.character(full_metadata$vialLabel)) #reorganize and filter to only desired participants

pre_dge <- edgeR::DGEList(counts = raw_counts_input) #standard dge processing
pre_dge <- edgeR::calcNormFactors(pre_dge) #standard dge processing

#here i add the custom covariates

added_covariates = COVARIATES_FILE %>%
  rbind(c("transcript-rna-seq", "Technical", "numerical", "T.cells", "blood")) %>%
  rbind(c("transcript-rna-seq", "Technical", "numerical", "B.cells", "blood")) %>%
  rbind(c("transcript-rna-seq", "Technical", "numerical", "NK.cells", "blood")) %>%
  rbind(c("transcript-rna-seq", "Technical", "numerical", "Monocytes", "blood")) %>%
  rbind(c("transcript-rna-seq", "Technical", "numerical", "Neutrophils", "blood"))

process_metadata = process_covariates(meta = full_metadata,
                                      selected_ome = desired_ome,
                                      tissue_input = tissue,
                                      custom_covariates = added_covariates)

.run_models(repo_local_dir = local_path,
            model_type = "acute",
            expression_object = pre_dge,
            process_metadata = process_metadata,
            tissue = tissue,
            ome = desired_ome,
            voom = TRUE,
            parallel = FALSE)
