
.convert_sex_diff_output = function(fit,
                                    metadata = NULL,
                                    formula = NULL,
                                    tissue = NULL,
                                    ome = NULL,
                                    output_dir = file.path(here(), "data/tmp/sex_differences_baseline/")){

  #subset to only the models with some type of actual contrast
  res_tissue = data.frame() #final output file
  contrast = "SexMale"
  res = variancePartition::topTable(fit, coef = contrast, number = Inf, p.value = 1, confint = TRUE)
  res_single = res %>%
    dplyr::mutate(
      logLik = fit$logLik,
      feature_id = rownames(.),
      assay = ome,
      contrast = "SexMale",
      full_model = formula,
      tissue = tissue
    ) %>%
    dplyr::rename(
      p_value = P.Value,
      adj_p_value = adj.P.Val
    ) %>%
    dplyr::select(assay,
                  feature_id,
                  logFC,
                  CI.L, CI.R,
                  logLik,
                  t,
                  AveExpr,
                  p_value, adj_p_value,
                  contrast,
                  full_model)

  #----here we now add the n per participant group the relevant metadata for
  res_tissue = rbind(res_tissue, res_single)
  res_tissue = res_tissue %>% dplyr::arrange(adj_p_value)

  if(grepl("metab", ome)){
    #because metab doesn't have direct gene connection so we just output what we have.
    write_with_path_name(res_tissue,
                         local_path = output_dir,
                         ome = ome,
                         tissue = tissue,
                         data_category = "DA",
                         data_details = "sex_diff_baseline_simple",
                         version = "1.31")
    return()
  }

  sex_diff_table = res_tissue %>%
    dplyr::select(-assay) %>%
    left_join(., MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE, by = "feature_id")

  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_annotation <- getBM(
    attributes = c("hgnc_symbol", "chromosome_name"),
    filters = "hgnc_symbol",
    values = sex_diff_table$gene_symbol,
    mart = ensembl
  )
  filt_gene_anno = gene_annotation %>%
    filter(chromosome_name %in% c(as.character(1:23), "X", "Y")) %>%
    distinct(hgnc_symbol, chromosome_name)

  sex_diff_table = sex_diff_table %>%
    left_join(.,filt_gene_anno, by = c("gene_symbol" = "hgnc_symbol")) %>%
    dplyr::select(feature_id, gene_symbol, chromosome_name, everything())

  #-------and then here we save all results
  write_with_path_name(sex_diff_table,
                       local_path = output_dir,
                       ome = ome,
                       tissue = tissue,
                       data_category = "DA",
                       data_details = "sex_diff_baseline_simple",
                       version = "1.31")
}

.match_ome_tissue_code = function(desired_ome, input_tissue){
  specific_output = MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
    dplyr::filter(ome == desired_ome, tissue == input_tissue)
  if(nrow(specific_output) == 0) stop("No ome matches the desired ome/tissue combo")
  tissue_name_output = specific_output[["tissue_code"]]
  return(tissue_name_output)
}

write_with_path_name = function(actual_data_object = NULL,
                                local_path = NULL,
                                ome = NULL,
                                tissue = NULL,
                                data_category = NULL,
                                data_details = NULL,
                                version = "1.2",
                                return_name_only = FALSE){
  file_type = ".txt"
  all_file_header = "human-precovid-sed-adu" #this is the base structure for all files within the phase.
  tissue_code = .match_ome_tissue_code(desired_ome = ome, input_tissue = tissue)
  file_name = paste(all_file_header, tissue_code, ome, data_category, data_details, sep = "_")
  file_name = paste0(local_path, file_name, "_v", version, file_type)
  if(return_name_only) {
    return(file_name)
  }else{
    write.table(actual_data_object, file = file_name, row.names = F, sep = '\t', quote = F)
  }
}


.run_sex_diff_dream = function(ome_data, tissue, desired_ome) {
  metadata = ome_data[[tissue]][[desired_ome]][["sample_metadata"]] %>%
    dplyr::filter(Timepoint == "pre_exercise")
  rownames(metadata) = metadata$vialLabel

  data_matrix = ome_data[[tissue]][[desired_ome]][["qc_norm"]] %>%
    dplyr::select(as.character(metadata$vialLabel))
  if (nrow(data_matrix) == 0) return(invisible(NULL))

  process_metadata = process_covariates(
    meta = metadata,
    selected_ome = desired_ome,
    tissue_input = tissue
  )
  process_metadata$metadata$Sex = relevel(process_metadata$metadata$Sex, ref = "Female")

  full_formula = "~ Sex"
  fit = variancePartition::dream(data_matrix, full_formula, process_metadata$metadata)
  fit = variancePartition::eBayes(fit)
  .convert_sex_diff_output(fit,
                           metadata = process_metadata$metadata,
                           formula = full_formula,
                           tissue = tissue,
                           ome = desired_ome
  )
}



process_covariates = function(meta,
                              selected_ome,
                              tissue_input,
                              include_technical = TRUE,
                              custom_covariates = NULL){
  covariates_return = list() #output list for the end part
  covariates_return[["original_meta"]] = meta

  if(!is.null(custom_covariates)){
    input_covariates = custom_covariates
  }else{
    input_covariates = COVARIATES_FILE
  }
  covariates = input_covariates %>%
    as.data.frame() %>%
    dplyr::filter(ome == selected_ome) %>%
    dplyr::filter(tissue == 'all' | tissue == tissue_input)

  num_cov = covariates %>% dplyr::filter(data_type == "numerical") #numerical covariates
  factor_cov = covariates %>% dplyr::filter(data_type == "factor")

  sel_meta = meta %>%
    dplyr::select(all_of(covariates$covariate)) %>%
    dplyr::mutate(across(all_of(num_cov$covariate), ~ scale(.) %>% as.numeric())) %>%
    dplyr::mutate(across(all_of(factor_cov$covariate), ~ as.factor(.) %>% droplevels())) %>%
    dplyr::mutate(group_timepoint = droplevels(interaction(randomGroupCode, Timepoint))) %>%
    dplyr::mutate(visit_group_timepoint = droplevels(interaction(visitcode, randomGroupCode, Timepoint)))

  technical_covs = covariates %>% filter(tech_or_design == "Technical")
  full_formula = names(sel_meta)[!names(sel_meta) %in% c("randomGroupCode", "Timepoint", "visitcode", "pid", "group_timepoint", "visit_group_timepoint")] #remove these from the character vector
  #the purpose of the design covariates section is to make a model.matrix()
  design_covs = c(full_formula[!full_formula %in% as.character(technical_covs$covariate)], "group_timepoint")
  if(!include_technical){ #remove for any modeling where some covariates have been regressed out
    full_formula = full_formula[!full_formula %in% technical_covs$covariate]
  }
  formula_string = paste(full_formula, collapse = " + ")
  #we readd group_timepoint first because of the way some contrast matrixes drop values in case of interactions w other levels of factors in the contrast matrixes
  formula_string_full = paste("~ 0 + group_timepoint + ", formula_string, "+ (1 | pid)")
  non_mixed_model = paste("~ 0 + group_timepoint + ", formula_string)

  #so i remove group_timepoint above and then make sure that it comes first because the order of the string can sometimes
  #actually change the contrast matrix formed and which columns are dropped in terms of the contrast comparisons.
  covariates_return[["technical_cov"]] = technical_covs
  covariates_return[["design_cov"]] = design_covs
  covariates_return[["full_formula"]] = formula_string_full
  covariates_return[["metadata"]] = sel_meta #so this is with all the tech/num cov in the correct format

  return(covariates_return)
}
