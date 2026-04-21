.generate_contrasts_acute = function(metadata){
  pre_contrast_expressions <- c()
  timepoints = unique(metadata$Timepoint)
  for(tp in timepoints){
    # message(tp)
    if (!tp == 'pre_exercise'){
      # meta_tp = metadata %>% filter(Timepoint == tp)
      # Chris: note - should just make contrast_Endur_Cntrl = paste0(contrast_Endur, " - ", contrast_Cntrls) at some point. This could definitely be refactored...
      contrast_Endur_Cntrl = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
      contrast_Resist_Cntrl = sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
      contrast_Endur_Resist = sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUResist.%s + group_timepointADUResist.pre_exercise", tp, tp)
      if (tp == 'during_20_min' | tp == 'during_40_min') {contrast_Resist_Cntrl = NULL; contrast_Endur_Resist = NULL; contrast_Resist = NULL} #resistance group doesn't get blood draws here
      pre_contrast_expressions = c(pre_contrast_expressions,
                                   contrast_Endur_Cntrl,
                                   contrast_Resist_Cntrl,
                                   contrast_Endur_Resist)
    }
  }
  return(pre_contrast_expressions)
}


plot_cbc_vs_cibersortx = function(blood_metadata, cbc_raw, filter_lt5 = TRUE) {
  cbc = cbc_raw %>%
    dplyr::mutate(
      Neutrophils = neutro_labr,
      Monocytes = mono_labr,
      Lymphocytes = lymp_labr
    ) %>%
    dplyr::select(pid, Neutrophils, Monocytes, Lymphocytes) %>%
    dplyr::mutate(pid = as.character(pid))

  if (filter_lt5) {
    cbc = cbc %>% dplyr::filter(Neutrophils > 5)
  }

  neutrophil_corr = blood_metadata %>%
    dplyr::filter(Timepoint == "pre_exercise") %>%
    dplyr::mutate(pid = as.character(pid)) %>%
    dplyr::select(pid, Neutrophils) %>%
    dplyr::rename(deconv_neutrophils = Neutrophils) %>%
    dplyr::inner_join(cbc %>% dplyr::select(pid, Neutrophils), by = "pid") %>%
    dplyr::rename(cbc_neutrophils = Neutrophils)

  r_val = cor(neutrophil_corr$deconv_neutrophils, neutrophil_corr$cbc_neutrophils,
              use = "complete.obs")

  filter_label = if (filter_lt5) "excluding 2 erroneous samples" else "all samples"

  neutrophil_corr %>%
    ggplot(aes(x = cbc_neutrophils, y = deconv_neutrophils)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5,
             label = paste0("r = ", round(r_val, 2)), size = 4) +
    labs(
      x = "CBC Neutrophils",
      y = "CIBERSORTx Neutrophils (%)",
      title = "Neutrophil estimates: CIBERSORTx vs. CBC (pre-exercise)",
      subtitle = filter_label
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      title = element_text(size = 13)
    )
}


squish_proportions = function(merged_proportions){
  #this cell type map should probably be turned into an excel
  #or sometype of dataframe or something in case we start getting
  #multiple
  cell_type_map = c(
    # T cells
    "T.cells.CD4.memory.activated" = "T cells",
    "T.cells.CD4.memory.resting" = "T cells",
    "T.cells.CD4.naive" = "T cells",
    "T.cells.CD8" = "T cells",
    "T.cells.follicular.helper" = "T cells",
    "T.cells.gamma.delta" = "T cells",
    "T.cells.regulatory..Tregs." = "T cells",

    # B cells
    "B.cells.memory" = "B cells",
    "B.cells.naive" = "B cells",
    "Plasma.cells" = "B cells",

    # NK cells
    "NK.cells.activated" = "NK cells",
    "NK.cells.resting" = "NK cells",

    # Monocytes
    "Monocytes" = "Monocytes",

    # Macrophages
    "Macrophages.M0" = "Macrophages",
    "Macrophages.M1" = "Macrophages",
    "Macrophages.M2" = "Macrophages",

    # Dendritic cells
    "Dendritic.cells.activated" = "Dendritic cells",
    "Dendritic.cells.resting" = "Dendritic cells",

    # Neutrophils
    #just to get the 100*...proportion stuff
    "Neutrophils" = "Neutrophils",

    # Mast cells
    "Mast.cells.activated" = "Mast cells",
    "Mast.cells.resting" = "Mast cells"#,
  )

  unique_groups = unique(unname(cell_type_map))
  meta_cols = setdiff(colnames(merged_proportions), names(cell_type_map))

  group_sums = purrr::map_dfc(unique_groups, function(grp) {
    cols = intersect(names(cell_type_map[cell_type_map == grp]), colnames(merged_proportions))
    tibble::tibble(!!grp := rowSums(dplyr::select(merged_proportions, dplyr::all_of(cols)), na.rm = TRUE))
  })

  summed_grouped = dplyr::bind_cols(merged_proportions[, meta_cols], group_sums) %>%
    pivot_longer(cols = dplyr::all_of(unique_groups),
                 names_to = "cell_group",
                 values_to = "proportion") %>%
    group_by(Mixture) %>%
    mutate(proportion = 100 * proportion / sum(abs(proportion))) %>%
    ungroup()

  #----------------------------------------------------
  #we check if there are any cells that are on avg < 1% of cells in general.
  #we just remove those from visualization because thats not really helpful
  cell_group_means = summed_grouped %>%
    group_by(cell_group) %>%
    summarise(cell_group_avg = mean(proportion), .groups = "drop") %>%
    filter(cell_group_avg > 1)

  summed_grouped = summed_grouped %>%
    filter(cell_group %in% cell_group_means$cell_group) %>%
    pivot_wider(names_from = "cell_group", values_from = "proportion")
  #----------------------------------------------------

  return(summed_grouped)
}
