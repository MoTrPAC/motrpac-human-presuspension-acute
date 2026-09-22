#!/usr/bin/env Rscript
# Extended Data 2 — Differential abundance supplement
#
# Panels:  ED2A  signed differential-abundance hit types by ome and tissue
#          ED2B  example single-feature trajectories per hit type
#          ED2C  cross-tissue overlap of differential features, one panel per ome
#          ED2D  NR4A1 and NFKBIZ transcript trajectories, all three tissues
#          ED2E  kinase and perturbation signatures enriched in the shared
#                phosphosites
#          ED2F  overlap between omes within each tissue
# Tables:  ST2f   phospho ORA (PTMsigDB)                    (written by ED2E)
#
# Needs consortium data access. ED2A, ED2C, ED2E and ED2F share one
# load_differential_analysis(epigen = TRUE) pass; the ATAC and methylCap DA
# tables are not in the data package and are downloaded into EPIGEN_QC_DIR on
# first use, then reused.
#
#   Rscript figures/landscape/ED2.R          every panel
#   Rscript figures/landscape/ED2.R ED2D     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(dplyr)
  library(ggh4x)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "da_overlap_helpers.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))

# ---- shared ----------------------------------------------------------------

# One load_differential_analysis(epigen = TRUE) pass behind ED2A, ED2C, ED2E and
# ED2F. Lazy and memoised rather than loaded at the top of the script: 7.7 GB
# should not be read to build a figure whose selected panels do not need it.
#
# Every contrast type, because ED2A cross-classifies a feature on three of them.
# cross_tissue_da() in lib/da_overlap_helpers.R makes this same load and then
# keeps only exercise_with_controls, so calling it here as well would read the
# tables a second time; exercise_da() below applies its filter to the load
# already in hand instead.
acute_da <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    cache <<- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = "all",
      selected_tissues = "all",
      epigen = TRUE,
      combine_with_featgene = TRUE,
      single_matrix = TRUE,
      repo_local_dir = epigen_qc_dir(),
      verbose = FALSE
    )

    # epigen = TRUE is a request, not a guarantee. load_DA_from_bucket() lists
    # the staging bucket with `gsutil ls -R` and, if that listing fails — an
    # expired credential is the common case — returns an empty list after a
    # message and nothing else, and the five package-shipped omes come back as
    # though nothing had happened. The same check cross_tissue_da() makes,
    # stated here because this is where the load is.
    missing_epigen <- setdiff(c("epigen-atac-seq", "epigen-methylcap-seq"),
                              unique(as.character(cache$assay)))
    if (length(missing_epigen) > 0) {
      stop("epigen = TRUE but the differential analysis came back without ",
           paste(missing_epigen, collapse = " and "),
           ".\n  The staging-bucket listing failed silently. Refresh the ",
           "credential with `gcloud auth login` and rebuild; the cached files ",
           "under EPIGEN_QC_DIR cannot be used without it.", call. = FALSE)
    }
    cache
  }
})

# What cross_tissue_da() returns, off the load above: contrast_short cut back to
# its left side, exercise_with_controls only. ED2C, ED2E and ED2F start here.
exercise_da <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- acute_da() %>%
        mutate(contrast_short = sub(" - .*", "", contrast_short)) %>%
        dplyr::filter(contrast_type == "exercise_with_controls")
    }
    cache
  }
})

# The FDR < 0.05 subset with the display columns the overlap panels group on.
# ED2C and ED2F use it whole; ED2E takes its prot-ph rows.
significant_features <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- significant_da(exercise_da())
    cache
  }
})

# ---- ED2A — signed differential-abundance hit types by ome and tissue -------

# Every feature x timepoint is cross-classified on three contrasts — exercise
# versus control, exercise alone, control alone — into a control-driven exercise
# effect, an exercise effect, or a non-exercise-induced change. Counts are drawn
# above zero for up-regulated features and below zero for down-regulated ones,
# with grey bars marking tissue x ome x timepoint cells that were never assayed.

ed2a <- function() {
  panel_init("ED2A")

  # Every ome but methylation, which this panel does not draw. The shared load
  # carries methylCap because ED2C and ED2F do.
  full_da_results_all_tissues <- acute_da() %>%
    dplyr::filter(assay != "epigen-methylcap-seq") %>%
    # All metabolomics platforms count as one ome, and Olink joins mass-spec
    # proteomics.
    mutate(assay = case_when(
      grepl("metab", assay) ~ "Metabolomics",
      grepl("prot-ol|prot-pr", assay) ~ "Proteomics",
      TRUE ~ assay
    )) %>%
    mutate(assay = recode(
      assay,
      "transcript-rna-seq" = "Transcriptomics",
      "prot-ph" = "Phosphoproteomics",
      "epigen-atac-seq" = "ATAC"
    ))

  # The six post-baseline timepoints this panel draws, in x-axis order, with the
  # short labels the axis carries. Both the counted data and the "not assayed"
  # grid are relabelled through this, so the two always agree.
  timepoint_labels <- c(
    "during_20_min"     = "D20M",
    "during_40_min"     = "D40M",
    "post_10_min"       = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr"     = "P3.5/4H",
    "post_24_hr"        = "P24H"
  )

  # Anything the map does not mention is left as it is.
  abbreviate_timepoint <- function(x) {
    x <- as.character(x)
    out <- unname(timepoint_labels[x])
    unmapped <- is.na(out)
    out[unmapped] <- x[unmapped]
    out
  }

  # One row per feature x timepoint x group with the three contrast p-values side
  # by side. The control-only p-value is shared across the two exercise groups, so
  # it is filled in from whichever row of the feature x timepoint block carries it.
  wide_by_contrasts <- full_da_results_all_tissues %>%
    dplyr::filter(!contrast_type %in% c("Endur_vs_Resist", "baseline")) %>%
    dplyr::mutate(direction = sign(logFC)) %>%
    dplyr::select(tissue, assay, randomGroupCode, Timepoint, contrast_type,
                  adj_p_value, feature_id, direction) %>%
    pivot_wider(names_from = c("contrast_type"), values_from = adj_p_value) %>%
    group_by(tissue, assay, Timepoint, feature_id) %>%
    mutate(control_only = ifelse(
      is.na(control_only),
      dplyr::first(na.omit(control_only)),
      control_only
    )) %>%
    ungroup() %>%
    dplyr::filter(!randomGroupCode == "Control") %>%
    mutate(Timepoint = factor(Timepoint, levels = names(timepoint_labels)))

  # Control effect: the feature is stable in the exercise arm and the comparison
  #   is significant only because the controls moved. A study without controls
  #   could not detect it.
  # Exercise effect: the exercise arm moves and the comparison against controls is
  #   significant, whichever way the controls went.
  # Non-exercise induced change: not a significant exercise-versus-control
  #   difference. A study without controls would call these hits.
  wide_by_contrasts_merge <- wide_by_contrasts %>%
    mutate(
      exercise_with_controls_cat = ifelse(exercise_with_controls <= 0.05, "Sig", "Non-sig"),
      exercise_no_controls_cat = ifelse(exercise_no_controls <= 0.05, "Sig", "Non-sig"),
      control_only_cat = ifelse(control_only <= 0.05, "Sig", "Non-sig")
    ) %>%
    mutate(
      category = case_when(
        (exercise_no_controls_cat == "Non-sig" & exercise_with_controls_cat == "Sig") ~ "Ctrl-driven exercise effects",
        (exercise_no_controls_cat == "Sig" & exercise_with_controls_cat == "Sig") ~ "Exercise effects",
        (control_only_cat == "Sig" & exercise_no_controls_cat == "Sig" & exercise_with_controls_cat == "Non-sig") ~ "Non-exercise induced changes",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::filter(category != "Other")

  # ---- plot ----

  category_colors <- c(
    "Non-exercise induced changes" = "#7570b3",
    "Exercise effects" = "brown4",
    "Ctrl-driven exercise effects" = "brown1"
  )

  wide_by_contrasts_merge$assay <- as.factor(wide_by_contrasts_merge$assay)
  wide_by_contrasts_merge$tissue <- as.factor(wide_by_contrasts_merge$tissue)
  wide_by_contrasts_merge$randomGroupCode <- as.factor(wide_by_contrasts_merge$randomGroupCode)

  # Counts per facet cell, signed so up-regulated features stack upward and
  # down-regulated ones downward.
  filtered_data <- wide_by_contrasts_merge %>%
    mutate(
      randomGroupCode = recode(randomGroupCode, "ADUEndur" = "EE", "ADUResist" = "RE"),
      Timepoint = abbreviate_timepoint(Timepoint),
      Timepoint = factor(Timepoint, levels = unname(timepoint_labels))
    ) %>%
    group_by(assay, tissue, randomGroupCode, Timepoint, category, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(y_value = ifelse(direction == 1, count, -count))

  # Row height for the grey "not assayed" bars: the tallest real bar in that ome.
  max_counts_by_assay <- filtered_data %>%
    summarise(max_count = max(abs(y_value)), .by = "assay")

  # Which tissue x ome x group x timepoint cells were run at all. Taken from the
  # unfiltered results so that a cell with no significant features is not confused
  # with a cell that was never assayed.
  existing_combinations <- full_da_results_all_tissues %>%
    distinct(tissue, assay, randomGroupCode, Timepoint) %>%
    mutate(
      Timepoint = abbreviate_timepoint(Timepoint),
      randomGroupCode = recode(randomGroupCode, "ADUEndur" = "EE", "ADUResist" = "RE")
    )

  # Blood carries every timepoint; muscle and adipose only the later ones.
  all_combinations <- rbind(
    expand.grid(
      Timepoint = unique(filtered_data$Timepoint[filtered_data$tissue == "blood"]),
      assay = unique(filtered_data$assay),
      tissue = "blood",
      randomGroupCode = unique(filtered_data$randomGroupCode)
    ),
    expand.grid(
      Timepoint = unique(filtered_data$Timepoint[filtered_data$tissue %in% c("muscle", "adipose")]),
      assay = unique(filtered_data$assay),
      tissue = c("muscle", "adipose"),
      randomGroupCode = unique(filtered_data$randomGroupCode)
    )
  )

  # The cells that need a grey background, one bar above and one below zero.
  missing_combinations <- dplyr::anti_join(
    all_combinations, existing_combinations,
    by = c("Timepoint", "assay", "tissue", "randomGroupCode")
  ) %>%
    left_join(max_counts_by_assay, by = "assay") %>%
    mutate(category = "Not assayed") %>%
    filter(!(tissue %in% c("adipose", "muscle") & Timepoint %in% c("D20M", "D40M", "P10M"))) %>%
    mutate(y_value = max_count) %>%
    tidyr::uncount(2, .id = "sign") %>%
    mutate(y_value = ifelse(sign == 1, y_value, -y_value))

  # ---- strip colors ----
  #
  # Looked up from the levels of the columns the facets are drawn from, which are
  # factors by this point. The x strips nest tissue over exercise group, so the
  # fills run tissue-by-tissue first and then group-by-group within each tissue;
  # the y strips are the omes. A name the palette does not carry is an error
  # rather than a silent grey strip — a renamed ome would otherwise lose its color
  # and nothing would say so.
  tissue_levels <- levels(filtered_data$tissue)
  assay_levels <- levels(filtered_data$assay)
  group_levels <- levels(filtered_data$randomGroupCode)

  lookup_colors <- function(palette, keys, what) {
    hits <- palette[keys]
    if (anyNA(hits)) {
      stop(what, " not in the package palette: ",
           paste(keys[is.na(hits)], collapse = ", "), call. = FALSE)
    }
    unname(hits)
  }

  tissue_fills <- lookup_colors(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS, tissue_levels, "tissue(s)")
  assay_fills <- lookup_colors(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS, assay_levels, "ome(s)")
  group_fills <- lookup_colors(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS,
    c("ADUEndur", "ADUResist")[seq_along(group_levels)], "exercise group(s)")

  # One fill per x strip, outer layer (tissue) first, then the group strips nested
  # under each tissue.
  x_fills <- c(tissue_fills, rep(group_fills, times = length(tissue_levels)))

  # Black or white label text, whichever reads against the fill behind it.
  text_on <- function(fills) {
    rgb_mat <- grDevices::col2rgb(fills) / 255
    luminance <- 0.299 * rgb_mat[1, ] + 0.587 * rgb_mat[2, ] + 0.114 * rgb_mat[3, ]
    ifelse(luminance < 0.55, "white", "black")
  }

  p <- ggplot(filtered_data, aes(x = Timepoint)) +
    geom_col(
      data = missing_combinations,
      aes(y = y_value),
      fill = "grey80", width = 0.8, alpha = 0.7
    ) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    geom_col(
      data = subset(filtered_data, category != "Non-exercise induced changes"),
      aes(y = y_value, fill = category),
      width = 0.3, position = "stack"
    ) +
    geom_col(
      data = subset(filtered_data, category == "Non-exercise induced changes"),
      aes(y = y_value, fill = category),
      width = 0.3, position = position_nudge(x = 0.3)
    ) +
    scale_fill_manual(values = c(category_colors, "Not assayed" = "grey80")) +
    coord_cartesian(clip = "off") +
    facet_grid2(
      assay ~ tissue + randomGroupCode,
      scales = "free",
      space = "free_x",
      strip = strip_nested(
        background_y = elem_list_rect(fill = assay_fills),
        text_y = elem_list_text(colour = text_on(assay_fills)),
        background_x = elem_list_rect(fill = x_fills),
        text_x = elem_list_text(colour = text_on(x_fills)),
        by_layer_x = FALSE,
        by_layer_y = FALSE
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.justification = "top",
      legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      legend.margin = margin(3, 3, 3, 3)
    )

  export_panel(p, "ED2A")
}

# ---- ED2B — example single-feature trajectories per hit type ----------------

# One representative feature per differential-abundance hit type, each drawn as
# the exercise-group mean with a 95% confidence interval across the acute
# timepoints. A point is filled black where that timepoint's adjusted p value is
# below 0.05. The four sub-plots sit two per row under one collected legend.

ed2b <- function() {
  panel_init("ED2B")
  export_panel(single_feature_plot("ED2B"), "ED2B")
}

# ---- ED2C — cross-tissue overlap of differential features, per ome ----------

# FIG2B's UpSet split by platform: six panels, each the three-tissue
# intersection of the features that ome called differential, each bar filled
# with that ome's own colour. It answers what FIG2B's stacked bars only imply,
# which is whether the cross-tissue overlap is driven by one platform.
#
# Laid out in OME_LEVELS order — ATAC, metabolomics, methylation,
# phosphoproteomics, proteomics, transcriptomics — so an ome sits in the same
# position here as it does in the stacked bars of FIG2B and ED4A.

ed2c <- function() {
  panel_init("ED2C")

  sig <- significant_features()

  TISSUE_SETS <- c("Muscle", "Blood", "Adipose")
  presence <- feature_presence_matrix(sig, by = "Tissue", sets = TISSUE_SETS)

  # The stripe colours are matched to TISSUE_SETS by position, and sort_sets is
  # off, so the row order is the one stated above rather than the size order the
  # data happens to produce. The legacy passed a fixed colour vector alongside
  # sort_sets = "ascending" and carried four different hardcoded orders across the
  # six sub-panels — the stripes were tuned by eye against one DA release.
  stripe_colors <- unname(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[tolower(TISSUE_SETS)]
  )

  # Title and bar colour per ome. The colour is looked up by the ome's own display
  # name, so it is the same colour FIG2B stacks that ome with.
  OME_PANELS <- list(
    list(ome = "Chromatin Accessibility (ATAC)",
         title = "Differentially Accessible Peaks",
         count_label = "Number of differential \npeaks"),
    list(ome = "Metabolomics",
         title = "DA metabolites",
         count_label = "Number of \ndifferential metabolites"),
    list(ome = "Methylation",
         title = "DA methylated regions",
         count_label = "Number of differential \nregions"),
    list(ome = "Phosphoproteomics",
         title = "DA Phosphosites",
         count_label = "Number of differential \nphosphosites"),
    list(ome = "Proteomics",
         title = "DA Proteins",
         count_label = "Number of differential \nproteins"),
    list(ome = "Transcriptomics",
         title = "DA transcripts",
         count_label = "Number of differential \ntranscripts")
  )

  ome_upset <- function(panel) {
    rows <- presence[presence$Ome == panel$ome, c(TISSUE_SETS), drop = FALSE]
    if (nrow(rows) == 0) {
      stop("no significant features for ", panel$ome, call. = FALSE)
    }

    fill <- unname(
      MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS[panel$ome]
    )
    if (is.na(fill)) {
      stop("ome not in HUMAN_OME_COLORS: ", panel$ome, call. = FALSE)
    }

    ComplexUpset::upset(
      rows,
      TISSUE_SETS,
      encode_sets = FALSE,
      width_ratio = 0.2,
      height_ratio = 1,
      set_sizes = FALSE,
      base_annotations = stats::setNames(
        list(intersection_size(counts = TRUE, text = list(size = 3), fill = fill)),
        panel$count_label
      ),
      themes = upset_modify_themes(
        list(intersections_matrix = theme(axis.title.x = element_blank()))
      ),
      sort_sets = FALSE,
      wrap = TRUE,
      stripes = stripe_colors
    ) +
      ggtitle(panel$title)
  }

  ome_plots <- lapply(OME_PANELS, ome_upset)

  combined <- patchwork::wrap_plots(ome_plots, nrow = 2)

  export_panel(function() print(combined), "ED2C")
}

# ---- ED2D — NR4A1 and NFKBIZ transcript trajectories ------------------------

# Two features, each drawn as the exercise-group mean with a 95% confidence
# interval across the acute timepoints, one facet per tissue. A point is filled
# black where that timepoint's adjusted p value is below 0.05. The two rows sit
# under one collected legend.

ed2d <- function() {
  panel_init("ED2D")
  export_panel(single_feature_plot("ED2D"), "ED2D")
}

# ---- ED2E — kinase and perturbation signatures in the shared phosphosites ----

# Over-representation against PTMsigDB of three phosphosite sets: muscle-specific,
# adipose-specific, and the sites called differential in both. The curated
# signatures are drawn as a heatmap of -log10(p), one column per set, with a grey
# cell where the enrichment does not clear the adjusted-p cut.
#
# THE PANEL IS DRAWN FROM THE VENDORED CSV, not from the computation below.
# landscape_figureED2_PTMSigDB_ORA_results_v2_Aug2026.csv is the published
# result; ST2f and the heatmap both come from it.
#
# The pass IS still run on every build, per tissue, and its input, background and
# set counts are reported against the CSV's. It is a check, not a source: nothing
# it produces reaches the figure or the table. Keeping it live means a divergence
# shows up in the build log rather than being discovered later.
#
# Site confidence in that check is read PER TISSUE, from
# *_PROT_PH_QC$feature_metadata rather than from HUMAN_FEATURE_TO_GENE. The two
# disagree by construction and the per-tissue table is the right one:
#
#   HUMAN_FEATURE_TO_GENE is keyed on (assay, feature_id) with no tissue column,
#   so its prot-ph rows are the UNION of the two tissues — 18,548 muscle plus
#   21,022 adipose over 7,865 shared, exactly 31,705. Localization confidence is
#   a per-tissue measurement, and 859 of those shared sites disagree between
#   muscle and adipose (782 muscle-only, 77 adipose-only). Having one row per
#   site, that table collapses every disagreement to FALSE, which undercounts
#   muscle's confident sites by up to 782. flanking_sequence does NOT diverge —
#   all 7,865 agree — because it is a property of the protein, not the assay.
#
# Per tissue, the check lands within 14 of the CSV on muscle's input (3,850 vs
# 3,836) but 572 over on its background, and 9 of the 48 drawn cells still cross
# the adjusted-p cut. Two questions are open for the follow-up: the exact
# confidence definition, and the CSV's adipose-and-muscle background of 15,818,
# which is the muscle background rather than the two-tissue intersection.
#
# PTMsigDB sets are keyed on 15-mer flanking sequences carrying a ";u"/";d"
# direction suffix. run_ORA() strips that suffix itself when every requested
# database is PTMSIGDB, so the input and background are bare flanking sequences —
# and the database argument has to be PTMSIGDB alone, or that branch never fires.
#
# This panel also writes ST2f.

ed2e <- function() {
  spec <- panel_init("ED2E")

  # ---- data ----

  all_da <- exercise_da()
  sig <- significant_features()

  phospho_da <- all_da |> filter(.data$assay == "prot-ph")
  phospho_sig <- sig |> filter(.data$assay == "prot-ph")

  presence <- feature_presence_matrix(phospho_sig, by = "Tissue",
                                      sets = c("Adipose", "Muscle"))

  # ---- per-tissue site confidence ----

  # One lookup per tissue: feature_id -> flanking sequence, restricted to the sites
  # that tissue localized confidently.
  confident_flanks <- lapply(c(adipose = "adipose", muscle = "muscle"), function(tis) {
    meta <- MotrpacHumanPreSuspensionData::load_qc(
      selected_tissues = tis, selected_omes = "prot-ph"
    )[[tis]][["prot-ph"]]$feature_metadata

    for (needed in c("feature_id", "flanking_sequence", "confident_site")) {
      if (!needed %in% names(meta)) {
        stop(tis, " prot-ph feature_metadata has no ", needed, " column",
             call. = FALSE)
      }
    }
    keep <- which(as.logical(meta$confident_site))
    stats::setNames(as.character(meta$flanking_sequence)[keep],
                    as.character(meta$feature_id)[keep])
  })

  flanking_of <- function(feature_ids, tissue) {
    out <- confident_flanks[[tissue]][as.character(feature_ids)]
    sort(unique(out[!is.na(out) & nzchar(out)]))
  }

  tissue_background <- lapply(c(adipose = "adipose", muscle = "muscle"), function(tis) {
    flanking_of(phospho_da$feature_id[as.character(phospho_da$tissue) == tis], tis)
  })

  # The membership sets. A site is only counted for a tissue where THAT tissue
  # localized it confidently, so a set spanning both tissues takes the sites both
  # agree on and is tested against the intersected background: a site can only be
  # called shared if it was confidently measured in both.
  PHOSPHO_ORA_SETS <- list(
    Muscle_only    = list(present = "Muscle",  absent = "Adipose",
                          tissues = "muscle"),
    Adipose_only   = list(present = "Adipose", absent = "Muscle",
                          tissues = "adipose"),
    Adipose_Muscle = list(present = c("Adipose", "Muscle"), absent = character(0),
                          tissues = c("adipose", "muscle"))
  )

  ora_results <- list()
  for (set_name in names(PHOSPHO_ORA_SETS)) {
    definition <- PHOSPHO_ORA_SETS[[set_name]]

    keep <- rep(TRUE, nrow(presence))
    for (tis in definition$present) keep <- keep & presence[[tis]] == 1
    for (tis in definition$absent)  keep <- keep & presence[[tis]] == 0
    member_ids <- presence$overlap_id[keep]

    input <- Reduce(intersect, lapply(definition$tissues, function(tis) {
      flanking_of(member_ids, tis)
    }))
    background <- Reduce(intersect, tissue_background[definition$tissues])

    stray <- setdiff(input, background)
    if (length(stray) > 0) {
      stop(set_name, ": ", length(stray), " foreground site(s) outside their own ",
           "background", call. = FALSE)
    }
    if (length(input) == 0) {
      stop(set_name, ": no confidently localized sites in the membership set",
           call. = FALSE)
    }

    ora_results[[set_name]] <- MotrpacHumanPreSuspensionAnalysis::run_ORA(
      input = input,
      background = background,
      database = "PTMSIGDB"
    )
  }

  computed_ora <- dplyr::bind_rows(ora_results, .id = "contrast")
  computed_ora$contrast <- factor(computed_ora$contrast,
                                  levels = names(PHOSPHO_ORA_SETS))

  # ---- the published result, which is what the panel uses ----

  ptmsigdb_csv <- panel_source(
    spec$figure_dir,
    "landscape_figureED2_PTMSigDB_ORA_results_v2_Aug2026.csv",
    env_var = "ED2E_PTMSIGDB_ORA_CSV"
  )

  # Comma-delimited, unlike most MoTrPAC files: an ORA result table rather than a
  # freeze file, with no sample-id columns for name repair to mangle.
  ora <- utils::read.csv(ptmsigdb_csv, check.names = FALSE, stringsAsFactors = FALSE)

  required <- c("contrast", "collection", "database", "set_id", "set", "set_short",
                "set_size", "set_size_DB", "size_ratio", "set_size_in_input",
                "input_size", "background_size", "p_value", "adj_p_value")
  missing_cols <- setdiff(required, names(ora))
  if (length(missing_cols) > 0) {
    stop("the PTMsigDB ORA table is missing column(s): ",
         paste(missing_cols, collapse = ", "),
         "\n  Read from: ", ptmsigdb_csv, call. = FALSE)
  }
  unexpected <- setdiff(unique(ora$contrast), names(PHOSPHO_ORA_SETS))
  if (length(unexpected) > 0) {
    stop("the ORA table carries contrasts this panel does not declare: ",
         paste(sort(unexpected), collapse = ", "), call. = FALSE)
  }

  # -log10, which is what this panel's legend says, and the scale FIG2Bii and FIG2F
  # use.
  ora$statistic_column <- -log10(ora$p_value)
  ora$contrast <- factor(ora$contrast, levels = names(PHOSPHO_ORA_SETS))

  write_st2f(ora)

  # ---- what the check says, for the log only ----

  for (set_name in levels(ora$contrast)) {
    mine <- computed_ora[computed_ora$contrast == set_name, ]
    theirs <- ora[ora$contrast == set_name, ]
    if (nrow(mine) == 0 || nrow(theirs) == 0) next
    message(sprintf(
      "        %-15s recomputed input %5d (published %5d)  background %6d (%6d)  sets %4d (%4d)",
      set_name, mine$input_size[1], theirs$input_size[1],
      mine$background_size[1], theirs$background_size[1],
      nrow(mine), nrow(theirs)))
  }

  # ---- the curated selection ----

  curated_ids <- unique(highlight_pathways("ED2E")$set_id)

  plot_df <- ora |>
    filter(as.character(.data$set_id) %in% curated_ids) |>
    mutate(set_short = sub("^PTMSIGDB_", "", .data$set_short)) |>
    select("set_short", "statistic_column", "adj_p_value", "contrast")

  found_ids <- intersect(curated_ids, as.character(ora$set_id))
  missing <- setdiff(curated_ids, found_ids)
  if (length(missing) > 0) {
    message(sprintf(
      "        ED2E: %d of %d curated signature(s) not returned by this run and not drawn: %s",
      length(missing), length(curated_ids), paste(sort(missing), collapse = ", ")))
  }
  if (length(found_ids) == 0) {
    stop("none of the curated signatures is in this run's ORA", call. = FALSE)
  }
  message(sprintf("        ED2E: drawing %d of %d curated signature(s)",
                  length(found_ids), length(curated_ids)))

  plot_df <- plot_df[order(plot_df$contrast), ]
  levels(plot_df$contrast) <- gsub("_", " ", gsub("_only", "",
                                                  levels(plot_df$contrast)))

  # ---- plot ----

  # enrichmap() finishes with ComplexHeatmap::draw(), whose newpage default starts
  # a fresh page on the device export_panel() already opened, leaving page 1 blank
  # and the heatmap on page 2. draw_args carries newpage = FALSE through, the same
  # way FIG3B does.
  heatmap <- function() {
    TMSig::enrichmap(
      as.data.frame(plot_df),
      set_column = "set_short",
      n_top = max(30L, length(found_ids)),
      statistic_column = "statistic_column",
      contrast_column = "contrast",
      padj_column = "adj_p_value",
      padj_legend_title = "adj p \n(background)",
      padj_fill = "grey80",
      colors = c("white", "#543483"),
      heatmap_args = list(
        name = "-log10(p)",
        na_col = "grey90",
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        border = TRUE,
        rect_gp = grid::gpar(col = "grey70", lwd = 2)
      ),
      draw_args = list(newpage = FALSE)
    )
  }

  export_panel(heatmap, "ED2E")
}

# ---- ED2F — overlap between omes within each tissue -------------------------

# The complement of ED2C: instead of asking which tissues share a feature, it
# asks which platforms agree about a gene inside one tissue. One UpSet per
# tissue, the sets being the omes, the unit being an Ensembl gene rather than a
# feature — a differential ATAC peak and a differential transcript are the same
# gene seen twice, and that is the whole point of the panel.
#
# Adipose over blood over muscle, one column.

ed2f <- function() {
  panel_init("ED2F")

  sig <- significant_features()

  # The set columns, in the order they are drawn. Olink is kept apart from
  # mass-spec proteomics here, unlike everywhere else in Figure 2: within one
  # tissue the two platforms are separate assays of the same gene, and collapsing
  # them would erase exactly the agreement this panel measures.
  # OME_LEVELS order, with Olink kept next to mass-spec proteomics because the two
  # are the same ome measured twice and the panel is about platform agreement.
  OME_SETS <- c(
    "ATAC"              = "epigen-atac-seq",
    "Methylation"       = "epigen-methylcap-seq",
    "Phosphoproteomics" = "prot-ph",
    "Proteomics"        = "prot-pr",
    "OLINK"             = "prot-ol",
    "Transcriptomics"   = "transcript-rna-seq"
  )

  # Every gene the DA tables can speak about, which is the universe the
  # indicator columns are built over.
  gene_of <- function(feature_ids, assays) {
    key <- paste(assays, feature_ids, sep = "\r")
    map <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
    map_key <- paste(as.character(map$assay), as.character(map$feature_id),
                     sep = "\r")
    as.character(map$ensembl_gene)[match(key, map_key)]
  }

  sig$ensembl_gene <- gene_of(as.character(sig$feature_id),
                              as.character(sig$assay))

  tissue_upset <- function(tissue_name) {
    rows <- sig[sig$Tissue == tissue_name &
                  !is.na(sig$ensembl_gene) &
                  nzchar(sig$ensembl_gene), , drop = FALSE]

    genes <- sort(unique(rows$ensembl_gene))
    if (length(genes) == 0) {
      stop("no significant gene-mapped features in ", tissue_name, call. = FALSE)
    }

    indicators <- data.frame(row.names = genes)
    for (set_name in names(OME_SETS)) {
      hit <- unique(rows$ensembl_gene[rows$assay == OME_SETS[[set_name]]])
      indicators[[set_name]] <- as.integer(genes %in% hit)
    }

    # An ome with no hits in this tissue would draw an empty set row. Dropping it
    # is what the legacy did; the difference is that the stripe colours are keyed
    # by name below, so dropping a set no longer shifts every colour after it.
    present_sets <- names(OME_SETS)[colSums(indicators) > 0]
    indicators <- indicators[, present_sets, drop = FALSE]

    # The set names are this panel's own labels, and two of them are not palette
    # keys: HUMAN_OME_COLORS has "Proteomics (Olink)" and "Chromatin Accessibility
    # (ATAC)" where the sets read "OLINK" and "ATAC". Mapped rather than renamed,
    # because the short labels are what fit the UpSet matrix.
    palette_key <- c(
      "ATAC"              = "ATAC",
      "Methylation"       = "Methylation",
      "Phosphoproteomics" = "Phosphoproteomics",
      "Proteomics"        = "Proteomics",
      "OLINK"             = "Proteomics (Olink)",
      "Transcriptomics"   = "Transcriptomics"
    )
    unmapped <- setdiff(present_sets, names(palette_key))
    if (length(unmapped) > 0) {
      stop("set with no palette key: ", paste(unmapped, collapse = ", "),
           "\n  Add it to palette_key.", call. = FALSE)
    }
    stripe_colors <- unname(
      MotrpacHumanPreSuspensionAnalysis::HUMAN_OME_COLORS[
        unname(palette_key[present_sets])]
    )
    if (anyNA(stripe_colors)) {
      stop("ome not in HUMAN_OME_COLORS: ",
           paste(unname(palette_key[present_sets])[is.na(stripe_colors)],
                 collapse = ", "), call. = FALSE)
    }

    tissue_fill <- unname(
      MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[tolower(tissue_name)]
    )

    ComplexUpset::upset(
      indicators,
      present_sets,
      encode_sets = FALSE,
      width_ratio = 0.2,
      height_ratio = 1,
      set_sizes = FALSE,
      base_annotations = list(
        "Number of \ndifferential genes" = intersection_size(
          counts = TRUE,
          fill = tissue_fill,
          color = "black",
          text = list(size = 4, color = "black")
        )
      ),
      themes = upset_modify_themes(
        list(intersections_matrix = theme(axis.title.x = element_blank()))
      ),
      sort_sets = FALSE,
      sort_intersections = "descending",
      wrap = TRUE,
      stripes = stripe_colors
    ) +
      ggtitle(toupper(tissue_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 16))
  }

  # ---- plot ----

  tissue_plots <- lapply(c("Adipose", "Blood", "Muscle"), tissue_upset)

  combined <- patchwork::wrap_plots(tissue_plots, ncol = 1)

  export_panel(function() print(combined), "ED2F")
}

# ---- ST2f — phospho ORA (PTMsigDB) ------------------------------------------

# The vendored CSV re-emitted under this repo's column schema, with
# statistic_column = -log10(p_value). The published result, not the
# recomputation ED2E logs beside it.
write_st2f <- function(ora) {
  export_table(
    as.data.frame(ora)[, table_spec("ST2f")$columns, drop = FALSE],
    "ST2f"
  )
}

run_panels(list(ED2A = ed2a, ED2B = ed2b, ED2C = ed2c, ED2D = ed2d,
                ED2E = ed2e, ED2F = ed2f))
