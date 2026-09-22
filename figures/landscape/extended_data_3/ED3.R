#!/usr/bin/env Rscript
# Extended Data 3 — Sex sensitivity and cellular deconvolution
#
# Panels:  ED3A  male vs female blood transcriptional response concordance
#          ED3B  PPARGC1A muscle transcript trajectory, split by sex
#          ED3C  CCN1 transcript trajectory in muscle and adipose, split by sex
#          ED3D  CIBERSORTx neutrophil estimates vs measured CBC neutrophils
#          ED3E  blood RNA-seq DA with and without cell-type covariates
# Tables:  ST2g   sex-specific DA estimate correlation   (written by ED3A)
#          ST2h   sex-specific differences enrichment    (written by ED3A)
#
# ---------------------------------------------------------------------------
# ED3D cannot be run from this repository, and its input cannot be released.
#
# It reads a CIBERSORTx deconvolution of MoTrPAC blood RNA-seq: 951 rows keyed
# by vial label and 22 cell-type proportions per sample. That is individual-level
# data: it is not distributed with this repository and may not be published.
# Running ED3D requires an approved MoTrPAC consortium data-access request and a
# local copy of the file at extended_data_3/sources/, or ED3D_CIBERSORTX_CSV
# pointing at one.
#
# ED3E is downstream of the same file. It reads a fit rather than the
# deconvolution itself, so it says which script to run rather than which file is
# missing — but analysis/02_celltype_da.R rebuilds its five covariates from that
# deconvolution and cannot be run here either.
#
# ED3A, ED3B and ED3C run without any of it.
# ---------------------------------------------------------------------------
#
# Needs consortium data access. ED3A reads the sex-stratified differential
# analysis, which neither data package carries: run
# `Rscript figures/landscape/analysis/01_sex_da.R` first. The fit is hours and is
# cached, so this is a one-time cost.
#
#   Rscript figures/landscape/ED3.R          every panel
#   Rscript figures/landscape/ED3.R ED3B     one panel

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "ED3_helpers.R"))

# ---- shared ----------------------------------------------------------------

# Nothing is shared here, and that is a judgement rather than an omission. The
# five panels read five different things: ED3A the sex-stratified DA directory
# plus load_differential_analysis() over two omes and every tissue; ED3B and
# ED3C load_qc() for different tissue sets; ED3D load_qc() for blood
# transcriptomics alone; ED3E load_differential_analysis() for blood
# transcriptomics alone. Two loads that differ in their arguments are two loads,
# so each stays inside its panel.

# ---- ED3A — male vs female blood response concordance ----------------------

# Exercise-induced logFC (post-exercise vs pre-exercise) is estimated separately
# in male and in female participants by a sex-stratified differential analysis,
# and the two estimates are plotted against one another for every blood
# transcript that was significant in the primary, sex-combined PreCAWG analysis.
# Rows are post-exercise timepoints, columns are exercise mode; each panel
# carries its Pearson correlation and OLS fit.
#
# Each timepoint row is labelled with its own ten transcripts, ranked by
# male-vs-female difference adjusted p-value at that timepoint. The ranking is
# made in resistance exercise (LABEL_SOURCE_GROUP) and the resulting labels are
# drawn in both exercise-mode columns, so a row names the same transcripts on
# both sides and the columns stay comparable.
#
# This panel also writes ST2g and ST2h. The table reports the same per-facet
# correlations over every tissue x ome the sex-stratified fit covers rather than
# the one plotted, and ST2h enriches the features that differ between the sexes;
# both come out of the one read of the fit below.

# Labels kept per timepoint row, and the exercise group whose male-vs-female
# difference p-value ranks them.
N_LABELS_PER_TIMEPOINT <- 10
LABEL_SOURCE_GROUP <- "ADUResist"

ed3a <- function() {
  skip <- skip_if_no_fit("ED3A")
  if (!is.null(skip)) return(skip)

  panel_init("ED3A")

  # ---- input ---------------------------------------------------------------
  #
  # The sex-stratified differential analysis is not in either data package.
  # analysis/01_sex_da.R fits it — model
  # `~ 0 + sex_group_timepoint + <covariates> + (1 | pid)` — and writes one
  # tab-delimited table per tissue x platform named
  # *_da_dream-sex_differences_*.txt into outputs/fits/sex_da/.
  # SEX_STRATIFIED_DA_DIR points there by default and takes an override: set it
  # to a directory of tables fitted elsewhere and the fit is bypassed. Either way
  # this panel only reads; it fits nothing.
  sex_da_dir <- Sys.getenv("SEX_STRATIFIED_DA_DIR", unset = "")

  # The same directory also collects DA runs from other models, so the file set is
  # pinned by the model name.
  #
  # Every tissue x ome is read, not just the one the panel plots. ST2g reports the
  # correlation for all of them, and reading the directory once is what keeps the
  # panel's own facet and the table's blood-transcriptomics row the same number.
  # The panel filters down to blood transcriptomics further below.
  sex_da_files <- list.files(
    sex_da_dir,
    pattern = "_da_dream-sex_differences_.*\\.txt$",
    full.names = TRUE
  )
  if (length(sex_da_files) == 0) {
    stop("no sex-differences DA tables under ", sex_da_dir,
         ". Run `Rscript figures/landscape/analysis/01_sex_da.R` to fit them, ",
         "or point SEX_STRATIFIED_DA_DIR at tables fitted elsewhere.",
         call. = FALSE)
  }
  if (!any(grepl("blood-rna", basename(sex_da_files)) &
           grepl("transcript-rna-seq", basename(sex_da_files)))) {
    stop("no blood transcript-rna-seq sex-differences DA table under ", sex_da_dir,
         " — that is the one this panel is drawn from.", call. = FALSE)
  }

  # ---- data ----------------------------------------------------------------

  read_sex_diff <- lapply(sex_da_files, function(f) {
    tissue <- dplyr::case_when(
      grepl("muscle", basename(f)) ~ "muscle",
      grepl("blood|plasma", basename(f)) ~ "blood",
      grepl("adipose", basename(f)) ~ "adipose"
    )
    read.csv(f, sep = "\t", check.names = FALSE) %>%
      dplyr::mutate(tissue = tissue)
  }) %>%
    dplyr::bind_rows()

  precawg_analysis <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    # Transcriptomics and metabolomics, every tissue: the omes and tissues the
    # sex-stratified fit covers, because pre_cawg_adj_p is what selects the
    # features every correlation in ST2g is taken over.
    selected_omes = c("transcript-rna-seq", "metab"),
    selected_tissues = "all",
    single_matrix = TRUE
  ) %>%
    dplyr::filter(contrast_type == "exercise_with_controls") %>%
    # No CI.L/CI.R. The released acute DA does not carry them: precovid-repro's
    # writer drops them deliberately, because variancePartition's topTable computes
    # them wrongly for a dream fit (df.total is a features x contrasts matrix and
    # `top` indexes it linearly, so every contrast gets the first contrast's df).
    #
    # Nothing is lost. This table is used only for its adj_p_value — which
    # transcripts were significant in the primary analysis — and the error bars
    # this panel draws come from the sex-stratified table, where
    # analysis/01_sex_da.R rebuilds the interval correctly per contrast.
    dplyr::select(tissue, assay, feature_id, logFC, adj_p_value,
                  randomGroupCode, Timepoint)

  timepoint_levels <- c("pre_exercise", "during_20_min", "during_40_min",
                        "post_10_min", "post_15_30_45_min", "post_3.5_4_hr",
                        "post_24_hr")

  # The contrast string is "<sex>.<group>.<timepoint> - <reference>", wrapped in
  # the model term name and in parentheses. Both the within-sex contrasts and the
  # male-vs-female difference contrasts are named this way, so both are parsed
  # here; a caller that does not need sex_comparison drops it in its own select().
  parse_sex_contrast <- function(df) {
    df %>%
      dplyr::mutate(
        contrast_short = gsub("sex_group_timepoint|[()|]", "", contrast),
        contrast_left = stringr::str_split_fixed(contrast_short, " - ", 2)[, 1],
        sex_comparison = stringr::str_split_fixed(contrast_left, "\\.", 3)[, 1],
        randomGroupCode = stringr::str_split_fixed(contrast_left, "\\.", 3)[, 2],
        Timepoint = sub("^[^.]+\\.[^.]+\\.", "", contrast_left)
      ) %>%
      dplyr::select(-contrast_left) %>%
      dplyr::mutate(Timepoint = factor(Timepoint, levels = timepoint_levels)) %>%
      dplyr::rename(platform = assay) %>%
      dplyr::mutate(assay = ifelse(grepl("metab", platform), "metab", platform))
  }

  # One row per sex x group x timepoint: the within-sex change from pre-exercise.
  sex_single_da_annotated <- read_sex_diff %>%
    dplyr::filter(!(grepl("Male", contrast) & grepl("Female", contrast))) %>%
    parse_sex_contrast() %>%
    dplyr::select(tissue, assay, platform, feature_id, logFC, CI.L, CI.R,
                  adj_p_value, sex_comparison, randomGroupCode, Timepoint) %>%
    dplyr::left_join(
      precawg_analysis %>%
        dplyr::select(tissue, assay, feature_id, Timepoint, randomGroupCode,
                      pre_cawg_adj_p = adj_p_value),
      # assay is part of the key. Without it a feature_id measured on more than one
      # ome joins across them, which cannot happen while the primary analysis is
      # loaded for one ome and silently does once it is loaded for two.
      by = c("tissue", "assay", "Timepoint", "randomGroupCode", "feature_id")
    )

  # One row per feature x group x timepoint, female and male estimates side by side.
  scatter_wide <- sex_single_da_annotated %>%
    dplyr::select(tissue, assay, platform, feature_id, randomGroupCode, Timepoint,
                  sex_comparison, logFC, CI.L, CI.R, adj_p_value, pre_cawg_adj_p) %>%
    tidyr::pivot_wider(
      names_from = sex_comparison,
      values_from = c(logFC, CI.L, CI.R, adj_p_value)
    ) %>%
    dplyr::filter(!is.na(logFC_Male), !is.na(logFC_Female)) %>%
    dplyr::mutate(
      Ome = dplyr::case_when(
        grepl("metab", assay)              ~ "Metabolomics",
        assay == "prot-ph"                 ~ "Phosphoproteomics",
        assay %in% c("prot-pr", "prot-ol") ~ "Proteomics",
        grepl("rna", assay)                ~ "Transcriptomics",
        TRUE ~ assay
      ),
      Ome = factor(Ome, levels = c("Transcriptomics", "Proteomics",
                                   "Phosphoproteomics", "Metabolomics")),
      Group = dplyr::recode(randomGroupCode,
                            "ADUEndur" = "EE", "ADUResist" = "RE",
                            "ADUControl" = "CON"),
      Group = factor(Group, levels = c("EE", "RE", "CON"))
    )

  # Correlations are over the features the primary analysis called significant,
  # computed per ome x tissue x group x timepoint.
  scatter_sig_features <- scatter_wide %>%
    dplyr::filter(!is.na(pre_cawg_adj_p) & pre_cawg_adj_p < 0.05, Group != "CON")

  scatter_cor <- scatter_sig_features %>%
    dplyr::group_by(Ome, tissue, Group, Timepoint) %>%
    dplyr::summarise(
      n_sig = dplyr::n(),
      Pearsons_R = round(stats::cor(logFC_Female, logFC_Male, method = "pearson",
                                    use = "complete.obs"), 3),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(c(Ome, tissue, Group, Timepoint), as.character))

  write_st2g(scatter_sig_features)

  tp_labels <- c(
    "post_10_min"       = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr"     = "P3.5/4H"
  )
  post_timepoints <- names(tp_labels)

  prettify_timepoint <- function(v) {
    factor(dplyr::recode(as.character(v), !!!tp_labels), levels = unname(tp_labels))
  }

  prettify_facets <- function(df) {
    df %>%
      dplyr::mutate(
        Timepoint = prettify_timepoint(Timepoint),
        Group = factor(Group, levels = c("EE", "RE"))
      )
  }

  # The panel is blood transcriptomics. scatter_wide now carries every tissue x
  # ome ST2g reports, so the assay is pinned here as well as the tissue.
  blood_post_raw <- scatter_wide %>%
    dplyr::filter(tissue == "blood",
                  assay == "transcript-rna-seq",
                  Group %in% c("EE", "RE"),
                  Timepoint %in% post_timepoints,
                  pre_cawg_adj_p < 0.05)

  blood_post <- prettify_facets(blood_post_raw)

  blood_post_cor <- scatter_cor %>%
    dplyr::filter(tissue == "blood", Group %in% c("EE", "RE"),
                  Timepoint %in% post_timepoints) %>%
    prettify_facets()

  # ---- label selection -----------------------------------------------------

  # The male-vs-female difference contrast, joined to the primary analysis, ranks
  # the transcripts to label.
  sex_direct_da_annotated <- read_sex_diff %>%
    dplyr::filter(grepl("Male", contrast) & grepl("Female", contrast)) %>%
    parse_sex_contrast() %>%
    dplyr::select(tissue, assay, platform, feature_id, logFC, CI.L, CI.R,
                  adj_p_value, randomGroupCode, Timepoint)

  precawg_and_m_vs_f_diff <- precawg_analysis %>%
    dplyr::right_join(sex_direct_da_annotated,
                      by = c("tissue", "assay", "feature_id", "randomGroupCode",
                             "Timepoint"),
                      suffix = c("_precawg", "_sex_diff"))

  write_st2h(precawg_and_m_vs_f_diff)

  # Candidates are restricted to what the label-source column actually draws, so a
  # row gets N_LABELS_PER_TIMEPOINT labels rather than N minus however many of its
  # top-ranked transcripts were filtered out of the scatter.
  # Timepoint as character on both sides of the semi_join: it is a factor here and
  # a character in precawg_and_m_vs_f_diff, whose join against the primary analysis
  # already coerced it.
  plotted_in_source_group <- blood_post_raw %>%
    dplyr::filter(randomGroupCode == LABEL_SOURCE_GROUP) %>%
    dplyr::distinct(feature_id, Timepoint) %>%
    dplyr::mutate(Timepoint = as.character(Timepoint))

  highlight_by_timepoint <- precawg_and_m_vs_f_diff %>%
    dplyr::filter(tissue == "blood",
                  assay == "transcript-rna-seq",
                  randomGroupCode == LABEL_SOURCE_GROUP,
                  Timepoint %in% post_timepoints) %>%
    dplyr::mutate(Timepoint = as.character(Timepoint)) %>%
    dplyr::semi_join(plotted_in_source_group, by = c("feature_id", "Timepoint")) %>%
    dplyr::group_by(Timepoint) %>%
    dplyr::slice_min(adj_p_value_sex_diff, n = N_LABELS_PER_TIMEPOINT,
                     with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(Timepoint = prettify_timepoint(Timepoint), feature_id)

  # ---- plot ----------------------------------------------------------------

  p <- make_sex_scatter("blood", "Transcriptomics", blood_post, blood_post_cor,
                        highlight_df = highlight_by_timepoint)
  if (is.null(p)) {
    stop("no blood transcriptomics rows survived the ED3A filters", call. = FALSE)
  }

  p <- p + labs(title = "Male vs. female transcriptional response to exercise in blood")

  export_panel(p, "ED3A")
}

#' Female-vs-male logFC scatter for one tissue x ome, faceted timepoint x group.
#'
#' @param tiss Tissue to draw.
#' @param ome_type Ome label to draw.
#' @param data Wide female/male logFC table, already filtered to the facets.
#' @param cor_df Per-facet correlation table, with n_sig and Pearsons_R.
#' The error bars and points are rasterised whatever their count —
#' `threshold = 0` — rather than only above maybe_rasterise()'s density cutoff.
#' That is what the legacy does (`ggrastr::rasterise(..., dpi = 600)` on the same
#' three layers), and this panel draws six facets of per-feature error bars in
#' both x and y, so the vector object count is several times the point count the
#' cutoff is measured against. Axes, strips, labels and the fit line stay vector.
#'
#' @param highlight_df Timepoint x feature_id table of what to label. Its
#'   Timepoint must already carry the prettified labels `data` is faceted on.
#'   Every feature listed for a timepoint is labelled in both group columns of
#'   that row, so the label set is not trimmed further here.
make_sex_scatter <- function(tiss, ome_type, data, cor_df, highlight_df = NULL) {
  df <- data %>% dplyr::filter(tissue == tiss, Ome == ome_type)
  if (nrow(df) == 0) return(NULL)

  cor_sub <- cor_df %>% dplyr::filter(tissue == tiss, Ome == ome_type)

  lim <- max(abs(c(df$logFC_Female, df$logFC_Male)), na.rm = TRUE) * 1.05

  lm_labels <- df %>%
    dplyr::group_by(Timepoint, Group) %>%
    dplyr::summarize(
      slope_val = tryCatch(stats::coef(stats::lm(logFC_Male ~ logFC_Female))[2],
                           error = function(e) NA_real_),
      intercept_val = tryCatch(stats::coef(stats::lm(logFC_Male ~ logFC_Female))[1],
                               error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      eq_label = dplyr::case_when(
        is.na(slope_val) ~ "",
        intercept_val >= 0 ~ paste0("y=", sprintf("%.2f", slope_val), "x+",
                                    sprintf("%.2f", intercept_val)),
        TRUE ~ paste0("y=", sprintf("%.2f", slope_val), "x",
                      sprintf("%.2f", intercept_val))
      )
    )

  label_df <- NULL
  if (!is.null(highlight_df) && nrow(highlight_df) > 0) {
    gene_map <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
      dplyr::select(feature_id, gene_symbol, refmet_name) %>%
      dplyr::distinct()

    # distinct() after the gene join: a feature_id mapping to more than one
    # gene_symbol row would otherwise draw its label twice in the same facet.
    label_df <- df %>%
      dplyr::semi_join(highlight_df, by = c("feature_id", "Timepoint")) %>%
      dplyr::left_join(gene_map, by = "feature_id") %>%
      dplyr::mutate(label = dplyr::coalesce(gene_symbol, refmet_name, feature_id)) %>%
      dplyr::distinct(feature_id, Timepoint, Group, .keep_all = TRUE)
  }

  # Timepoint rows and exercise-group columns take the canonical palettes.
  # Timepoint labels may be raw ("post_3.5_4_hr") or prettified ("P3.5/4H");
  # both resolve to the same colour.
  tp_pal <- MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  tp_label_to_raw <- c(
    "Pre" = "pre_exercise",
    "D20" = "during_20_min", "D20M" = "during_20_min",
    "D40" = "during_40_min", "D40M" = "during_40_min",
    "P10" = "post_10_min", "P10M" = "post_10_min",
    "P1545" = "post_15_30_45_min", "P15-45M" = "post_15_30_45_min",
    "P35" = "post_3.5_4_hr", "P3.5/4H" = "post_3.5_4_hr",
    "P24" = "post_24_hr", "P24H" = "post_24_hr"
  )
  group_pal <- .exercise_group_palette()

  resolve_tp <- function(v) {
    v <- as.character(v)
    if (v %in% names(tp_pal)) return(unname(tp_pal[v]))
    if (v %in% names(tp_label_to_raw)) return(unname(tp_pal[tp_label_to_raw[v]]))
    "grey85"
  }
  resolve_grp <- function(v) {
    v <- as.character(v)
    if (v %in% names(group_pal)) unname(group_pal[v]) else "grey85"
  }

  row_vals <- levels(droplevels(as.factor(df$Timepoint)))
  col_vals <- levels(droplevels(as.factor(df$Group)))
  row_cols <- vapply(row_vals, resolve_tp, character(1))
  col_cols <- vapply(col_vals, resolve_grp, character(1))

  strip_spec <- .themed_strip(col_cols, row_cols)

  p <- ggplot(df, aes(x = logFC_Female, y = logFC_Male)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "grey30", linewidth = 0.5) +
    maybe_rasterise(
      geom_errorbar(
        aes(ymin = CI.L_Male, ymax = CI.R_Male),
        width = 0, alpha = 0.3, linewidth = 0.35, color = "grey40"
      ),
      n = nrow(df), panel = "ED3A", threshold = 0
    ) +
    maybe_rasterise(
      geom_errorbarh(
        aes(xmin = CI.L_Female, xmax = CI.R_Female),
        height = 0, alpha = 0.3, linewidth = 0.35, color = "grey40"
      ),
      n = nrow(df), panel = "ED3A", threshold = 0
    ) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.7) +
    maybe_rasterise(
      geom_point(alpha = 0.8, size = 1.5, color = "grey30"),
      n = nrow(df), panel = "ED3A", threshold = 0
    ) +
    geom_text(
      data = cor_sub,
      aes(label = paste0("r=", Pearsons_R, " (n=", n_sig, ")")),
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.2,
      size = 5, fontface = "italic",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = lm_labels,
      aes(label = eq_label),
      x = -Inf, y = Inf,
      hjust = -0.05, vjust = 2.4,
      size = 5, fontface = "italic",
      inherit.aes = FALSE
    ) +
    coord_fixed(ratio = 1, xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggh4x::facet_grid2(Timepoint ~ Group, strip = strip_spec) +
    labs(
      title = paste0(tiss, " - ", ome_type),
      x     = "logFC (Female participants)",
      y     = "logFC (Male participants)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 9, face = "bold"),
      axis.title       = element_text(size = 15),
      axis.text        = element_text(size = 12),
      legend.position  = "bottom"
    )

  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(x = logFC_Female, y = logFC_Male, label = label),
      size = 3.5,
      fontface = "bold",
      bg.color = "white",
      bg.r = 0.15,
      max.overlaps = 30,
      segment.size = 0.3,
      segment.color = "grey40",
      inherit.aes = FALSE
    )
  }

  p
}

# ---- ST2g — the correlation summary ----------------------------------------

# The same per-facet numbers the panel annotates its scatter with, over every
# tissue x ome the sex-stratified fit covers rather than the one plotted. slope
# and intercept are the OLS fit of male on female logFC, which is the line the
# panel draws; they are computed here rather than in the plotting helper so the
# table and the annotation cannot come from two different fits.
write_st2g <- function(scatter_sig_features) {
  #' n, Pearson r and the OLS fit of male on female logFC for one set of features.
  #'
  #' One definition used twice, so the per-cell rows and the pooled row cannot end
  #' up computed differently. An OLS fit needs two complete pairs; fewer gives NA
  #' rather than an error.
  cor_summary <- function(df) {
    complete <- !is.na(df$logFC_Female) & !is.na(df$logFC_Male)
    fit <- if (sum(complete) >= 2) {
      stats::coef(stats::lm(df$logFC_Male[complete] ~ df$logFC_Female[complete]))
    } else {
      c(NA_real_, NA_real_)
    }
    data.frame(
      n_sig = nrow(df),
      Pearsons_R = round(stats::cor(df$logFC_Female, df$logFC_Male,
                                    method = "pearson", use = "complete.obs"), 3),
      slope = unname(fit[2]),
      intercept = unname(fit[1])
    )
  }

  # Ordered by the ome factor's own levels rather than alphabetically:
  # Transcriptomics before Metabolomics, which is the order the omes are given in
  # everywhere else and the order the reference sheet uses. Timepoint likewise —
  # both are factors here, and sorting them as character would put post_10_min
  # before pre_exercise.
  st2g_cells <- scatter_sig_features %>%
    dplyr::arrange(.data$Ome, .data$tissue, .data$Group, .data$Timepoint) %>%
    dplyr::group_by(.data$Ome, .data$tissue, .data$Group, .data$Timepoint) %>%
    dplyr::group_modify(~ cor_summary(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(c("Ome", "tissue", "Group", "Timepoint"),
                                as.character))

  # One pooled row, last: every significant feature across every tissue, ome,
  # group and timepoint at once. It answers a different question from the cells
  # above — whether male and female responses agree overall — so it is a row of
  # the table rather than a summary of the rows above it, and it is not their
  # average.
  st2g_overall <- cbind(
    data.frame(Ome = "All", tissue = "All", Group = "All", Timepoint = "All"),
    cor_summary(scatter_sig_features)
  )

  st2g <- as.data.frame(dplyr::bind_rows(st2g_cells, st2g_overall))

  export_table(st2g[, table_spec("ST2g")$columns, drop = FALSE], "ST2g")
}

# ---- ST2h — enrichment of the sex-different features ------------------------

# The one cell where the male-vs-female difference is large enough to enrich:
# blood transcriptomics, resistance arm, 3.5/4 hours. A feature qualifies when it
# was significant in the primary analysis and also differs between the sexes.
#
# The two cuts are not the same number, and the difference is deliberate. The
# primary analysis is well powered here - 3,773 of 20,526 features clear FDR 0.05
# - so 0.05 is what selects a responding feature. The sex-stratified contrast is
# the sensitivity analysis, fitted on half the samples per cell, and only 21
# features clear FDR 0.05 in it against 111 at 0.10. Holding both at 0.05 leaves
# 18 input genes, which is too few for an over-representation test to say
# anything stable; 0.10 on the sex-difference cut gives 95.
#
# Relaxing the primary cut instead would change nothing - the same 18 genes -
# because the sex-difference contrast is the only binding constraint.
#
# Note this departs from the table legend, which describes the inputs as
# differing between the sexes at FDR < 0.05. The legend needs updating with the
# table.
ST2H_PRECAWG_FDR <- 0.05
ST2H_SEX_DIFF_FDR <- 0.10

write_st2h <- function(precawg_and_m_vs_f_diff) {
  # The background is every blood transcript carrying a gene symbol, not only the
  # tested ones, which is what makes the over-representation a statement about the
  # measured transcriptome.
  st2h_annotated <- precawg_and_m_vs_f_diff %>%
    dplyr::left_join(
      MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
        dplyr::filter(assay == "transcript-rna-seq") %>%
        dplyr::select(feature_id, assay, gene_symbol),
      by = c("feature_id", "assay")
    )

  st2h_blood_rna <- st2h_annotated %>%
    dplyr::filter(tissue == "blood", assay == "transcript-rna-seq")

  st2h_input <- st2h_blood_rna %>%
    dplyr::filter(Timepoint == "post_3.5_4_hr",
                  randomGroupCode == "ADUResist",
                  adj_p_value_precawg < ST2H_PRECAWG_FDR,
                  adj_p_value_sex_diff < ST2H_SEX_DIFF_FDR)

  if (nrow(st2h_input) == 0) {
    stop("ST2h: no blood transcriptomics feature is significant in both the ",
         "primary analysis and the male-vs-female contrast at RE / 3.5-4 h, so ",
         "there is nothing to enrich.", call. = FALSE)
  }

  st2h <- MotrpacHumanPreSuspensionAnalysis::run_ORA(
    input = as.character(unique(st2h_input$gene_symbol)),
    background = as.character(unique(st2h_blood_rna$gene_symbol))
  ) %>%
    as.data.frame()

  message(sprintf(
    "[ST2h] %d input gene(s) at precawg FDR < %.2f and sex-difference FDR < %.2f, %d background, %d set(s) returned",
    dplyr::n_distinct(st2h_input$gene_symbol), ST2H_PRECAWG_FDR, ST2H_SEX_DIFF_FDR,
    dplyr::n_distinct(st2h_blood_rna$gene_symbol), nrow(st2h)))

  export_table(st2h[, table_spec("ST2h")$columns, drop = FALSE], "ST2h")
}

# ---- ED3B — PPARGC1A muscle transcript trajectory, split by sex ------------

# Mean normalized abundance of every muscle feature mapping to PPARGC1A, with
# its 95% t confidence interval, at each acute timepoint. Three series per
# panel: female, male, and the pooled overall mean. Rows are the muscle tissue x
# assay combinations carrying the gene, columns are the exercise groups.
# Concordant female and male trajectories show the reported response is not an
# artefact of the female-majority cohort.
#
# The plot is one entry of the single-feature catalog in single_feature_plots.R.
# Unlike the rest of the catalog it is drawn by sex_differences_single_feature()
# from helpers/ED3.R, not by plot_single_feature().

ed3b <- function() {
  panel_init("ED3B")
  export_panel(single_feature_plot("ED3B"), "ED3B")
}

# ---- ED3C — CCN1 transcript trajectory in muscle and adipose ---------------

# Mean normalized abundance of the CCN1 transcript, with its 95% t confidence
# interval, at each acute timepoint in muscle and adipose. Three series per
# panel: female, male, and the pooled overall mean; columns are the exercise
# groups. The feature is named by its Ensembl id rather than by the gene symbol
# so that only the RNA-seq transcript is drawn - the symbol also matches ATAC
# peaks and methylation sites, which would add rows the panel is not about.

ed3c <- function() {
  panel_init("ED3C")
  export_panel(single_feature_plot("ED3C"), "ED3C")
}

# ---- ED3D — CIBERSORTx neutrophils vs measured CBC neutrophils -------------

# Blood bulk RNA-seq is deconvolved with CIBERSORTx against the LM22 leukocyte
# signature; the 22 subsets are collapsed into coarse populations and rescaled to
# percentages per sample. Pre-exercise neutrophil percentages are then joined by
# participant to the screening complete blood count and plotted against each
# other, with an OLS fit and the Pearson correlation. Two CBC rows carrying
# implausibly low neutrophil values (3.6 and 4.1) are excluded first.
#
# One point per participant: the CBC table also carries a follow-up visit, and
# taking both would plot 47 participants twice — see plot_cbc_vs_cibersortx().
#
# The deconvolution it reads is individual-level data, is not in this repository
# and cannot be released; see the note at the top of this file.

ed3d <- function() {
  skip <- skip_if_gated(
    "ED3D",
    path = "extended_data_3/sources/CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv",
    env_var = "ED3D_CIBERSORTX_CSV",
    input = paste("a CIBERSORTx deconvolution of MoTrPAC blood RNA-seq carrying",
                  "951 rows keyed by vial label and 22 cell-type proportions per",
                  "sample"),
    access = "an approved MoTrPAC consortium data-access request"
  )
  if (!is.null(skip)) return(skip)

  spec <- panel_init("ED3D")

  # ---- inputs ----
  #
  # The CIBERSORTx deconvolution result: one row per sample vial label in a
  # Mixture column, one column per LM22 leukocyte subset, plus the P-value,
  # Correlation and RMSE run-quality columns. It is the download from a run of the
  # CIBERSORTx web service against the LM22 signature matrix; the service has no R
  # implementation, so nothing here regenerates it. The run parameters it was
  # produced under are recorded in docs/external_dependencies.md — they decide the
  # numbers, and two runs of the same service over the same mixtures do not agree
  # if they differ.
  #
  # analysis/02_celltype_da.R reads this same file for the covariates it fits
  # ED3E's table with. One copy, so the fractions in the model and the fractions
  # in this figure cannot drift apart.

  cibersortx_csv <- panel_source(spec$figure_dir,
                                 "CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv",
                                 env_var = "ED3D_CIBERSORTX_CSV")

  # check.names = FALSE keeps the CIBERSORTx column labels as written ("T cells
  # CD8"); squish_proportions() matches them through make.names().
  cibersort_results <- read.csv(cibersortx_csv, check.names = FALSE)

  # ---- data ----

  # P-value/Correlation/RMSE describe the fit of each mixture, not its
  # composition, and must not be summed into a cell population. Both spellings are
  # accepted so a CSV read with check.names = TRUE also works.
  cibersort_merge_df <- cibersort_results |>
    select(-any_of(c("P-value", "P.value", "Correlation", "RMSE"))) |>
    squish_proportions() |>
    rename(vialLabel = Mixture) |>
    mutate(vialLabel = as.character(vialLabel))

  blood_qc <- load_qc(selected_tissues = "blood",
                      selected_omes = "transcript-rna-seq")

  blood_metadata <- blood_qc$blood$`transcript-rna-seq`$sample_metadata |>
    mutate(vialLabel = as.character(vialLabel)) |>
    select(vialLabel, pid, Timepoint, randomGroupCode) |>
    inner_join(cibersort_merge_df, by = "vialLabel") |>
    select(-any_of("Eosinophils"))

  cbc_raw <- MotrpacHumanPreSuspensionData::cln_raw_local_lab_results$data

  # ---- plot ----

  p <- plot_cbc_vs_cibersortx(blood_metadata, cbc_raw, filter_lt5 = TRUE)

  export_panel(p, "ED3D")
}

# ---- ED3E — blood RNA-seq DA with and without cell-type covariates ---------

# Per-transcript exercise effects from the published blood RNA-seq differential
# analysis are plotted against effects from a re-fit that additionally carries
# the z-scored coarse cell-type fractions (T cells, B cells, NK cells,
# monocytes, neutrophils) as numeric covariates. One point per transcript,
# faceted by exercise group and timepoint, against the line of identity with the
# Pearson correlation annotated per facet. The metric plotted is logFC.
#
# The input this panel reads is fitted by analysis/02_celltype_da.R, not carried
# by either data package — and that fit rebuilds its covariates from the gated
# deconvolution, so it cannot be run here either.

ed3e <- function() {
  # skip_if_no_fit() reads the manifest's override_env_var, so a table fitted
  # elsewhere and named in ED3E_DA_WITH_CELL_TYPES satisfies the input on its own.
  skip <- skip_if_no_fit("ED3E")
  if (!is.null(skip)) return(skip)

  panel_init("ED3E")

  # ---- inputs ----
  #
  # The blood RNA-seq differential analysis re-fit with cell-type covariates: the
  # five z-scored CIBERSORTx cell fractions added to the covariate table and the
  # standard acute dream model re-run on the blood RSEM gene counts. Tab-delimited
  # and ~330 MB.
  #
  # analysis/02_celltype_da.R fits it and writes it into outputs/fits/celltype_da/
  # under the freeze naming convention. CELLTYPE_DA_DIR points there by default.
  # ED3E_DA_WITH_CELL_TYPES overrides the lookup with a single table fitted
  # elsewhere — the frozen v1.3 artifact, say — and bypasses the fit. Either way
  # this panel only reads; it fits nothing.

  da_with_cell_types_path <- Sys.getenv("ED3E_DA_WITH_CELL_TYPES", unset = "")

  if (nzchar(da_with_cell_types_path)) {
    if (!file.exists(da_with_cell_types_path)) {
      stop("ED3E_DA_WITH_CELL_TYPES names a file that is not there: ",
           da_with_cell_types_path, call. = FALSE)
    }
  } else {
    celltype_da_dir <- Sys.getenv("CELLTYPE_DA_DIR", unset = "")

    # The same directory also collects the covariate table the fit was built from,
    # so the file set is pinned by the model name and by the tissue x assay this
    # panel plots — not by anything that merely ends in .txt.
    hits <- list.files(celltype_da_dir,
                       pattern = "_da_dream-acute-cell_types_.*\\.txt$",
                       full.names = TRUE)
    hits <- hits[grepl("blood-rna", basename(hits)) &
                   grepl("transcript-rna-seq", basename(hits))]

    if (length(hits) == 0) {
      stop("no blood transcript-rna-seq cell-type DA table under ", celltype_da_dir,
           ". Run `Rscript figures/landscape/analysis/02_celltype_da.R` to fit it, ",
           "or point ED3E_DA_WITH_CELL_TYPES at a table fitted elsewhere.",
           call. = FALSE)
    }
    if (length(hits) > 1) {
      # Two tables differing only in version is two different fits, and picking
      # silently between them is how a panel ends up labelled with the wrong one.
      stop("several blood transcript-rna-seq cell-type DA tables under ",
           celltype_da_dir, ": ", paste(basename(hits), collapse = ", "),
           " — remove one, or name the one you want in ED3E_DA_WITH_CELL_TYPES",
           call. = FALSE)
    }
    da_with_cell_types_path <- hits[[1]]
  }

  message("[ED3E] cell-type DA: ", basename(da_with_cell_types_path))

  # fread reads a four-column subset of a ~330 MB tab-delimited table; read.csv
  # cannot do a partial read and would parse all fourteen columns to reach four.
  finished_da_with_cell_types <- data.table::fread(
    da_with_cell_types_path,
    sep = "\t",
    select = c("feature_id", "z.std", "logFC", "contrast"),
    data.table = FALSE,
    check.names = FALSE
  )

  # ---- data ----

  existing_motrpac_da <- load_differential_analysis(
    selected_omes = "transcript-rna-seq",
    selected_tissues = "blood",
    single_matrix = TRUE
  ) |>
    select(feature_id, z.std, logFC, contrast)

  combined_da <- finished_da_with_cell_types |>
    select(feature_id, z.std, logFC, contrast) |>
    left_join(
      existing_motrpac_da,
      by = c("feature_id", "contrast"),
      suffix = c(".cell_type", ".no_cell_type")
    ) |>
    left_join(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER,
              by = "contrast") |>
    filter(contrast_type == "exercise_with_controls")

  # ---- plot ----

  p <- plot_covariate_corr(combined_da, "logFC")

  export_panel(p, "ED3E")
}

# Facet strips are labelled with the short timepoint and group names and filled
# with the package colours for the same levels; only the keys are renamed, the
# colours themselves are the package constants.
# The short labels every other acute panel carries — ED3C draws them from
# helpers/ED3.R, ED3A and ED2A from their own copies of the same scheme — listed
# in sampling order, which is also the facet column order.
ED3E_TIMEPOINT_SHORT <- c(
  "pre_exercise"      = "Pre",
  "during_20_min"     = "D20M",
  "during_40_min"     = "D40M",
  "post_10_min"       = "P10M",
  "post_15_30_45_min" = "P15-45M",
  "post_3.5_4_hr"     = "P3.5/4H",
  "post_24_hr"        = "P24H"
)

ED3E_GROUP_SHORT <- c(
  "ADUEndur"   = "EE",
  "ADUResist"  = "RE",
  "ADUControl" = "CON"
)

# Rename through a lookup, leaving anything the map does not mention alone.
.abbreviate <- function(x, map) {
  hit <- match(x, names(map))
  ifelse(is.na(hit), x, unname(map)[hit])
}

# White strip text on dark fills, black on light fills.
.contrast_text <- function(fills) {
  rgb_mat <- grDevices::col2rgb(fills) / 255
  luminance <- 0.299 * rgb_mat[1, ] + 0.587 * rgb_mat[2, ] + 0.114 * rgb_mat[3, ]
  ifelse(luminance < 0.55, "white", "black")
}

#' Scatter one DA metric with cell-type covariates against the same metric
#' without them, faceted group x timepoint.
plot_covariate_corr <- function(da, metric) {
  timepoint_short <- ED3E_TIMEPOINT_SHORT
  group_short <- ED3E_GROUP_SHORT

  group_colors <- MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS
  names(group_colors) <- .abbreviate(names(group_colors), group_short)

  timepoint_colors <- MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS
  names(timepoint_colors) <- .abbreviate(names(timepoint_colors), timepoint_short)

  x_col <- paste0(metric, ".no_cell_type")
  y_col <- paste0(metric, ".cell_type")

  plot_df <- da |>
    filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]),
           randomGroupCode %in% c("ADUEndur", "ADUResist")) |>
    mutate(
      group = .abbreviate(as.character(randomGroupCode), group_short),
      timepoint = .abbreviate(as.character(Timepoint), timepoint_short)
    ) |>
    mutate(
      group = factor(group, levels = intersect(c("EE", "RE", "CON"), group)),
      timepoint = factor(
        timepoint,
        levels = unname(timepoint_short[timepoint_short %in% timepoint])
      )
    )

  corr_labels <- plot_df |>
    group_by(group, timepoint) |>
    summarise(
      r = stats::cor(.data[[x_col]], .data[[y_col]], use = "complete.obs"),
      .groups = "drop"
    ) |>
    mutate(label = paste0("R = ", round(r, 2)))

  col_fills <- unname(timepoint_colors[levels(plot_df$timepoint)])
  row_fills <- unname(group_colors[levels(plot_df$group)])

  ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    maybe_rasterise(geom_point(size = 0.4, alpha = 0.25),
                    n = nrow(plot_df), panel = "ED3E") +
    geom_abline(slope = 1, intercept = 0, color = "firebrick",
                linetype = "dashed", linewidth = 0.6) +
    geom_text(data = corr_labels,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.15, vjust = 1.5,
              size = 3.5, inherit.aes = FALSE) +
    ggh4x::facet_grid2(
      group ~ timepoint,
      strip = ggh4x::strip_themed(
        background_x = lapply(col_fills, function(f) {
          element_rect(fill = f, color = "grey20")
        }),
        background_y = lapply(row_fills, function(f) {
          element_rect(fill = f, color = "grey20")
        }),
        text_x = lapply(.contrast_text(col_fills), function(cl) {
          element_text(color = cl, size = 11, face = "bold")
        }),
        text_y = lapply(.contrast_text(row_fills), function(cl) {
          element_text(color = cl, size = 11, face = "bold")
        })
      )
    ) +
    labs(
      x = paste0(metric, " (no cell type covariate)"),
      y = paste0(metric, " (with cell type covariate)"),
      title = paste0(
        "Effect of adding cell type covariates on blood RNA-seq DA (", metric, ")"
      )
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      title = element_text(size = 13)
    )
}

run_panels(list(ED3A = ed3a, ED3B = ed3b, ED3C = ed3c, ED3D = ed3d, ED3E = ed3e))
