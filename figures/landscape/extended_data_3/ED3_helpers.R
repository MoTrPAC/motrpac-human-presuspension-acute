# helpers/ED3.R — shared pieces of Extended Data 3.
#
# Two groups of helpers, both read by more than one place:
#
#   sex sensitivity        `.text_on` is used by the coloured facet strips of
#                          both the single-feature panels (ED3B, ED3C) and the
#                          concordance scatter (ED3A); `.sex_feature_strip_spec`
#                          and `sex_differences_single_feature` are shared by
#                          ED3B and ED3C, and by single_feature_plots.R.
#   deconvolution          `squish_proportions` is read by ED3D and by
#                          analysis/02_celltype_da.R, which is the point: the
#                          covariates the fit uses have to be the same numbers
#                          ED3D plots, and two copies of the LM22 ->
#                          coarse-population map is exactly how they would stop
#                          being.
#
# Nothing here reads or writes a file. A panel gets a ggplot back and hands it
# to export_panel().

# ============================================================================
# Sex sensitivity — ED3A, ED3B, ED3C
# ============================================================================

# ---- facet strips ----------------------------------------------------------

#' Black or white text, whichever reads on a given background colour.
#'
#' @param hex A single colour, as hex or as an R colour name.
.text_on <- function(hex) {
  rgb <- grDevices::col2rgb(hex) / 255
  lum <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  if (lum > 0.55) "black" else "white"
}

#' A ggh4x strip spec with per-strip background fills and readable text.
#'
#' Each strip takes its own background colour, and its label is set in black or
#' white, whichever reads against that colour.
#'
#' @param col_cols Fill colour per column strip, in strip order.
#' @param row_cols Fill colour per row strip, in strip order.
#' @param text_size Strip text size. NULL leaves it at the theme default.
.themed_strip <- function(col_cols, row_cols, text_size = NULL) {
  strip_rect <- function(fill) {
    ggplot2::element_rect(fill = fill, colour = "grey20")
  }
  strip_text <- function(fill) {
    ggplot2::element_text(colour = .text_on(fill), face = "bold",
                          size = text_size)
  }

  ggh4x::strip_themed(
    background_x = lapply(col_cols, strip_rect),
    background_y = lapply(row_cols, strip_rect),
    text_x = lapply(col_cols, strip_text),
    text_y = lapply(row_cols, strip_text)
  )
}

#' The canonical exercise-group palette, keyed by the short group labels.
#'
#' HUMAN_EXERCISE_GROUP_COLORS keys on randomGroupCode; every ED3 panel labels
#' the groups EE / RE / CON. The colours themselves are the package constants.
.exercise_group_palette <- function() {
  pal <- MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS
  names(pal) <- dplyr::recode(names(pal),
                              "ADUEndur" = "EE", "ADUResist" = "RE",
                              "ADUControl" = "CON")
  pal
}

#' Coloured facet strips for the single-feature trajectory panels.
#'
#' Exercise-group columns take the canonical group palette and tissue rows the
#' canonical tissue palette. When every panel shares one ome, that ome is
#' dropped from the row strips - it is the same on every row - and returned in
#' the title instead, and the rows key on tissue alone.
#'
#' @param plot_data The data frame the panel is drawn from.
#' @param label_map Named vector mapping randomGroupCode to the strip label.
#' @param sc Text scale factor.
#'
#' @return A list with the ggh4x strip spec, the facet formula, and a function
#'   that builds the plot title from a feature label.
.sex_feature_strip_spec <- function(plot_data, label_map, sc) {

  group_pal <- .exercise_group_palette()
  tissue_pal <- MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS

  resolve_grp <- function(v) {
    v <- dplyr::recode(as.character(v), !!!label_map)
    if (v %in% names(group_pal)) unname(group_pal[v]) else "grey85"
  }
  resolve_tissue <- function(v) {
    key <- tolower(as.character(v))
    if (key %in% names(tissue_pal)) unname(tissue_pal[key]) else "grey85"
  }

  # The single-feature panels key the ome on `assay`; the logFC panels use
  # `platform`.
  ome_col <- if ("assay" %in% names(plot_data)) plot_data$assay else plot_data$platform
  omes_present <- unique(ifelse(grepl("metab", ome_col), "metab", as.character(ome_col)))
  single_ome <- length(omes_present) == 1
  row_var <- if (single_ome) "tissue" else "tissue_assay"

  col_levels <- levels(droplevels(as.factor(plot_data$randomGroupCode)))
  row_levels <- levels(droplevels(as.factor(plot_data[[row_var]])))
  col_cols <- vapply(col_levels, resolve_grp, character(1))
  row_cols <- if (single_ome) {
    vapply(row_levels, resolve_tissue, character(1))
  } else {
    rep("grey85", length(row_levels))
  }

  list(
    strip = .themed_strip(col_cols, row_cols, text_size = 12 * sc),
    facet_formula = stats::as.formula(paste0(row_var, " ~ randomGroupCode")),
    title = function(feature_label) {
      if (single_ome) paste0(feature_label, " — ", omes_present) else feature_label
    }
  )
}

# ---- single-feature trajectory ---------------------------------------------

#' Mean normalized abundance over time for one feature, female / male / overall.
#'
#' Pulls the qc_norm matrices for the requested tissues, keeps the rows matching
#' the feature, and summarises the mean and its 95% t confidence interval within
#' each tissue x assay x sex x exercise group x timepoint cell. The overall
#' (both-sex) series is summarised from the same sample-level values, so it is
#' the pooled mean rather than an average of the two sex means.
#'
#' @param feature A feature_id, gene symbol, or RefMet name.
#' @param selected_tissues "all", or any subset of adipose / blood / muscle.
#' @param scale_factor Multiplies every text, line and point size.
#' @param include_legend Whether to draw the female / male / overall legend.
#' @param legend_position Passed to ggplot2::theme(legend.position).
#' @param verbose Whether to warn when the feature maps to several labels.
#'
#' @return A ggplot.
sex_differences_single_feature <- function(feature,
                                           selected_tissues = c("all", "adipose",
                                                                "blood", "muscle"),
                                           scale_factor = 1,
                                           include_legend = TRUE,
                                           legend_position = "right",
                                           verbose = TRUE) {

  selected_tissues <- match.arg(selected_tissues, several.ok = TRUE)
  if ("all" %in% selected_tissues) {
    selected_tissues <- MotrpacHumanPreSuspensionAnalysis::tissue_available_list()
  }

  feature_info <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    dplyr::filter(tolower(feature_id) == tolower(feature) |
                    tolower(gene_symbol) == tolower(feature) |
                    tolower(refmet_name) == tolower(feature)) %>%
    dplyr::semi_join(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE,
                     by = "gene_symbol") %>%
    dplyr::mutate(feature_id = dplyr::case_when(
      !is.na(refmet_name) ~ refmet_name,
      TRUE ~ feature_id
    ))

  gene_symbol_options <- as.character(unique(feature_info$gene_symbol))
  refmet_options <- as.character(unique(feature_info$refmet_name))
  label_options <- union(gene_symbol_options, refmet_options)
  label_options <- label_options[!is.na(label_options)]
  feature_label <- label_options[1]

  if (length(label_options) > 1 && verbose) {
    message("several gene symbols or refmet names match this input; the first ",
            "is used as the label")
  }

  qc_data <- MotrpacHumanPreSuspensionData::load_qc(
    selected_tissues = selected_tissues
  )

  pheno_data <- MotrpacHumanPreSuspensionData::pheno$data %>%
    dplyr::select(vialLabel, Sex, randomGroupCode, Timepoint)

  individual_data <- lapply(names(qc_data), function(tissue_name) {
    lapply(names(qc_data[[tissue_name]]), function(ome_name) {
      mat <- qc_data[[tissue_name]][[ome_name]][["qc_norm"]]
      if (is.null(mat)) return(NULL)

      feature_rows <- which(rownames(mat) %in% feature_info$feature_id)
      if (length(feature_rows) == 0) return(NULL)

      mat[feature_rows, , drop = FALSE] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("feature_id") %>%
        tidyr::pivot_longer(-feature_id, names_to = "vialLabel",
                            values_to = "value") %>%
        dplyr::inner_join(pheno_data, by = "vialLabel") %>%
        dplyr::mutate(tissue = tissue_name, assay = ome_name)
    }) %>%
      dplyr::bind_rows()
  }) %>%
    dplyr::bind_rows()

  if (nrow(individual_data) == 0) {
    stop("no data matches '", feature, "'; it must appear in the feature_id or ",
         "gene_symbol column of HUMAN_FEATURE_TO_GENE, or - for metabolites - ",
         "in refmet_name", call. = FALSE)
  }

  # Mean sample value and its 95% t confidence interval within one set of
  # grouping columns. A single-observation cell has zero degrees of freedom, so
  # qt() returns NaN and that point carries no interval.
  summarise_series <- function(group_cols) {
    individual_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::summarize(
        Mean = mean(value, na.rm = TRUE),
        SD = stats::sd(value, na.rm = TRUE),
        Count = sum(!is.na(value)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        SE = SD / sqrt(Count),
        CI_95 = stats::qt((1 + 0.95) / 2, Count - 1),
        CI_low = Mean - CI_95 * SE,
        CI_high = Mean + CI_95 * SE
      )
  }

  # One summary cell is a tissue x assay x feature x group x timepoint. The two
  # sex series split each cell further by Sex, next to the feature it belongs to.
  cell_cols <- c("tissue", "assay", "feature_id", "randomGroupCode", "Timepoint")
  sex_cell_cols <- append(cell_cols, "Sex", after = match("feature_id", cell_cols))

  sex_summary <- summarise_series(sex_cell_cols) %>%
    dplyr::mutate(series = Sex)

  # Summarised from the same sample-level values as the two sex series, so this
  # is the pooled mean and not the average of the two sex means.
  overall_summary <- summarise_series(cell_cols) %>%
    dplyr::mutate(series = "Overall")

  # Listed in sampling order, which is also the x-axis order.
  timepoint_recode <- c(
    "pre_exercise"      = "Pre",
    "during_20_min"     = "D20M",
    "during_40_min"     = "D40M",
    "post_10_min"       = "P10M",
    "post_15_30_45_min" = "P15-45M",
    "post_3.5_4_hr"     = "P3.5/4H",
    "post_24_hr"        = "P24H"
  )
  timepoint_levels <- unname(timepoint_recode)

  plot_data <- dplyr::bind_rows(sex_summary, overall_summary) %>%
    dplyr::filter(randomGroupCode %in% c("ADUEndur", "ADUResist", "ADUControl")) %>%
    dplyr::mutate(
      tissue = stringr::str_to_sentence(tissue),
      assay = ifelse(grepl("metab", assay), "metab", assay),
      Timepoint = dplyr::recode(as.character(Timepoint), !!!timepoint_recode),
      Timepoint = factor(Timepoint, levels = timepoint_levels),
      tissue_assay = stringr::str_c(tissue, " ", assay),
      series = factor(series, levels = c("Female", "Male", "Overall"))
    )

  series_colors <- c(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_SEX_COLORS[c("Female", "Male")],
    "Overall" = "black"
  )
  label_map <- c("ADUControl" = "CON", "ADUEndur" = "EE", "ADUResist" = "RE")

  sc <- scale_factor * 0.7

  strip_spec <- .sex_feature_strip_spec(plot_data, label_map, sc)

  g <- ggplot(plot_data, aes(x = Timepoint, color = series, group = series)) +
    geom_line(aes(y = Mean), linewidth = 0.5 * sc) +
    geom_point(aes(y = Mean), size = 1.8 * sc) +
    geom_errorbar(
      aes(ymin = CI_low, ymax = CI_high),
      width = 0.3 * sc,
      linewidth = 0.4 * sc,
      alpha = 0.6
    ) +
    scale_color_manual(values = series_colors) +
    ggh4x::facet_grid2(strip_spec$facet_formula,
                       scales = "free_y",
                       labeller = labeller(randomGroupCode = label_map),
                       strip = strip_spec$strip) +
    ggtitle(strip_spec$title(feature_label)) +
    ylab("log2(normalized value)") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1, size = 9 * sc, color = "black"),
      axis.text.y = element_text(size = 10 * sc, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 11 * sc),
      plot.title = element_text(size = 12 * sc, margin = margin(0, 0, 0, 0)),
      strip.text.x = element_text(size = 12 * sc,
                                  margin = margin(0.11 * sc, 0, 0.11 * sc, 0, "cm")),
      strip.text.y = element_text(size = 12 * sc,
                                  margin = margin(0, 0.11 * sc, 0, 0.11 * sc, "cm")),
      legend.text = element_text(size = 9 * sc, margin = margin(0, 0, 0, 0)),
      legend.spacing.x = unit(0.01 * sc, "in"),
      legend.spacing.y = unit(0.1 * sc, "in"),
      legend.box.spacing = unit(0.01 * sc, "in"),
      legend.key.size = unit(0.1 * sc, "in"),
      legend.title = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      axis.ticks = element_line(linewidth = 0.3 * sc),
      panel.grid = element_line(linewidth = 0.3 * sc),
      panel.border = element_rect(linewidth = 0.3 * sc),
      plot.margin = margin(0.05 * sc, 0.05 * sc, 0.05 * sc, 0.05 * sc, "in")
    )

  if (include_legend && !is.null(legend_position)) {
    g <- g + theme(legend.position = legend_position)
  }

  g
}

# ============================================================================
# Cellular deconvolution — ED3D, ED3E and analysis/02_celltype_da.R
# ============================================================================
#
# The deconvolution itself happens outside this repository: CIBERSORTx is a web
# service (cibersortx.stanford.edu) with no R implementation, so its result CSV
# is an external input. It is individual-level data, is not distributed here and
# may not be published; see the note at the top of ED3.R. What lives here is
# everything downstream of that CSV that more than one place in the analysis
# needs: collapsing the 22 LM22 leukocyte subsets into the coarse populations the
# analysis reasons about, and the CIBERSORTx-vs-CBC neutrophil scatter.
#
# Every call in this half is namespace-qualified so it behaves the same whether
# or not the sourcing script has attached dplyr, tidyr or ggplot2.

# ---- LM22 subset -> coarse population --------------------------------------

#' Mapping from the 22 LM22 leukocyte subsets to eight coarse populations.
#'
#' LM22 resolves subsets whose marker profiles are highly correlated (resting vs
#' activated dendritic cells, M0/M1/M2 macrophages), which is finer than any
#' clinical measurement the estimates are compared against and is where a linear
#' deconvolution is least stable. Summing the collinear subsets into a coarser
#' population is what the analysis reports.
#'
#' Names are the LM22 column labels in `make.names()` form, so a CIBERSORTx CSV
#' read with either `check.names = TRUE` ("T.cells.CD8") or `check.names = FALSE`
#' ("T cells CD8") matches. Eosinophils are deliberately absent: they are not
#' assigned to a coarse population and therefore pass through as a carried
#' column rather than contributing to the rescaled totals.
DECONV_LM22_CELL_TYPE_MAP <- c(
  # T cells
  "T.cells.CD4.memory.activated" = "T cells",
  "T.cells.CD4.memory.resting"   = "T cells",
  "T.cells.CD4.naive"            = "T cells",
  "T.cells.CD8"                  = "T cells",
  "T.cells.follicular.helper"    = "T cells",
  "T.cells.gamma.delta"          = "T cells",
  "T.cells.regulatory..Tregs."   = "T cells",

  # B cells
  "B.cells.memory" = "B cells",
  "B.cells.naive"  = "B cells",
  "Plasma.cells"   = "B cells",

  # NK cells
  "NK.cells.activated" = "NK cells",
  "NK.cells.resting"   = "NK cells",

  # Monocytes
  "Monocytes" = "Monocytes",

  # Macrophages
  "Macrophages.M0" = "Macrophages",
  "Macrophages.M1" = "Macrophages",
  "Macrophages.M2" = "Macrophages",

  # Dendritic cells
  "Dendritic.cells.activated" = "Dendritic cells",
  "Dendritic.cells.resting"   = "Dendritic cells",

  # Neutrophils; mapped to itself so it goes through the same rescaling as the
  # collapsed populations and comes out as a percentage.
  "Neutrophils" = "Neutrophils",

  # Mast cells
  "Mast.cells.activated" = "Mast cells",
  "Mast.cells.resting"   = "Mast cells"
)

#' Collapse LM22 subset fractions into coarse cell populations.
#'
#' Sums the subsets of each coarse population, rescales every sample to sum to
#' 100, and drops populations averaging 1% or less across samples. Columns that
#' are not LM22 subsets (the `Mixture` sample identifier, and Eosinophils, which
#' `DECONV_LM22_CELL_TYPE_MAP` does not assign) are carried through untouched;
#' they are excluded from the rescaling denominator, so their share is
#' redistributed across the populations that remain.
#'
#' @param merged_proportions Data frame of CIBERSORTx fractions with a `Mixture`
#'   column and one column per LM22 subset. The CIBERSORTx run-quality columns
#'   (P-value, Correlation, RMSE) must already be dropped.
#' @return One row per mixture: the carried columns, then one column per
#'   retained population holding its percentage of the sample.
squish_proportions <- function(merged_proportions) {
  if (!"Mixture" %in% colnames(merged_proportions)) {
    stop("squish_proportions() needs the CIBERSORTx 'Mixture' column; ",
         "call it before renaming that column to vialLabel", call. = FALSE)
  }

  subset_key <- make.names(colnames(merged_proportions))
  is_subset <- subset_key %in% names(DECONV_LM22_CELL_TYPE_MAP)
  unique_groups <- unique(unname(DECONV_LM22_CELL_TYPE_MAP))

  group_sums <- vapply(unique_groups, function(grp) {
    members <- names(DECONV_LM22_CELL_TYPE_MAP)[DECONV_LM22_CELL_TYPE_MAP == grp]
    cols <- colnames(merged_proportions)[subset_key %in% members]
    rowSums(as.data.frame(merged_proportions)[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(merged_proportions)))
  group_sums <- as.data.frame(group_sums, check.names = FALSE)

  summed_grouped <- dplyr::bind_cols(
    merged_proportions[, !is_subset, drop = FALSE],
    group_sums
  ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(unique_groups),
      names_to = "cell_group",
      values_to = "proportion"
    ) |>
    dplyr::group_by(Mixture) |>
    dplyr::mutate(proportion = 100 * proportion / sum(abs(proportion))) |>
    dplyr::ungroup()

  # A population averaging 1% or less of the sample is below the resolution the
  # deconvolution can be trusted at, and is dropped rather than plotted.
  retained <- summed_grouped |>
    dplyr::group_by(cell_group) |>
    dplyr::summarise(cell_group_avg = mean(proportion), .groups = "drop") |>
    dplyr::filter(cell_group_avg > 1) |>
    dplyr::pull(cell_group)

  summed_grouped |>
    dplyr::filter(cell_group %in% retained) |>
    tidyr::pivot_wider(names_from = "cell_group", values_from = "proportion")
}

# ---- CIBERSORTx vs clinical CBC --------------------------------------------

#' Scatter of CIBERSORTx neutrophil estimates against measured CBC neutrophils.
#'
#' One point per joined pre-exercise sample: the clinical complete blood count
#' neutrophil value on x, the CIBERSORTx neutrophil percentage on y, with an OLS
#' fit and the Pearson correlation annotated. The correlation is computed after
#' the `filter_lt5` exclusion and after the join, so it describes exactly the
#' points drawn. Both quantities are percentages of the same samples, so the two
#' axes share one range and one scale and the panel is square.
#'
#' The CBC table carries two visits — ADU_SCP (screening) and ADU_FUP
#' (follow-up) — and 47 participants appear in both, with values that differ
#' because the draws are months apart. The join key is `pid` alone, so keeping
#' both would plot those participants twice and weight them double in the
#' correlation. Only the screening visit is kept: it is one row per participant
#' by construction and it is the baseline panel.
#'
#' @param blood_metadata Blood RNA-seq sample metadata joined to collapsed
#'   CIBERSORTx proportions; needs `pid`, `Timepoint` and `Neutrophils`.
#' @param cbc_raw The `cln_raw_local_lab_results$data` table; needs `visit_code`.
#' @param filter_lt5 Drop CBC rows whose neutrophil value is 5 or less. Two rows
#'   carry implausibly low clinical neutrophil values (3.6 and 4.1).
plot_cbc_vs_cibersortx <- function(blood_metadata, cbc_raw, filter_lt5 = TRUE) {
  cbc <- cbc_raw |>
    dplyr::filter(.data$visit_code == "ADU_SCP") |>
    dplyr::mutate(
      Neutrophils = neutro_labr,
      Monocytes   = mono_labr,
      Lymphocytes = lymp_labr
    ) |>
    dplyr::select(pid, Neutrophils, Monocytes, Lymphocytes) |>
    dplyr::mutate(pid = as.character(pid))

  # One row per participant is what makes `pid` a valid join key. If a future
  # release repeats a screening visit, stop rather than silently double-count.
  if (anyDuplicated(cbc$pid)) {
    stop("the screening CBC carries more than one row for ",
         sum(duplicated(cbc$pid)), " participant(s); `pid` is no longer a key",
         call. = FALSE)
  }

  if (filter_lt5) {
    cbc <- cbc |> dplyr::filter(Neutrophils > 5)
  }

  neutrophil_corr <- blood_metadata |>
    dplyr::filter(Timepoint == "pre_exercise") |>
    dplyr::mutate(pid = as.character(pid)) |>
    dplyr::select(pid, Neutrophils) |>
    dplyr::rename(deconv_neutrophils = Neutrophils) |>
    dplyr::inner_join(cbc |> dplyr::select(pid, Neutrophils), by = "pid") |>
    dplyr::rename(cbc_neutrophils = Neutrophils)

  r_val <- stats::cor(neutrophil_corr$deconv_neutrophils,
                      neutrophil_corr$cbc_neutrophils,
                      use = "complete.obs")

  # Both axes are neutrophil percentages of the same samples, so they are drawn
  # over one shared range at ratio 1. How far a point sits from the identity
  # line is then how far the estimate sits from the measurement, rather than a
  # property of two independently scaled axes.
  axis_lims <- range(c(neutrophil_corr$cbc_neutrophils,
                       neutrophil_corr$deconv_neutrophils), na.rm = TRUE)

  neutrophil_corr |>
    ggplot2::ggplot(ggplot2::aes(x = cbc_neutrophils, y = deconv_neutrophils)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5,
                      label = paste0("R = ", round(r_val, 2)), size = 4) +
    ggplot2::coord_fixed(ratio = 1, xlim = axis_lims, ylim = axis_lims) +
    ggplot2::labs(
      x = "CBC Neutrophils",
      y = "CIBERSORTx Neutrophils (%)",
      title = "Neutrophil estimates: CIBERSORTx vs. CBC (pre-exercise)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      title      = ggplot2::element_text(size = 13)
    )
}
