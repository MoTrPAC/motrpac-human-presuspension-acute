#!/usr/bin/env Rscript
# Figure 3 — Resistance vs endurance exercise
#
# Panels:  FIG3A    logFC RE vs EE per tissue x ome, with UpSet bars
#          FIG3B    muscle transcriptomics pathway enrichment, RE vs EE
#                                          (also writes ST3a, ST3b and ST3c)
#          FIG3C    features responding in opposite directions to EE and RE
#          FIG3D    CAR 6:0, BCKDHB, CBFA2T3 and LRRC14B trajectories
#          FIG3EFG  VO2peak association and acute response, EE-specific,
#                   RE-specific and shared transcripts
# Tables:  ST3a     EE muscle-only ORA features     (written by FIG3B)
#          ST3b     RE muscle-only ORA features     (written by FIG3B)
#          ST3c     all-muscle RE/EE-only ORA       (written by FIG3B)
#
# Needs consortium data access. FIG3EFG additionally reads the baseline clinical
# x omics fit, which is not in this repository and is not cheap: run
#
#   Rscript figures/landscape/analysis/03_clinical_omics.R
#
# once and it is cached. Without it FIG3EFG reports SKIPPED and the other four
# panels still build.
#
#   Rscript figures/landscape/FIG3.R          every panel
#   Rscript figures/landscape/FIG3.R FIG3B    one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tidyr)
  # Attached, not just namespace-qualified: load_differential_analysis()
  # resolves the per-tissue DA objects by name off the search path.
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "..", "lib", "single_feature_helpers.R"))
source(file.path(here, "..", "single_feature_plots.R"))
source(file.path(here, "FIG3_ED4_helpers.R"))

# ---- shared ----------------------------------------------------------------

# FIG3A and FIG3B are two readings of one RE-vs-EE feature table, built with the
# same arguments, so one load serves both. Lazy and memoised rather than loaded
# at the top of the script: the DA tables should not be read to build FIG3D.
re_vs_ee_features <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- re_vs_ee_feature_table()
    cache
  }
})

# opposite_ee_re_table() is NOT shared with ED4F: FIG3C calls it with
# selected_tissues = "all" and ED4F with "muscle", so the two are different
# loads, and they are in different scripts besides.

# ---- FIG3A — logFC resistance vs endurance, per tissue x ome ---------------

# Top row: for every feature, logFC of the resistance arm against its control
# plotted against logFC of the endurance arm against its control, one facet per
# tissue. Features significant in the direct RE-vs-EE contrast are drawn as
# filled circles coloured by which single-group contrasts they are also
# significant in; the rest are grey. Bottom row: an UpSet plot per tissue over
# the same three significance calls, with the intersection bars coloured to
# match.

fig3a <- function() {
  panel_init("FIG3A")

  categorized <- re_vs_ee_features()

  # "Only EE" / "Only RE" take the exercise-group palette; the two categories
  # that are not a single exercise group get their own colours.
  category_colors <- c(
    "Only EE"     = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[["ADUEndur"]],
    "Only RE"     = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[["ADUResist"]],
    "Both Groups" = "brown",
    "Neither"     = "#86C5D8"
  )

  # Paint one UpSet intersection in a category colour. Only the matrix dots and
  # the size bar are recoloured; every other component keeps its ComplexUpset
  # default.
  color_intersection <- function(sets, category) {
    upset_query(
      intersect = sets,
      color = category_colors[[category]],
      fill  = category_colors[[category]],
      only_components = c("intersections_matrix", "Intersection size")
    )
  }

  not_re_ee_hit <- dplyr::filter(categorized, sig_contrast != "Sig")
  re_ee_hit     <- dplyr::filter(categorized, sig_contrast == "Sig")

  scatter <- ggplot() +
    # Diagonal and axes, so the quadrants read without a grid.
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    maybe_rasterise(
      geom_point(
        data = not_re_ee_hit,
        aes(x = logFC_ADUEndur, y = logFC_ADUResist),
        color = "grey", size = 0.5, shape = 20
      ),
      nrow(not_re_ee_hit), "FIG3A"
    ) +
    maybe_rasterise(
      geom_point(
        data = re_ee_hit,
        aes(x = logFC_ADUEndur, y = logFC_ADUResist, fill = Significant_in),
        color = "black", size = 1.5, shape = 21, stroke = 0.4
      ),
      nrow(re_ee_hit), "FIG3A"
    ) +
    scale_fill_manual(values = category_colors) +
    labs(
      x = "logFC (EE vs Control)",
      y = "logFC (RE vs Control)",
      title = NULL
    ) +
    facet_grid(~ tissue) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 12)
    )

  # One UpSet per tissue, in the tissue order the feature table is sorted into.
  upset_list <- lapply(unique(categorized$tissue), function(tissue_i) {
    df_tiss <- categorized |>
      dplyr::filter(tissue == tissue_i) |>
      dplyr::select(feature_id, assay, sig_contrast, sig_ee_single,
                    sig_re_single) |>
      dplyr::distinct() |>
      dplyr::mutate(dplyr::across(
        c(sig_contrast, sig_ee_single, sig_re_single), ~ . == "Sig"
      )) |>
      dplyr::rename(
        EE_vs_RE    = sig_contrast,
        EE_vs_Ctrl  = sig_ee_single,
        RE_vs_Ctrl  = sig_re_single
      ) |>
      # Drop features significant in none of the three comparisons.
      dplyr::filter(EE_vs_RE | EE_vs_Ctrl | RE_vs_Ctrl)

    # upset_query() errors on an intersection a tissue has no members for, so
    # the triple intersection is only queried where it exists.
    custom_colors <- list(
      color_intersection(c("EE_vs_RE", "EE_vs_Ctrl"), "Only EE"),
      color_intersection(c("EE_vs_RE", "RE_vs_Ctrl"), "Only RE"),
      color_intersection("EE_vs_RE", "Neither")
    )

    if (tissue_i == "muscle" || tissue_i == "blood") {
      # Appended by length rather than at a fixed index, so adding an entry to
      # the list above cannot silently overwrite this one.
      custom_colors[[length(custom_colors) + 1L]] <- color_intersection(
        c("EE_vs_RE", "EE_vs_Ctrl", "RE_vs_Ctrl"), "Both Groups"
      )
    }

    ComplexUpset::upset(
      df_tiss,
      intersect = c("EE_vs_RE", "EE_vs_Ctrl", "RE_vs_Ctrl"),
      # Intersections are listed explicitly and neither they nor the sets are
      # re-sorted, so every tissue uses the same column order.
      sort_intersections = FALSE,
      sort_sets = FALSE,
      intersections = list(
        c("EE_vs_RE", "EE_vs_Ctrl", "RE_vs_Ctrl"),
        c("RE_vs_Ctrl", "EE_vs_RE"),
        c("EE_vs_Ctrl", "EE_vs_RE"),
        "EE_vs_RE",
        "RE_vs_Ctrl",
        "EE_vs_Ctrl"
      ),
      name = tissue_i,
      encode_sets = FALSE,
      min_size = 3,
      set_sizes = FALSE,
      base_annotations = list(
        # The assay fill legend is dropped here so it is not repeated per tissue.
        "Intersection size" = intersection_size(counts = TRUE) +
          guides(fill = "none")
      ),
      matrix = intersection_matrix(
        geom = geom_point(shape = "square", size = 3.5),
        segment = geom_segment(linetype = "dotted")
      ),
      queries = custom_colors
    )
  })

  combined_upset <- patchwork::wrap_plots(upset_list, nrow = 1)

  final_plot <- scatter / combined_upset +
    plot_layout(heights = c(1.5, 1.2))

  export_panel(function() print(final_plot), "FIG3A")
}

# ---- FIG3B — muscle transcriptomics pathway enrichment, RE vs EE -----------

# Over-representation analysis of the features that are significant in the
# RE-vs-EE contrast and in exactly one exercise arm. ORA is run per tissue x
# assay x arm against that tissue x assay's own measured genes as background;
# only muscle transcriptomics clears the 50-feature minimum. The top 15
# signatures per arm are drawn as a bubble heatmap, coloured by -log10 p-value
# and shaded by BH-adjusted significance.
#
# This panel also writes three supplementary sub-tables, all of them by-products
# of the loop below rather than separate computations:
#
#   ST3a  the muscle transcriptomics features significant in EE only — the ORA's
#         own input list for that arm
#   ST3b  the same for RE
#   ST3c  the accumulated ORA results the panel plots the top 15 of
#
# Writing them anywhere else would mean a second pass over the same fit, and the
# panel plots a slice of ST3c, so the two could disagree about a p-value.

fig3b <- function() {
  panel_init("FIG3B")

  categorized <- re_vs_ee_features()

  # ORA per tissue x assay x arm. Metabolomics and phosphoproteomics are
  # skipped: their signature databases (RefMet, PTMsigDB) are not gene sets.
  all_ora_results <- data.frame()

  # The per-arm input lists ST3a and ST3b report, captured as the loop builds
  # them. Muscle transcriptomics only, which is what those two sub-tables are.
  ora_inputs <- list()

  for (tissue_i in MotrpacHumanPreSuspensionAnalysis::tissue_available_list()) {
    assays_i <- categorized |>
      dplyr::filter(tissue == tissue_i) |>
      dplyr::pull(assay) |>
      unique()

    for (assay_i in assays_i) {
      if (assay_i == "metab" || assay_i == "prot-ph") next

      for (group_i in c("Only EE", "Only RE")) {
        contrast_i <- if (group_i == "Only EE") "EE" else "RE"

        # Every gene measured in this tissue x assay is the ORA background.
        all_results <- categorized |>
          dplyr::filter(tissue == tissue_i, assay == assay_i) |>
          dplyr::left_join(
            MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE,
            by = c("feature_id", "assay")
          ) |>
          dplyr::filter(!is.na(gene_symbol))

        sig_results <- all_results |>
          dplyr::filter(Significant_in == group_i, adj_p_value < 0.1)

        message(paste(tissue_i, group_i, assay_i, nrow(sig_results)))

        # Too few hits to over-represent anything.
        if (nrow(sig_results) < 50) next

        if (tissue_i == "muscle" && assay_i == "transcript-rna-seq") {
          ora_inputs[[contrast_i]] <- sig_results
        }

        significant <- as.character(unique(sig_results$gene_symbol))
        all_feats   <- as.character(unique(all_results$gene_symbol))

        ora_results <- MotrpacHumanPreSuspensionAnalysis::run_ORA(
          input = significant,
          background = all_feats
        )

        all_ora_results <- rbind(
          all_ora_results,
          dplyr::mutate(
            ora_results,
            tissue = tissue_i,
            assay = assay_i,
            contrast = contrast_i
          )
        )
      }
    }
  }

  all_ora_results$contrast <- as.factor(all_ora_results$contrast)

  write_st3a(ora_inputs)
  write_st3b(ora_inputs)
  write_st3c(all_ora_results)

  # Top 15 signatures per arm, then the database prefix stripped off the label.
  ora_figure <- all_ora_results |>
    dplyr::filter(tissue == "muscle") |>
    dplyr::group_by(contrast) |>
    dplyr::arrange(adj_p_value, .by_group = TRUE) |>
    dplyr::slice_head(n = 15) |>
    dplyr::ungroup() |>
    dplyr::filter(adj_p_value < 0.05) |>
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      set_short = truncate_set_label(
        stringr::str_remove(set, paste0("^", database, "_"))
      )
    )

  # Rows are drawn in the order they were selected, arm by arm, rather than
  # alphabetically.
  ordered_sets <- dplyr::pull(ora_figure, set_short)

  modality_colors <- c(
    "EE" = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[["ADUEndur"]],
    "RE" = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[["ADUResist"]]
  )

  dotplot <- function() {
    TMSig::enrichmap(
      x = ora_figure,
      n_top = 30L,
      set_column = "set_short",
      statistic_column = "neg_log10_p",
      contrast_column = "contrast",
      padj_column = "adj_p_value",
      colors = c("white", "#543483"),
      heatmap_args = list(
        column_title = "Muscle \n Transcriptomics",
        row_order = ordered_sets,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        na_col = "grey55",
        top_annotation = HeatmapAnnotation(
          Modality = c("EE", "RE"),
          col = list(Modality = modality_colors),
          show_annotation_name = TRUE
        ),
        heatmap_legend_param = list(
          title = latex2exp::TeX("$\\bf{$-log$_{10}($P-Value$)}$")
        )
      ),
      # enrichmap finishes with ComplexHeatmap::draw(), whose newpage default
      # starts a fresh page on the device export_panel() already opened - leaving
      # page 1 blank and the heatmap on page 2, which is the page Illustrator
      # discards.
      draw_args = list(newpage = FALSE)
    )
  }

  # cairo_pdf: the adjusted-p legend labels enrichmap draws use U+2265 ("≥"),
  # which base pdf() cannot encode under a Latin-1 Helvetica.
  export_panel(dotplot, "FIG3B")
}

# ---- ST3a, ST3b — the ORA input lists --------------------------------------

# One sub-table per arm, same schema. entrez_gene and gene_symbol are already on
# the frame: the ORA background join above attached them, and they are what the
# gene sets were matched on.
write_ora_input <- function(ora_inputs, arm, table_id) {
  input_rows <- ora_inputs[[arm]]
  if (is.null(input_rows)) {
    stop(table_id, ": no muscle transcriptomics input list for the ", arm,
         " arm — the ORA loop skipped it, so there is nothing to report.",
         call. = FALSE)
  }
  export_table(
    as.data.frame(input_rows)[, table_spec(table_id)$columns, drop = FALSE],
    table_id
  )
}

write_st3a <- function(ora_inputs) write_ora_input(ora_inputs, "EE", "ST3a")
write_st3b <- function(ora_inputs) write_ora_input(ora_inputs, "RE", "ST3b")

# ---- ST3c — the accumulated ORA --------------------------------------------

# Muscle only, matching what the panel plots and what the table reports. The
# other tissue x assay combinations never clear the 50-feature minimum, so this
# filter is a statement of scope rather than a reduction.
write_st3c <- function(all_ora_results) {
  st3c <- all_ora_results |>
    dplyr::filter(tissue == "muscle") |>
    # z.std is the two-sided p converted to a z score by the inverse normal cdf.
    # A p small enough to round to zero sends qnorm to Inf, so it is capped at 15 —
    # beyond that the score is reporting the floating-point floor, not the data.
    dplyr::mutate(
      z.std = stats::qnorm(1 - .data$p_value / 2),
      z.std = ifelse(is.infinite(.data$z.std), 15, .data$z.std),
      contrast = as.character(.data$contrast)
    ) |>
    as.data.frame()

  export_table(st3c[, table_spec("ST3c")$columns, drop = FALSE], "ST3c")
}

# ---- FIG3C — features responding in opposite directions to EE and RE -------

# Every feature that is significantly UP in one exercise arm and significantly
# DOWN in the other at the same timepoint, plotted as its endurance logFC
# against its resistance logFC. Point outline is the tissue, fill is the ome, and
# the labelled points are those furthest from agreement.
#
# The diagonal is what the panel is about: a feature on it responds the same way
# to both arms, and everything here is off it by construction — the second and
# fourth quadrants are the whole population.
#
# Ported from Opposite_EE_vs_RE.Rmd, which drew this as
# `opposite_ee_vs_re_all_per_tissue.pdf`. The frame and the label rule are in
# helpers/FIG3_ED4.R, shared with ED4F, which is the same Rmd's muscle-only
# version.

fig3c = function() {
  panel_init("FIG3C")

  opposite = opposite_ee_re_table(selected_tissues = "all")

  tissue_colors = lookup_colors(
  MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS,
  sort(unique(as.character(opposite$tissue))), "tissue"
  )

  ome_scale = ome_fill_scale(opposite$assay)

  # color and fill are mapped at the top level, not inside geom_point: the label
  # layer below recolours itself with after_scale(fill), which can only see an
  # aesthetic it inherits. The zero lines are drawn with inherit.aes = FALSE so
  # they do not pick up either scale.
  p = ggplot(opposite, aes(x = .data$logFC_Endur, y = .data$logFC_Resist,
                          color = .data$tissue, fill = .data$assay)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80",
               inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80",
               inherit.aes = FALSE) +
    geom_point(size = 3, shape = 21, stroke = 1) +
    geom_text_repel(
      aes(label = .data$plot_label, color = after_scale(fill)),
      size = 3.5, box.padding = 0.3, point.padding = 0.2,
      max.overlaps = Inf, max.iter = 10000, force = 2, force_pull = 1,
      min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6,
      show.legend = FALSE
    ) +
    scale_color_manual(name = "Tissue", values = tissue_colors) +
    scale_fill_manual(name = "Ome", values = ome_scale$colors,
                      labels = ome_scale$labels) +
    labs(
      title = "Differential EE vs RE responses",
      x = "log2FC endurance",
      y = "log2FC resistance"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 14))

  export_panel(p, "FIG3C")
}

# ---- FIG3D — CAR 6:0, BCKDHB, CBFA2T3 and LRRC14B trajectories -------------

# Four features whose responses separate between the resistance and endurance
# exercise groups, each drawn as the exercise-group mean with a 95% confidence
# interval across the acute timepoints. A point is filled black where that
# timepoint's adjusted p value is below 0.05. The four sub-plots sit two per row
# under one collected legend.
#
# The plot is one entry of the single-feature catalog in single_feature_plots.R.

fig3d <- function() {
  panel_init("FIG3D")
  export_panel(single_feature_plot("FIG3D"), "FIG3D")
}

# ---- FIG3EFG — VO2peak association and acute response ----------------------

# One figure in place of the three dumbbell panels e, f and g. Muscle transcripts
# whose baseline abundance is associated with VO2peak (FDR < 0.05 in the
# analysis/03_clinical_omics.R fit) and that respond to acute exercise in a
# category's pattern, capped at the 70 rows VO2MAX_MAX_ROWS allows per category.
# Columns are the three categories, rows the three post-exercise timepoints;
# within a block, one row per transcript, the EE and RE z statistics joined by a
# segment and the VO2peak association tiled beside it. One legend, one z axis and
# one VO2peak colour scale for all nine blocks.
#
# The categories are declared in VO2MAX_PANELS in helpers/FIG3_ED4.R, which is
# where the drawing lives.

fig3efg <- function() {
  # The manifest marks CLINICAL_OMICS_TABLE must_exist for this panel, so an
  # unrun fit is reported with the command that produces it rather than dying on
  # a missing path. A fit produced elsewhere therefore has to be named in
  # CLINICAL_OMICS_TABLE, not only in FIG3_CLINICAL_OMICS_TABLE, for this panel
  # to run.
  skip <- skip_if_no_fit("FIG3EFG")
  if (!is.null(skip)) return(skip)

  panel_init("FIG3EFG")

  export_panel(function() print(vo2max_dumbbell_combined()), "FIG3EFG")
}

run_panels(list(FIG3A = fig3a, FIG3B = fig3b, FIG3C = fig3c, FIG3D = fig3d,
                FIG3EFG = fig3efg))
