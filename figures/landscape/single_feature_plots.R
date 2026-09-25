#!/usr/bin/env Rscript
# single_feature_plots.R — every single-feature trajectory plot, in one place.
#
# A single-feature plot is one feature's abundance over the acute timepoints:
# the exercise-group mean with a 95% confidence interval, a point filled black
# where that timepoint's adjusted p value clears the threshold.
#
#   source()d   defines SINGLE_FEATURE_PLOTS and single_feature_plot(), and
#               draws nothing. This is what the figure scripts do: FIG2.R's
#               fig2d() is export_panel(single_feature_plot("FIG2D"), "FIG2D").
#
#   run         builds every entry through the export contract, so the sixteen
#               can be rebuilt without running ten figure scripts:
#
#                 Rscript figures/landscape/single_feature_plots.R
#                 Rscript figures/landscape/single_feature_plots.R ED2B FIG7C
#
# ED3B and ED3C are drawn by sex_differences_single_feature() from
# extended_data_3/ED3_helpers.R; every other entry uses plot_single_feature().

# Sourced by a figure script as well as run directly, so the root is taken from
# landscape_root() when panel_export.R is already loaded and from this file's
# own path when it is not.
if (!exists("landscape_root", mode = "function")) {
  here <- dirname(sub("^--file=", "",
                      grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  source(file.path(here, "lib", "panel_export.R"))
}
source(file.path(landscape_root(), "lib", "single_feature_helpers.R"))
source(file.path(landscape_root(), "lib", "highlights.R"))

# The caller attaches the Analysis package: plot_single_feature() calls
# load_differential_analysis(), which resolves its DA objects through the search
# path. panel_init() does it for panels, the direct-run block for the rest.
plot_single_feature <- MotrpacHumanPreSuspensionAnalysis::plot_single_feature

#' plot_single_feature() on one feature named in config/highlights.json.
#'
#' The feature, tissues and omes come from the file; plotting options are
#' passed through.
highlighted_feature_plot <- function(panel, name, ...) {
  hl <- highlight_feature(panel, name)
  plot <- plot_single_feature(feature = hl$feature,
                              selected_tissues = hl$tissues,
                              selected_omes = hl$omes, ...)
  if (!is.null(hl$title)) plot <- plot + ggplot2::ggtitle(hl$title)
  plot
}

# ---- The catalog -----------------------------------------------------------
#
#   panel       the id in config/panel_map.json, or NA.
#   figure      where it sits in the manuscript.
#   title       what it shows.
#   build       a function of no arguments returning the ggplot or patchwork.
#
# Sizes are not here: export_panel() reads each panel's from the manifest.

SINGLE_FEATURE_PLOTS <- list(

  # ---- Figure 2 ------------------------------------------------------------

  FIG2D = list(
    panel = "FIG2D",
    figure = "Figure 2D",
    title = "MYC and AREG transcript trajectories, all three tissues",
    build = function() {
      assemble_feature_panel(
        list(
          highlighted_feature_plot("FIG2D", "myc", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("FIG2D", "areg", color_time_labels = FALSE,
                                   include_legend = TRUE)
        ),
        ncol = 1,
        legend_position = "right"
      )
    }
  ),

  # ---- Extended Data 2 -----------------------------------------------------

  ED2B = list(
    panel = "ED2B",
    figure = "Extended Data 2B",
    title = "Example single-feature trajectories per hit type",
    build = function() {
      assemble_feature_panel(
        list(
          highlighted_feature_plot("ED2B", "ppargc1a", color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("ED2B", "per2", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("ED2B", "alloisoleucine",
                                   color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("ED2B", "hspa4", color_time_labels = FALSE,
                                   include_legend = FALSE)
        ),
        ncol = 2,
        legend_position = "right"
      )
    }
  ),

  ED2D = list(
    panel = "ED2D",
    figure = "Extended Data 2D",
    title = "NR4A1 and NFKBIZ transcript trajectories, all three tissues",
    build = function() {
      assemble_feature_panel(
        list(
          highlighted_feature_plot("ED2D", "nr4a1", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("ED2D", "nfkbiz", color_time_labels = FALSE,
                                   include_legend = TRUE)
        ),
        ncol = 1,
        legend_position = "right"
      )
    }
  ),

  # ---- Extended Data 3 -----------------------------------------------------
  #
  # ED3_helpers.R is sourced per entry, not at the top of the file:
  # nothing else in the catalog needs it.

  ED3B = list(
    panel = "ED3B",
    figure = "Extended Data 3B",
    title = "PPARGC1A muscle transcript trajectory, split by sex",
    build = function() {
      source(file.path(landscape_root(), "extended_data_3", "ED3_helpers.R"))
      hl <- highlight_feature("ED3B", "ppargc1a")
      sex_differences_single_feature(hl$feature, selected_tissues = hl$tissues)
    }
  ),

  ED3C = list(
    panel = "ED3C",
    figure = "Extended Data 3C",
    title = "CCN1 transcript trajectory in muscle and adipose, split by sex",
    build = function() {
      source(file.path(landscape_root(), "extended_data_3", "ED3_helpers.R"))
      # CCN1 is named by Ensembl id: the symbol also matches ATAC peaks and
      # methylation sites, which would add rows.
      hl <- highlight_feature("ED3C", "ccn1")
      sex_differences_single_feature(hl$feature, selected_tissues = hl$tissues)
    }
  ),

  # ---- Figure 3 ------------------------------------------------------------

  FIG3D = list(
    panel = "FIG3D",
    figure = "Figure 3D",
    title = "Feature trajectories diverging between RE and EE",
    build = function() {
      assemble_feature_panel(
        list(
          highlighted_feature_plot("FIG3D", "car_6_0", color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG3D", "bckdhb", color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG3D", "cbfa2t3", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("FIG3D", "lrrc14b", color_time_labels = FALSE,
                                   include_legend = TRUE)
        ),
        ncol = 2,
        legend_position = "bottom"
      )
    }
  ),

  # ---- Figure 4 ------------------------------------------------------------

  FIG4I = list(
    panel = "FIG4I",
    figure = "Figure 4I",
    title = "TAMALIN transcript trajectory in adipose, blood and muscle",
    build = function() {
      # color_time_labels and include_legend stay at their defaults, as in the
      # legacy chunk.
      assemble_feature_panel(
        list(
          highlighted_feature_plot("FIG4I", "adipose"),
          highlighted_feature_plot("FIG4I", "blood"),
          highlighted_feature_plot("FIG4I", "muscle")
        ),
        ncol = 3,
        legend_position = "right"
      )
    }
  ),

  # ---- Extended Data 4 -----------------------------------------------------

  ED4D = list(
    panel = "ED4D",
    figure = "Extended Data 4D",
    title = "HSPH1 muscle transcript trajectory",
    build = function() {
      highlighted_feature_plot("ED4D", "hsph1", color_time_labels = FALSE,
                               include_legend = TRUE)
    }
  ),

  ED4E = list(
    panel = "ED4E",
    figure = "Extended Data 4E",
    title = "Diverging feature grid, muscle",
    build = function() {
      assemble_feature_panel(
        list(
          highlighted_feature_plot("ED4E", "fhod1", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("ED4E", "myog", color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("ED4E", "zbtb4", color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("ED4E", "nol3", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("ED4E", "myot", color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("ED4E", "lncrna", color_time_labels = FALSE,
                                   include_legend = FALSE)
        ),
        ncol = 2,
        legend_position = "bottom"
      )
    }
  ),

  # ---- Figure 6 ------------------------------------------------------------

  FIG6F = list(
    panel = "FIG6F",
    figure = "Figure 6F",
    title = "Lactic acid across the three tissues",
    build = function() {
      # color_time_labels and include_legend stay at their defaults, as in the
      # legacy chunk.
      highlighted_feature_plot("FIG6F", "lactic_acid")
    }
  ),

  # ---- Extended Data 8 -----------------------------------------------------

  ED8B = list(
    panel = "ED8B",
    figure = "Extended Data 8B",
    title = "TFEB transcript across the three tissues",
    build = function() {
      highlighted_feature_plot("ED8B", "tfeb")
    }
  ),

  ED8D = list(
    panel = "ED8D",
    figure = "Extended Data 8D",
    title = "HSP90AA1 transcript across the three tissues",
    build = function() {
      highlighted_feature_plot("ED8D", "hsp90aa1")
    }
  ),

  ED8F = list(
    panel = "ED8F",
    figure = "Extended Data 8F",
    title = "Leucine and Leu-Ile in blood and muscle",
    build = function() {
      # The legacy wrote these as two PDFs to be placed side by side; one
      # panel here, as FIG5G is.
      assemble_feature_panel(
        list(
          highlighted_feature_plot("ED8F", "leucine"),
          highlighted_feature_plot("ED8F", "leu_ile")
        ),
        ncol = 2,
        legend_position = "right"
      )
    }
  ),

  # ---- Figure 7 ------------------------------------------------------------
  #
  # CCN1 also appears as ED3C, drawn the other way: sex-split, muscle and
  # adipose only.

  FIG7C = list(
    panel = "FIG7C",
    figure = "Figure 7C",
    title = "CCN1 across muscle and adipose transcript, blood Olink and muscle protein",
    build = function() {
      # legend_position on the protein plot has no effect once the guides are
      # collected; it is what the legacy chunk passed.
      assemble_feature_panel(
        list(
          highlighted_feature_plot("FIG7C", "muscle_transcript",
                                   color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG7C", "adipose_transcript",
                                   color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG7C", "blood_olink",
                                   color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG7C", "muscle_protein",
                                   color_time_labels = FALSE,
                                   include_legend = TRUE,
                                   legend_position = "bottom")
        ),
        ncol = 2,
        legend_position = "bottom"
      )
    }
  ),

  FIG7B = list(
    panel = "FIG7B",
    figure = "Figure 7B",
    title = "CX3CL1 (fractalkine) across muscle and adipose transcript and Olink",
    build = function() {
      assemble_feature_panel(
        list(
          highlighted_feature_plot("FIG7B", "muscle_transcript",
                                   color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG7B", "adipose_transcript",
                                   color_time_labels = FALSE,
                                   include_legend = FALSE),
          highlighted_feature_plot("FIG7B", "olink", color_time_labels = FALSE,
                                   include_legend = TRUE)
        ),
        ncol = 3,
        legend_position = "bottom"
      )
    }
  ),

  # ---- Extended Data 9 -----------------------------------------------------

  ED9C = list(
    panel = "ED9C",
    figure = "Extended Data 9C",
    title = "FLNA phosphosite S1533, muscle and adipose",
    build = function() {
      # The same residue is quantified under two accessions, muscle on the FLNA
      # isoform-2 accession and adipose on the canonical one, so "all" tissues
      # resolves to one facet per plot.
      assemble_feature_panel(
        list(
          highlighted_feature_plot("ED9C", "muscle_s1533",
                                   color_time_labels = FALSE,
                                   include_legend = TRUE),
          highlighted_feature_plot("ED9C", "adipose_s1533",
                                   color_time_labels = FALSE,
                                   include_legend = TRUE)
        ),
        ncol = 1,
        legend_position = "bottom"
      )
    }
  )
)

# ---- Reading the catalog ---------------------------------------------------

#' Build one single-feature plot by catalog name.
#'
#' This is what a panel script calls. It returns the plot and writes nothing;
#' the panel script hands the result to export_panel().
#'
#' @param name A name in SINGLE_FEATURE_PLOTS. Every entry is a panel, so the
#'   name IS the panel id.
single_feature_plot <- function(name) {
  entry <- single_feature_entry(name)
  entry$build()
}

#' The catalog entry for one name, with a message naming the alternatives when
#' there is no such entry.
single_feature_entry <- function(name) {
  if (!name %in% names(SINGLE_FEATURE_PLOTS)) {
    stop(
      "unknown single-feature plot '", name, "'. Known plots: ",
      paste(names(SINGLE_FEATURE_PLOTS), collapse = ", "),
      call. = FALSE
    )
  }
  SINGLE_FEATURE_PLOTS[[name]]
}

#' The catalog names that carry a panel id, in catalog order.
single_feature_panel_names <- function() {
  Filter(function(name) !is.na(SINGLE_FEATURE_PLOTS[[name]]$panel),
         names(SINGLE_FEATURE_PLOTS))
}

# ---- Running the whole catalog ---------------------------------------------

if (sys.nframe() == 0) {
  # Attached, not merely loaded: plot_single_feature() calls
  # load_differential_analysis(), which resolves its lazy-loaded DA objects by
  # name through the search path.
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(MotrpacHumanPreSuspensionData)
    library(MotrpacHumanPreSuspensionAnalysis)
  })

  # run_panels() is the repo's runner, so these sixteen report exactly as the
  # figure scripts' panels do. Keyed by panel id rather than catalog name: the
  # two agree for every entry that has a panel.
  entries <- single_feature_panel_names()
  panels <- lapply(entries, function(name) {
    force(name)
    function() {
      panel_init(SINGLE_FEATURE_PLOTS[[name]]$panel)
      export_panel(single_feature_plot(name), SINGLE_FEATURE_PLOTS[[name]]$panel)
    }
  })
  names(panels) <- vapply(entries, function(n) SINGLE_FEATURE_PLOTS[[n]]$panel,
                          character(1))

  run_panels(panels)
}
