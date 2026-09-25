# single_feature_helpers.R — the layout the single-feature plots are assembled
# with.
#
# In lib/ rather than in a figure script because its one caller,
# single_feature_plots.R, spans ten figures. A helper used by one figure lives
# in that figure's script; this one has no figure to belong to.
#
# Most trajectories are drawn by
# MotrpacHumanPreSuspensionAnalysis::plot_single_feature(), which the catalog
# calls directly; ED3B and ED3C use sex_differences_single_feature(). What is
# left here is the assembly: several of those plots laid out in a grid under one
# collected legend.
#
# Nothing here attaches a package; every call is namespace-qualified.

# plot_single_feature() scales its text, line and legend sizes by
# scale_factor * 0.7. The assembled legend below is spaced in the same units so
# that a sub-plot and the shared legend stay in proportion.
SINGLE_FEATURE_SCALE <- 0.7

#' Lay trajectory plots out in a grid under one shared legend.
#'
#' @param plots list of ggplots, in sub-panel order (filled row by row).
#' @param ncol number of columns in the grid.
#' @param legend_position where the collected legend is drawn.
assemble_feature_panel <- function(plots,
                                   ncol,
                                   legend_position = c("right", "bottom")) {
  legend_position <- match.arg(legend_position)
  sc <- SINGLE_FEATURE_SCALE

  patchwork::wrap_plots(plots) +
    patchwork::plot_layout(ncol = ncol, guides = "collect") &
    ggplot2::theme(
      legend.position = legend_position,
      legend.spacing.x = grid::unit(0.01 * sc, "in"),
      legend.spacing.y = grid::unit(0.1 * sc, "in"),
      legend.box.spacing = grid::unit(0.01 * sc, "in"),
      legend.key.size = grid::unit(0.1 * sc, "in")
    )
}
