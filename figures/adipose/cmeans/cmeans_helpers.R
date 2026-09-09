# Helpers for adipose_cmeans.Rmd: the trajectory panel, the cluster enrichment
# bubble heatmap, and the combined Fig 4C/4D-style layout of the two.
#
# The bubble heatmap uses the same grammar as plot_cluster_enrichment() /
# TMSig::enrichmap(): fill is -log10(p), radius is -log10(BH p) scaled to the
# row maximum, and a black outline marks BH p < padj_cutoff.

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(grid)
library(tibble)
library(tidyr)

# Marker set names -> the labels used in the manuscript figure.
CELLTYPE_LABELS = c(Adip_1 = "Adipocyte-1",
                    Adip_2 = "Adipocyte-2",
                    LAM = "LAMΦ",
                    Mast = "Mast cell",
                    NK.T = "NK/T cell",
                    Pre_Ad = "Pre-adipocyte",
                    Resident = "Resident MΦ",
                    Stem = "Stem cell",
                    Vascular = "Vascular cell")

MODALITY_LABELS = c(Endur = "EE", Resist = "RE")

# SRON iridescent scheme, the palette plot_cmeans() uses.
# https://personal.sron.nl/~pault/
CMEANS_PALETTE = c("#fefbe9", "#fcf7d5", "#f5f3c1", "#eaf0b5", "#ddecbf",
                   "#d0e7ca", "#c2e3d2", "#b5ddd8", "#a8d8dc", "#9bd2e1",
                   "#8dcbe4", "#81c4e7", "#7bbce7", "#7eb2e4", "#88a5dd",
                   "#9398d2", "#9b8ac4", "#9d7db2", "#9a709e", "#906388",
                   "#805770", "#684957", "#46353a")

# One faceted ggplot rather than plot_cmeans()' grid of standalone plots.
# facet_wrap() is what the manuscript figure uses: it draws the timepoint axis
# only under panels with nothing below them, keeps every panel rectangle the
# same size with no space reserved on the rows that have no axis, and carries a
# single legend and y-axis title. The data steps below mirror plot_cmeans().
#
# The legend is placed inside the empty grid slot left by the last row.
# One modality at a time; plot_cmeans() itself covers the two-modality case.
plot_cmeans_facets = function(FCM,
                              modality = c("Endur", "Resist"),
                              keep_clusters,
                              min_membership = 0,
                              common_ylim = TRUE,
                              ncol = 4L,
                              legend_position = c(0.68, 0)) {
  if (!inherits(FCM, "fclust"))
    stop("`FCM` must be an object of class 'fclust'. See plot_cmeans().")

  modality = match.arg(modality)
  min_membership = max(0, min(min_membership, 1, na.rm = TRUE))

  if (missing(keep_clusters))
    keep_clusters = as.integer(rownames(FCM[["centers"]]))

  mem = FCM[["membership"]][, keep_clusters, drop = FALSE]
  centers = FCM[["centers"]][keep_clusters, , drop = FALSE]
  colnames(mem) = rownames(centers) = paste0("Cluster ", colnames(mem))
  cluster_levels = rownames(centers)

  strip_contrast = function(nms) {
    paste0("timepoint_", sub("group_timepointADU([^ ]+) .*", "\\1", nms))
  }
  colnames(centers) = strip_contrast(colnames(centers))
  zmat = FCM[["input"]]
  colnames(zmat) = strip_contrast(colnames(zmat))

  centers = centers[, grepl(modality, colnames(centers)), drop = FALSE]
  zmat = zmat[, grepl(modality, colnames(zmat)), drop = FALSE]

  # Pre-exercise is the zero every contrast is measured against, so it is added
  # rather than read out of the matrices.
  add_pre_exercise = function(m) {
    pre = matrix(0, nrow(m), 1L,
                 dimnames = list(rownames(m),
                                 sprintf("timepoint_%s.pre_exercise", modality)))
    cbind(pre, m)
  }
  zmat = add_pre_exercise(zmat)
  centers = add_pre_exercise(centers)

  tidy_timepoint = function(d) {
    d %>%
      mutate(timepoint = sub("^[^.]+\\.(.*)$", "\\1", timepoint),
             timepoint = gsub("(?<=\\d)_(?=\\d)", "/", timepoint, perl = TRUE),
             timepoint = gsub("_", " ", timepoint))
  }

  z_df = zmat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(cols = -feature, names_to = "timepoint",
                        names_pattern = "timepoint_(.*)", values_to = "z") %>%
    tidy_timepoint()
  timepoint_levels = unique(z_df$timepoint)

  cluster_assignment = data.frame(feature = names(FCM[["cluster"]]),
                                  cluster = paste0("Cluster ", FCM[["cluster"]]))

  cluster_df = mem %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    tidyr::pivot_longer(cols = -feature, names_to = "cluster",
                        values_to = "membership") %>%
    filter(membership >= min_membership) %>%
    inner_join(cluster_assignment, by = c("feature", "cluster")) %>%
    mutate(cluster = factor(cluster, levels = cluster_levels))

  # Ordering by membership makes the confident features draw last, on top.
  df = inner_join(z_df, cluster_df, by = "feature") %>%
    arrange(cluster, membership) %>%
    mutate(feature = factor(feature, levels = unique(feature)),
           timepoint = factor(timepoint, levels = timepoint_levels)) %>%
    droplevels.data.frame() %>%
    mutate(cluster = factor(cluster, levels = cluster_levels))

  center_df = centers %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cluster") %>%
    tidyr::pivot_longer(cols = -cluster, names_to = "timepoint",
                        names_pattern = "timepoint_(.*)", values_to = "center") %>%
    tidy_timepoint() %>%
    mutate(timepoint = factor(timepoint, levels = timepoint_levels),
           cluster = factor(cluster, levels = cluster_levels)) %>%
    filter(cluster %in% unique(df$cluster))

  ylim = if (common_ylim) range(df$z) else c(NA_real_, NA_real_)

  ggplot(df, aes(x = timepoint, y = z)) +
    geom_line(aes(color = membership, group = feature)) +
    geom_line(aes(x = timepoint, y = center, group = cluster),
              data = center_df, inherit.aes = FALSE, color = "black") +
    facet_wrap(~ cluster, ncol = ncol) +
    scale_x_discrete(name = NULL, expand = expansion(add = 0.2)) +
    scale_y_continuous(name = "Scaled Z-Score", limits = ylim) +
    scale_color_gradientn(name = "Membership Probability",
                          colors = CMEANS_PALETTE,
                          values = seq(0, 1, length.out = length(CMEANS_PALETTE)),
                          limits = c(0, 1),
                          breaks = seq(0, 1, 0.2),
                          guide = guide_colourbar(title.position = "top",
                                                  barwidth = unit(1.7, "in"),
                                                  barheight = unit(0.18, "in"))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing.x = unit(10, "pt"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = rel(1.1),
                                    margin = margin(t = 3)),
          axis.title.y = element_text(color = "black", face = "bold",
                                      margin = margin(r = 8)),
          # Upright, not angled: at 45 degrees a label's horizontal reach is
          # ~0.8in against ~0.3in of tick spacing, so the labels overrun each
          # other inside their own panel. Only shorter timepoint names would
          # make an angle fit at this panel width.
          axis.text.x = element_text(size = rel(0.9), color = "black",
                                     angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(color = "black"),
          legend.position = "inside",
          legend.position.inside = legend_position,
          legend.direction = "horizontal",
          legend.title = element_text(hjust = 0, size = rel(1.1)),
          legend.text = element_text(size = rel(1)),
          # theme_bw() fills the legend background white, which erases the
          # timepoint labels the inset legend sits under.
          legend.background = element_rect(fill = NA, colour = NA),
          legend.key = element_rect(fill = NA, colour = NA))
}

# `x` needs block, set_label, cluster, p_value and adj_p_value; block becomes
# the row slice. Returns the heatmap and its significance legend so the caller
# can draw them on a page of its own choosing.
build_cluster_bubble_heatmap = function(x, column_title = NULL,
                                        padj_cutoff = 0.05,
                                        cell_size = unit(19, "pt")) {
  fontsize = 0.9 * as.numeric(convertUnit(cell_size, "pt"))
  # Legends are set smaller than the heatmap labels: much above this the TeX
  # subscript in the colour bar title descends onto the first tick.
  legend_fontsize = 13.5

  x = x %>%
    mutate(cluster = as.integer(as.character(cluster))) %>%
    filter(.by = c(block, set_label), any(adj_p_value < padj_cutoff))

  rows = x %>%
    distinct(block, set_label) %>%
    arrange(block, set_label)
  if (anyDuplicated(rows$set_label))
    stop("`set_label` must be unique across blocks.")

  clusters = sort(unique(x$cluster))
  stat_mat = matrix(NA_real_, nrow(rows), length(clusters),
                    dimnames = list(rows$set_label, clusters))
  padj_mat = stat_mat
  idx = cbind(match(x$set_label, rows$set_label), match(x$cluster, clusters))
  stat_mat[idx] = -log10(x$p_value)
  padj_mat[idx] = x$adj_p_value

  # Radius carries the adjusted p-value, scaled within each row so that every
  # set's strongest cluster is drawn at full size.
  radius_mat = -log10(padj_mat)
  radius_mat = sweep(radius_mat, 1L, apply(radius_mat, 1L, max, na.rm = TRUE), FUN = "/")
  radius_mat = radius_mat * (0.95 - 0.2) + 0.2

  col_fun = colorRamp2(c(0, max(stat_mat, na.rm = TRUE)), c("white", "#543483"))

  layer_fun = function(j, i, xx, yy, w, h, fill) {
    sig = pindex(padj_mat, i, j) < padj_cutoff
    grid.rect(x = xx, y = yy, width = w, height = h,
              gp = gpar(col = "grey85",
                        fill = ifelse(is.na(sig), "grey95",
                                      ifelse(sig, "grey", "white"))))
    grid.circle(x = xx, y = yy,
                r = pindex(radius_mat, i, j) / 2 * cell_size,
                gp = gpar(col = ifelse(is.na(sig) | !sig, NA, "black"),
                          fill = fill))
  }

  ht = Heatmap(
    matrix = stat_mat,
    col = col_fun,
    name = "log_p",
    na_col = "grey95",
    rect_gp = gpar(type = "none"),
    layer_fun = layer_fun,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    row_split = factor(rows$block, levels = unique(rows$block)),
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title_rot = 90,
    row_gap = unit(6, "pt"),
    border = TRUE,
    column_names_side = "top",
    column_names_rot = 90,
    column_names_centered = TRUE,
    column_title = column_title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    height = cell_size * nrow(stat_mat),
    width = cell_size * ncol(stat_mat),
    heatmap_legend_param = list(
      title = latex2exp::TeX("$\\bf{$-log$_{10}($P-Value$)}$"),
      title_gp = gpar(fontsize = legend_fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = legend_fontsize),
      border = "black",
      legend_height = unit(40, "mm"),
      grid_width = unit(legend_fontsize, "pt")
    )
  )

  padj_legend = Legend(
    at = 1:2,
    title = "BH Adjusted\nP-Value",
    labels = paste(c("<", "≥"), padj_cutoff),
    title_gp = gpar(fontsize = legend_fontsize, fontface = "bold"),
    labels_gp = gpar(fontsize = legend_fontsize),
    legend_gp = gpar(fill = c("grey", "white")),
    grid_height = unit(legend_fontsize, "pt"),
    grid_width = unit(legend_fontsize, "pt"),
    border = "black",
    nrow = 2L
  )

  # Measured rather than a fixed allowance, so the legend column stops leaving
  # slack when cell_size changes.
  in_width = function(x, ...) as.numeric(convertUnit(max_text_width(x, ...), "in"))
  legend_width = in_width(c("-log10(P-Value)", "BH Adjusted"),
                          gp = gpar(fontsize = legend_fontsize,
                                    fontface = "bold")) + 0.7

  height = as.numeric(convertUnit(nrow(stat_mat) * cell_size, "in")) + 1.2
  width = as.numeric(convertUnit(ncol(stat_mat) * cell_size, "in")) +
    in_width(rownames(stat_mat), gp = gpar(fontsize = fontsize)) +
    legend_width + 0.4

  list(ht = ht, padj_legend = padj_legend, height = height, width = width)
}

draw_bubble_heatmap = function(parts, newpage = TRUE) {
  draw(parts$ht,
       annotation_legend_list = list(parts$padj_legend),
       merge_legends = TRUE,
       align_heatmap_legend = "heatmap_top",
       align_annotation_legend = "heatmap_top",
       newpage = newpage)
}

# Trajectory grid on the left, enrichment heatmap on the right, under one title
# spanning both, as in Fig 4C/4D.
save_combined_panel = function(traj_plot, heatmap_parts, filename, title,
                               traj_width = 6.5, traj_height = 6.8,
                               title_height = 0.5) {
  width = traj_width + heatmap_parts$width
  body_height = max(traj_height, heatmap_parts$height)
  height = body_height + title_height

  # cairo_pdf() so the ">=" in the legend renders on Mac.
  cairo_pdf(filename = filename, height = height, width = width, onefile = FALSE)
  grid.newpage()

  pushViewport(viewport(y = 1, height = unit(title_height, "in"), just = "top"))
  grid.text(title, gp = gpar(fontsize = 16, fontface = "bold"))
  popViewport()

  pushViewport(viewport(y = 0, height = unit(body_height, "in"),
                        just = "bottom"))
  # Both halves are top-aligned and sized to their own content, so the cluster
  # panels are not stretched to fill whatever height the heatmap needs.
  pushViewport(viewport(x = 0, y = 1,
                        width = traj_width / width,
                        height = unit(min(traj_height, body_height), "in"),
                        just = c("left", "top")))
  print(traj_plot, newpage = FALSE)
  popViewport()
  # Top-aligned rather than centred, so the heatmap starts level with the first
  # row of cluster trajectories.
  pushViewport(viewport(x = 1, y = 1,
                        width = heatmap_parts$width / width,
                        height = unit(min(heatmap_parts$height,
                                          body_height), "in"),
                        just = c("right", "top")))
  draw_bubble_heatmap(heatmap_parts, newpage = FALSE)
  popViewport()
  popViewport()

  invisible(dev.off())
}
