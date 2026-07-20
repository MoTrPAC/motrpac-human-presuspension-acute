# Helper functions for n_documentation_features.Rmd
# All output figures and files should be written to the outputs/ subdirectory.

plot_ome_stats_histogram = function(group_stats_df, title, subtitle = NULL, vline = NULL) {
  group_stats_df %>%
    ggplot(aes(x = Count, fill = assay)) +
    geom_histogram(
      binwidth = 1,
      color = "white",
      linewidth = 0.25,
      alpha = 0.85,
      position = "stack"
    ) +
    scale_fill_manual(values = HUMAN_OME_COLORS, name = "Assay") +
    { if (!is.null(vline)) geom_vline(xintercept = vline, linetype = "dashed", color = "grey30", linewidth = 0.6) } +
    scale_x_continuous(
      limits = c(0, 75),
      breaks = seq(0, 75, by = 5),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Observed samples per group × timepoint for each feature",
      y = "Number of feature × group × timepoint observations"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      plot.margin = margin(10, 16, 10, 10)
    )
}

plot_group_stats_histogram = function(group_stats_df, title, subtitle = NULL, vline = NULL, facet_vars = NULL) {
  group_labels = c(
    "ADUControl" = "Control",
    "ADUEndur"   = "Endurance",
    "ADUResist"  = "Resistance"
  )
  group_stats_df %>%
    mutate(
      Group = group_labels[randomGroupCode],
      Group = factor(Group, levels = c("Control", "Endurance", "Resistance"))
    ) %>%
    ggplot(aes(x = Count, fill = Group)) +
    geom_histogram(
      binwidth = 1,
      color = "white",
      linewidth = 0.25,
      alpha = 0.85,
      position = "identity"
    ) +
    scale_fill_manual(values = HUMAN_EXERCISE_GROUP_COLORS, name = "Group", guide = "none") +
    { if (!is.null(vline)) geom_vline(xintercept = vline, linetype = "dashed", color = "grey30", linewidth = 0.6) } +
    scale_x_continuous(
      limits = c(0, 75),
      breaks = seq(0, 75, by = 5),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Observed samples per group × timepoint for each feature",
      y = "Number of feature × group × timepoint observations"
    ) +
    { if (!is.null(facet_vars)) facet_wrap(facet_vars, scales = "free_y") } +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      plot.margin = margin(10, 16, 10, 10)
    )
}

plot_sig_nonsig_histogram = function(sig_df, nonsig_df, title, subtitle = NULL,
                                     vline = NULL, facet_vars = vars(tissue, ome)) {
  # Pre-aggregate rather than letting stat_bin do it. stat_bin emits a rectangle
  # for every bin x fill combination, including empty ones; those become
  # zero-height stroked rects in the PDF that Illustrator shows as stray lines.
  facet_cols = if (is.null(facet_vars)) character(0) else sapply(facet_vars, rlang::as_name)

  plot_data = dplyr::bind_rows(
    dplyr::mutate(sig_df, Category = "Significant"),
    dplyr::mutate(nonsig_df, Category = "Non-significant")
  ) %>%
    dplyr::mutate(
      Category = factor(Category, levels = c("Significant", "Non-significant"))
    ) %>%
    dplyr::count(dplyr::across(dplyr::all_of(c(facet_cols, "Category", "Count"))),
                 name = "n_obs") %>%
    dplyr::filter(n_obs > 0)

  x_max = max(plot_data$Count, na.rm = TRUE)

  strip_labels = c(
    "prot-ph" = "Phosphoproteomics",
    "prot-pr" = "Proteomics",
    "adipose" = "Adipose",
    "muscle" = "Muscle",
    "blood" = "Blood"
  )
  relabel_strip = function(x) {
    x = as.character(x)
    out = unname(strip_labels[x])
    ifelse(is.na(out), x, out)
  }

  ggplot(plot_data, aes(x = Count, y = n_obs, fill = Category)) +
    # geom_col with no border: a stroked bar keeps its stroke width when the
    # figure is scaled down in Illustrator, so thin bars collapse into slivers.
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(
      values = c("Significant" = "#E64B35", "Non-significant" = "#ADB5BD"),
      name = NULL
    ) +
    { if (!is.null(vline)) geom_vline(xintercept = vline, linetype = "dashed", color = "grey30", linewidth = 0.6) } +
    scale_x_continuous(
      # bars are centred on integer counts, so pad by half a bin on each side
      limits = c(-0.5, x_max + 0.5),
      breaks = scales::breaks_width(5),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Observed samples per group × timepoint for each feature",
      y = "Number of feature × group × timepoint observations"
    ) +
    { if (!is.null(facet_vars)) facet_wrap(facet_vars, scales = "free_y",
                                           labeller = labeller(.default = relabel_strip)) } +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      strip.text = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 13),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      plot.margin = margin(10, 16, 10, 10),
      legend.position = "bottom"
    )
}
