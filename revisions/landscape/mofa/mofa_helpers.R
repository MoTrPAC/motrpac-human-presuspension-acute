HUMAN_TP_LEVELS = c("Pre", "D20M", "D40M", "P10M", "P15-45M", "P3.5/4H", "P24H")

recode_timepoint = function(tp) {
  dplyr::case_when(
    tp == "pre_exercise"       ~ "Pre",
    tp == "during_20_min"      ~ "D20M",
    tp == "during_40_min"      ~ "D40M",
    tp == "post_10_min"        ~ "P10M",
    tp == "post_15_30_45_min"  ~ "P15-45M",
    tp == "post_3.5_4_hr"      ~ "P3.5/4H",
    tp == "post_24_hr"         ~ "P24H",
    TRUE ~ tp
  )
}

recode_modality = function(rg) {
  dplyr::case_when(
    rg == "ADUEndur"   ~ "EE",
    rg == "ADUResist"  ~ "RE",
    rg == "ADUControl" ~ "CON",
    TRUE ~ rg
  )
}


#' Plot individual feature weightings for selected factors, faceted by view
#'
#' Extracts the top \code{n_features} features (by absolute weight) from each
#' view for each requested factor and draws a lollipop/dot plot.  Positive
#' weights point right and negative weights point left so the sign is
#' immediately readable.
#'
#' @param mofa_trained Trained MOFA2 object.
#' @param factors Integer vector of factors to display (default 1:3).
#' @param views Character vector of view names to include.  NULL (default) uses
#'   all views.
#' @param n_features Number of top features per view per factor (by |weight|).
#' @param view_colors Named character vector mapping view names to colors.  NULL
#'   uses a built-in palette.
#' @param png_path Path for PNG export at 600 dpi.  NA (default) skips export.
#'
#' @return A ggplot2 object (invisibly when \code{png_path} is set).
plot_feature_weights = function(mofa_trained, factors = 1:3, views = NULL,
                                n_features = 10, view_colors = NULL,
                                group = NULL, png_path = NA) {
  all_views = MOFA2::views_names(mofa_trained)
  if (is.null(views)) views = all_views
  views = intersect(views, all_views)
  if (length(views) == 0) stop("No matching views found in the MOFA object.")

  weights_list = MOFA2::get_weights(mofa_trained, views = views, as.data.frame = TRUE)
  # get_weights returns a data.frame with columns: feature, factor, view, value
  factor_labels = paste0("Factor", factors)
  df = weights_list %>%
    filter(factor %in% factor_labels, view %in% views) %>%
    group_by(factor, view) %>%
    arrange(desc(abs(value))) %>%
    slice_head(n = n_features) %>%
    ungroup() %>%
    mutate(
      feature = sub(paste0("_(",  paste(toupper(views), collapse="|"), ")$"), "",
                    feature, ignore.case = TRUE),
      factor_num = as.integer(sub("Factor", "", factor)),
      factor_label = paste0("Factor ", factor_num),
      factor_label = factor(factor_label,
                            levels = paste0("Factor ", sort(unique(factor_num)))),
      view = factor(view, levels = views)
    )

  # Order features within each facet by weight for readability
  df = df %>%
    arrange(factor_label, view, value) %>%
    mutate(feature_order = paste0(factor_label, "__", view, "__", feature),
           feature_order = factor(feature_order, levels = unique(feature_order)))

  if (is.null(view_colors)) {
    default_pal = HUMAN_OME_COLORS
    view_colors = default_pal[intersect(names(default_pal), views)]
    missing = setdiff(views, names(view_colors))
    if (length(missing) > 0) {
      extra = scales::hue_pal()(length(missing))
      names(extra) = missing
      view_colors = c(view_colors, extra)
    }
  }

  p = ggplot2::ggplot(df, ggplot2::aes(x = value, y = feature_order, color = view)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = value, y = feature_order, yend = feature_order),
      linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "gray40",
                        linewidth = 0.3) +
    ggplot2::facet_grid(view ~ factor_label, scales = "free_y", space = "free_y") +
    ggplot2::scale_color_manual(values = view_colors, guide = "none") +
    ggplot2::scale_y_discrete(labels = function(x) {
      lbl <- sub(".*__", "", x)
      ifelse(nchar(lbl) > 18, paste0(substr(lbl, 1, 18), "..."), lbl)
    }) +
    ggplot2::labs(
      title = if (!is.null(group)) paste0(group, " — Top feature weights per factor") else "Top feature weights per factor",
      x = "Weight",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(face = "bold", size = 18),
      strip.text.y = ggplot2::element_text(face = "bold", angle = 0, size = 18),
      axis.text.y  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 17),
      plot.title   = ggplot2::element_text(size = 19, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (!is.na(png_path)) {
    n_cols = length(factors)
    n_rows = length(views)
    ggplot2::ggsave(png_path, plot = p, dpi = 600,
                    width = 10, height = 1.5 + 0.2 * n_features * n_rows,
                    units = "in")
    cat("Saved:", png_path, "\n")
    invisible(p)
  } else {
    p
  }
}
