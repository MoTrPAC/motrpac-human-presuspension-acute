# FIG4_ED5.R — the CAMERA-PR bubble heatmap and the fuzzy c-means objects, as
# Figure 4 and Extended Data 5 read them.
#
# One file because both figures read both halves:
#   camera_enrich_heatmap()   FIG4A, FIG4B, ED5A
#   the cmeans_* functions    FIG4C-FIG4F, ED5B-ED5E
#
# Sourced after lib/panel_export.R, which defines landscape_root().
#
# Everything hand-set lives in figure_4/cmeans.env and is read through
# cmeans_clusters() and cmeans_setting() below. Nothing here hardcodes a cluster
# number.

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(dplyr)
})

# ============================================================================
# CAMERA-PR
# ============================================================================
#
# FIG4A, FIG4B and ED5A are all one call to
# MotrpacHumanPreSuspensionAnalysis::plot_enrich_heatmap() in the legacy, over
# CAMERA_RESULTS, differing only in ome and selection:
#
#   FIG4A  transcript-rna-seq, twelve curated sets
#   FIG4B  prot-pr, eight curated sets
#   ED5A   metab, the six most significant sets per tissue and contrast
#
# That function cannot be called here. It opens its own cairo_pdf() on the
# filename it is given and closes it with on.exit(), so it writes a PDF instead
# of returning a drawing, and export_panel() has nothing to export. Its page size
# is computed from the data and is not a value the manifest could carry either.
#
# camera_enrich_heatmap() below is that function with the device removed. The
# body is otherwise the upstream one: the same filters, the same annotation
# frame, the same column split, the same size arithmetic, and the same
# TMSig::enrichmap() call. Three changes, all forced:
#
#   - grDevices::cairo_pdf() and its on.exit(dev.off()) are gone. The caller
#     hands the returned $draw to export_panel(), which opens the device the
#     manifest declares.
#   - draw_args gains newpage = FALSE. enrichmap() ends in ComplexHeatmap::draw(),
#     whose default starts a fresh page on the device export_panel() already
#     opened, leaving page 1 blank and the heatmap on page 2. ED2E and FIG3B
#     carry the same argument for the same reason.
#   - the computed width and height are returned rather than consumed, so the
#     caller passes them to export_panel(). A panel whose row count decides its
#     page size cannot take the size from the manifest.
#
# The three unexported helpers it reaches for (.add_contrast_labels(),
# .contrast_colors(), .enrich_heatmap_color_function) are called through ::: on
# the installed package rather than copied, so the labels, the contrast colours
# and the z-score colour ramp stay whatever the pinned version says they are.
#
# Provenance: MotrpacHumanPreSuspensionAnalysis R/plot_enrich_heatmap.R, at the
# version config/required_packages.tsv pins. The legacy
# callers are figures/landscape/CAMERA-PR/curated_CAMERA-PR_heatmaps.R.

#' The CAMERA-PR bubble heatmap, as a drawing and the page it needs.
#'
#' @param set_ids Character set ids to draw, or NULL to take the n_top most
#'   significant per tissue and contrast.
#' @param selected_ome One of the omes CAMERA_RESULTS carries.
#' @param selected_tissues Tissues to include. With more than one, only sets
#'   tested in at least two of them are kept.
#' @param contrast_type Which contrast family to draw.
#' @param padj_cutoff The adjusted-p cut an asterisk marks.
#' @param n_top Sets per tissue and contrast when set_ids is NULL.
#' @param zscore_colors Length-2 ramp for the smallest and largest z-scores.
#'
#' @return list(draw, width, height). draw is zero-argument, for export_panel().
camera_enrich_heatmap <- function(set_ids = NULL,
                                  selected_ome,
                                  selected_tissues = c("adipose", "blood", "muscle"),
                                  contrast_type = "exercise_with_controls",
                                  padj_cutoff = 0.05,
                                  n_top = 6L,
                                  zscore_colors = c("#3366ff", "darkred")) {
  pkg <- asNamespace("MotrpacHumanPreSuspensionAnalysis")

  x <- MotrpacHumanPreSuspensionAnalysis::CAMERA_RESULTS

  required_cols <- c("tissue", "assay", "contrast", "contrast_type",
                     "set_id", "set", "set_short", "p_value", "adj_p_value",
                     "z.std")
  missing_cols <- setdiff(required_cols, colnames(x))
  if (length(missing_cols) > 0) {
    stop("CAMERA_RESULTS is missing column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  x <- x %>%
    dplyr::filter(.data$contrast_type == !!contrast_type,
                  .data$assay == selected_ome,
                  .data$tissue %in% selected_tissues) %>%
    droplevels.data.frame()

  selected_tissues <- intersect(selected_tissues, unique(as.character(x$tissue)))

  # With more than one tissue, a set has to have been tested in at least two of
  # them to be comparable across the columns.
  if (length(selected_tissues) > 1L) {
    x <- dplyr::filter(x, .by = "set", length(unique(.data$tissue)) > 1L)
  }
  x <- dplyr::filter(x, .by = "set", any(.data$adj_p_value < padj_cutoff))

  if (nrow(x) == 0) {
    stop("no sets are significant at padj_cutoff = ", padj_cutoff,
         " for ome ", selected_ome, call. = FALSE)
  }

  if (is.null(set_ids)) {
    # p_value rather than adj_p_value for the ordering, to avoid ties.
    set_ids <- x %>%
      dplyr::filter(.data$adj_p_value < padj_cutoff) %>%
      dplyr::slice_min(.data$p_value, by = c("tissue", "contrast"), n = n_top) %>%
      dplyr::pull("set_id") %>%
      unique() %>%
      as.character()
  } else {
    set_ids <- sprintf("%05d", as.numeric(gsub("[^[:digit:]]", "", set_ids)))
    set_ids <- unique(set_ids)
    found <- intersect(set_ids, as.character(x$set_id))
    if (length(found) == 0) {
      stop("none of the curated set ids survives this run's CAMERA-PR for ome ",
           selected_ome, call. = FALSE)
    }
    if (length(found) < length(set_ids)) {
      message(sprintf(
        "        %d of %d curated set(s) not returned by this run and not drawn: %s",
        length(set_ids) - length(found), length(set_ids),
        paste(sort(setdiff(set_ids, found)), collapse = ", ")))
    }
  }

  x <- dplyr::filter(x, .data$set_id %in% set_ids)

  contrast_df <- pkg$.add_contrast_labels() %>%
    dplyr::filter(.data$contrast %in% levels(x$contrast)) %>%
    droplevels.data.frame()

  column_df <- dplyr::distinct(x, .data$tissue, .data$contrast) %>%
    dplyr::mutate(tissue = factor(.data$tissue, levels = selected_tissues)) %>%
    dplyr::arrange(.data$contrast, .data$tissue) %>%
    dplyr::left_join(contrast_df, by = "contrast") %>%
    dplyr::rename(modality = "anno_group") %>%
    dplyr::mutate(modality = factor(.data$modality, levels = c("EE", "RE")))

  anno_df <- dplyr::select(column_df,
                           Tissue = "tissue",
                           Modality = "modality",
                           Timepoint = "contrast_labels")

  anno_col <- list(
    Tissue = MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[selected_tissues],
    Modality = stats::setNames(c("#d95f02", "#1b9e77"), c("EE", "RE")),
    Timepoint = pkg$.contrast_colors()[levels(anno_df$Timepoint)]
  )

  # A contrast in one tissue is a different column from the same contrast in
  # another.
  column_df <- column_df %>%
    dplyr::mutate(contrast2 = paste(.data$tissue, .data$contrast),
                  contrast2 = factor(.data$contrast2, levels = unique(.data$contrast2)))

  if (length(selected_tissues) == 1L) {
    anno_df$Tissue <- NULL
    anno_col["Tissue"] <- NULL
  }

  show_column_names <- FALSE
  if (identical(contrast_type, "baseline")) {
    anno_df$Timepoint <- NULL
    anno_col["Timepoint"] <- NULL
    show_column_names <- TRUE
  }
  if (!contrast_type %in% c("exercise_with_controls", "exercise_no_controls")) {
    anno_df$Modality <- NULL
    anno_col["Modality"] <- NULL
  }

  top_annotation <- NULL
  column_split <- NULL

  if (length(anno_col) > 0) {
    if (!is.null(anno_df[["Tissue"]]) && !is.null(anno_df[["Modality"]])) {
      column_split <- anno_df %>%
        dplyr::mutate(column_split = paste(.data$Tissue, .data$Modality),
                      row_order = seq_len(dplyr::n())) %>%
        dplyr::arrange(.data$Tissue, .data$Modality) %>%
        dplyr::mutate(column_split = factor(.data$column_split,
                                            levels = unique(.data$column_split))) %>%
        dplyr::arrange(.data$row_order) %>%
        dplyr::pull("column_split")
    } else if (!is.null(anno_df[["Tissue"]])) {
      column_split <- anno_df[["Tissue"]]
    } else if (!is.null(anno_df[["Modality"]])) {
      column_split <- anno_df[["Modality"]]
    }

    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      df = anno_df,
      col = anno_col,
      which = "column",
      border = TRUE,
      gap = grid::unit(2, "pt"),
      annotation_name_gp = grid::gpar(fontsize = 0.9 * 14),
      annotation_legend_param = list(
        border = TRUE,
        title_gp = grid::gpar(fontsize = 0.9 * grid::unit(14, "pt"),
                              fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 0.9 * grid::unit(14, "pt"))
      )
    )
  }

  # One 14 pt cell per set, plus room for the annotation and the legends.
  n_sets <- length(unique(x[["set"]]))
  height <- grid::convertUnit(n_sets * grid::unit(14, "pt"), "in")
  height <- max(as.numeric(height), 4 + 2 * (1L - show_column_names)) +
    3 + show_column_names

  row_label_width <- ComplexHeatmap::max_text_width(
    text = x[["set_short"]],
    gp = grid::gpar(fontsize = 0.9 * 14)
  )
  row_label_width <- grid::convertUnit(row_label_width, "in")
  width <- grid::convertUnit(nrow(anno_df) * grid::unit(14, "pt"), "in")
  width <- as.numeric(width + row_label_width) + 3

  extended_range <- TMSig::extendRangeNum(x[["z.std"]], nearest = 0.1)
  breaks <- if (all(extended_range <= 0)) {
    c(extended_range[1], 0)
  } else if (all(extended_range >= 0)) {
    c(0, extended_range[2])
  } else {
    c(extended_range[1], 0, extended_range[2])
  }

  # Row clustering fails when a set is measured in too few contrasts to cluster
  # on, so it is turned off for the same case upstream turns it off for.
  cluster_rows <- x %>%
    dplyr::count(.data$set_short) %>%
    dplyr::pull("n") %>%
    {all(. >= 0.5 * max(.))}

  x <- x %>%
    dplyr::mutate(contrast2 = paste(.data$tissue, .data$contrast),
                  contrast2 = factor(.data$contrast2,
                                     levels = levels(column_df$contrast2)),
                  set_short = factor(.data$set_short,
                                     levels = sort(unique(.data$set_short))))

  # Few rows leave the legends taller than the heatmap; the padding keeps the
  # page from growing to fit them.
  draw_args <- list(newpage = FALSE)
  if (length(unique(x$set)) <= 10L) {
    draw_args$padding <- grid::unit(c(185, 0, 0, 0), "pt")
  } else if (length(unique(x$set)) <= 20L) {
    draw_args$padding <- grid::unit(c(70, 0, 0, 0), "pt")
  }

  draw <- function() {
    TMSig::enrichmap(
      x = as.data.frame(x),
      n_top = Inf,
      set_column = "set_short",
      statistic_column = "z.std",
      contrast_column = "contrast2",
      padj_column = "adj_p_value",
      padj_cutoff = padj_cutoff,
      plot_sig_only = TRUE,
      heatmap_color_fun = pkg$.enrich_heatmap_color_function,
      colors = zscore_colors,
      padj_legend_title = "BH Adjusted\nP-Value",
      draw_args = draw_args,
      heatmap_args = list(
        na_col = "grey95",
        rect_gp = grid::gpar(fill = "white", col = "grey85"),
        cluster_rows = cluster_rows,
        column_split = column_split,
        row_labels = latex2exp::TeX(levels(x$set_short)),
        column_labels = column_df$contrast_labels,
        show_column_names = show_column_names,
        column_names_side = "top",
        column_title_gp = grid::gpar(fontsize = 0),
        top_annotation = top_annotation,
        heatmap_legend_param = list(
          title = "Z-Score",
          at = breaks,
          labels = breaks
        )
      )
    )
  }

  list(draw = draw, width = width, height = height, n_sets = n_sets)
}

# ============================================================================
# Fuzzy c-means
# ============================================================================
#
# Eight panels are drawn from FCM_CLUSTERS and FCM_ORA: FIG4C's and FIG4E's
# sparkline grids, FIG4D's cluster bar-and-dot summary and ED5B's full grid of
# it, the three cross-tissue heatmaps (FIG4F, ED5D, ED5E) and ED5C's UpSet. The
# legacy Rmd built the same three frames separately in each chunk, which is how
# its cluster numbers came to disagree between chunks after the re-clustering.
#
# Everything below is namespace-qualified, so a figure script that sources this
# file inherits no attached packages from the c-means half.

# ---- the hand-set numbers --------------------------------------------------

# figure_4/cmeans.env is read on load, so a panel gets the cluster numbers by
# sourcing this helper and nothing else. It is a shell file because it is a list
# of settings and reads like one; parsed here rather than run through a shell so
# an R session needs no subprocess.
#
# Values already exported in the environment WIN, which is what makes a one-off
# override possible:
#   CMEANS_VEGF_BLOOD=4,5 Rscript figures/landscape/ED5.R ED5C
.cmeans_load_env <- function() {
  path <- file.path(landscape_root(), "figure_4", "cmeans.env")
  if (!file.exists(path)) {
    stop("cmeans.env is missing from config/", call. = FALSE)
  }
  for (line in readLines(path, warn = FALSE)) {
    line <- trimws(line)
    if (!startsWith(line, "export CMEANS_")) next
    body <- sub("^export ", "", line)
    name <- sub("=.*$", "", body)
    if (nzchar(Sys.getenv(name, unset = ""))) next          # already set: leave it
    value <- sub("^[^=]*=", "", body)
    value <- sub('^"', "", sub('"$', "", value))
    # "${CMEANS_X:-default}" -> default
    value <- sub("^\\$\\{[A-Za-z0-9_]+:-", "", value)
    value <- sub("\\}$", "", value)
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(NULL)
}

.cmeans_load_env()


#' A cluster list from figure_4/cmeans.env.
#'
#' @param name The suffix after CMEANS_, e.g. "VEGF_BLOOD".
#' @param numeric Return integers rather than the character labels the FCM
#'   objects key on.
#' @returns The parsed comma-separated list, in the order written.
cmeans_clusters <- function(name, numeric = FALSE) {
  var <- paste0("CMEANS_", name)
  raw <- Sys.getenv(var, unset = "")
  if (!nzchar(raw)) {
    stop(var, " is not set. It is defaulted in figure_4/cmeans.env, which this ",
         "file reads on load.", call. = FALSE)
  }
  parts <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    stop(var, " is set but empty: '", raw, "'", call. = FALSE)
  }
  if (numeric) {
    out <- suppressWarnings(as.integer(parts))
    if (anyNA(out)) {
      stop(var, " is not a comma-separated list of integers: '", raw, "'",
           call. = FALSE)
    }
    return(out)
  }
  parts
}

source(file.path(landscape_root(), "lib", "highlights.R"))

#' The signature the c-means panels are built around, and its display label,
#' from config/highlights.json.
cmeans_pathway <- function() highlight_pathways("cmeans")$set
cmeans_pathway_label <- function() highlight_pathways("cmeans")$label

#' A single scalar from figure_4/cmeans.env.
cmeans_setting <- function(name) {
  var <- paste0("CMEANS_", name)
  raw <- Sys.getenv(var, unset = "")
  if (!nzchar(raw)) {
    stop(var, " is not set; see figure_4/cmeans.env.", call. = FALSE)
  }
  raw
}

# ---- tissues ---------------------------------------------------------------

# The order the panels list tissues in, and the columns of the heatmap strip.
CMEANS_TISSUES <- c("adipose", "blood", "muscle")

# ---- the timepoint axis ----------------------------------------------------

# The abbreviations the c-means panels label with, and the order they run in.
# One definition: the legacy Rmd wrote this map out twice, at lines 69 and 249,
# and the two had already diverged in which timepoints they carried.
CMEANS_TIMEPOINTS <- c(
  "pre_exercise"      = "Pre",
  "post_10_min"       = "P10M",
  "post_15_30_45_min" = "P15-45M",
  "post_3.5_4_hr"     = "P3.5/4H",
  "post_24_hr"        = "P24H"
)

#' Blood's post-10-min and post-15/30/45-min groups, merged.
#'
#' A blood trajectory peaking at either is the same response, and keeping them
#' apart split one cluster's header label across two columns. Muscle and adipose
#' are untouched. Controlled by CMEANS_BLOOD_MERGE_P10_INTO_P15.
#'
#' @param tissue,timepoint Character vectors of the same length.
cmeans_merge_blood_timepoints <- function(tissue, timepoint) {
  if (!identical(toupper(Sys.getenv("CMEANS_BLOOD_MERGE_P10_INTO_P15", "TRUE")),
                 "TRUE")) {
    return(timepoint)
  }
  merged <- timepoint == "P10M" & tissue == "blood"
  timepoint[merged] <- "P15-45M"
  timepoint
}

# ---- the cluster objects ---------------------------------------------------

#' Cluster centroids, tidied: one row per tissue x cluster x arm x timepoint.
#'
#' FCM_CLUSTERS is a list of `fclust` objects, one per tissue, and its centroid
#' matrix is keyed by the full contrast expression:
#'
#'   "group_timepointADUEndur.post_24_hr - group_timepointADUEndur.pre_exercise
#'    - group_timepointADUControl.post_24_hr + group_timepointADUControl.pre_exercise"
#'
#' Every panel needs the arm and the timepoint out of that, and the legacy
#' parsed it inline in each chunk. Parsed once here.
#'
#' The tissues do NOT share a timepoint axis: blood carries post_10_min where
#' muscle and adipose start at post_15_30_45_min, and adipose is clustered at
#' k = 13 against 12 for the other two. A panel that assumes a common grid is
#' wrong, so the frame is returned long and each panel states its own order.
#'
#' @returns A data.frame: tissue, cluster (character), arm ("EE"/"RE"),
#'   timepoint, timepoint_short, value.
cmeans_centroids <- function() {
  if (!"package:MotrpacHumanPreSuspensionAnalysis" %in% search()) {
    stop("cmeans_centroids() needs MotrpacHumanPreSuspensionAnalysis on the ",
         "search path: declare it in the panel's data_packages.", call. = FALSE)
  }
  clusters <- MotrpacHumanPreSuspensionAnalysis::FCM_CLUSTERS

  rows <- lapply(names(clusters), function(tissue) {
    centers <- clusters[[tissue]]$centers
    if (is.null(centers)) {
      stop("FCM_CLUSTERS$", tissue, " has no centers matrix", call. = FALSE)
    }
    # The left-hand term of the contrast is the arm and timepoint being tested.
    lhs <- sub(" - .*$", "", colnames(centers))
    arm <- sub("^group_timepointADU([A-Za-z]+)\\..*$", "\\1", lhs)
    timepoint <- sub("^group_timepointADU[A-Za-z]+\\.", "", lhs)

    known_arm <- c("Endur" = "EE", "Resist" = "RE")
    if (!all(arm %in% names(known_arm))) {
      stop(tissue, ": centroid column with an unrecognised arm: ",
           paste(unique(arm[!arm %in% names(known_arm)]), collapse = ", "),
           call. = FALSE)
    }
    unknown_tp <- setdiff(timepoint, names(CMEANS_TIMEPOINTS))
    if (length(unknown_tp) > 0) {
      stop(tissue, ": centroid column with an unrecognised timepoint: ",
           paste(sort(unique(unknown_tp)), collapse = ", "),
           "\n  Add it to CMEANS_TIMEPOINTS.", call. = FALSE)
    }

    observed <- data.frame(
      tissue = tissue,
      cluster = rep(rownames(centers), times = ncol(centers)),
      arm = rep(unname(known_arm[arm]), each = nrow(centers)),
      timepoint = rep(timepoint, each = nrow(centers)),
      value = as.vector(centers),
      stringsAsFactors = FALSE
    )

    # Every centroid is a delta-delta against pre-exercise, so pre-exercise is
    # zero by construction and carries no column in the matrix. It is still a
    # point on the trajectory, and the one every cluster starts from: without it
    # a sparkline begins at its first post-exercise value and a rise off
    # baseline is indistinguishable from a fall onto it.
    #
    # The legacy cbind()s a zero column per arm for exactly this reason. Adding
    # it also gives muscle a fourth x position, which is what lets all four of
    # its selected clusters take a distinct label in FIG4E.
    baseline <- expand.grid(
      cluster = rownames(centers),
      arm = unname(known_arm[unique(arm)]),
      stringsAsFactors = FALSE
    )
    baseline$tissue <- tissue
    baseline$timepoint <- "pre_exercise"
    baseline$value <- 0

    rbind(observed, baseline[, names(observed)])
  })

  out <- do.call(rbind, rows)
  out$timepoint_short <- unname(CMEANS_TIMEPOINTS[out$timepoint])
  # NOT merged here. The blood P10M/P15-45M merge belongs to the GROUPING that
  # picks a cluster's palette, not to the x axis: collapsing two timepoints into
  # one axis position gives a cluster two values at the same x, which doubles
  # back the sparkline and draws its label twice. The legacy applied the merge
  # to max_abs_contrast alone, and so does cmeans_cluster_groups().
  out
}

#' Hard cluster assignment per feature, with the assay split back out.
#'
#' fclust names each feature "<assay> <feature_id>"; the panels join on
#' feature_id and colour by assay, so the two are separated here rather than in
#' each caller.
#'
#' @returns A data.frame: tissue, assay, feature_id, cluster (character).
cmeans_assignments <- function() {
  clusters <- MotrpacHumanPreSuspensionAnalysis::FCM_CLUSTERS
  rows <- lapply(names(clusters), function(tissue) {
    assignment <- clusters[[tissue]]$cluster
    labels <- names(assignment)
    if (is.null(labels)) {
      stop("FCM_CLUSTERS$", tissue, "$cluster has no feature names",
           call. = FALSE)
    }
    data.frame(
      tissue = tissue,
      assay = sub(" .*$", "", labels),
      feature_id = sub("^[^ ]+ ", "", labels),
      cluster = as.character(unname(assignment)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' The over-representation results behind the cluster annotations.
cmeans_ora <- function() {
  if (!"package:MotrpacHumanPreSuspensionAnalysis" %in% search()) {
    stop("cmeans_ora() needs MotrpacHumanPreSuspensionAnalysis on the search ",
         "path: declare it in the panel's data_packages.", call. = FALSE)
  }
  MotrpacHumanPreSuspensionAnalysis::FCM_ORA
}

#' The features of one molecular signature.
#'
#' The legacy Rmd carried this as `.choose_pathway_features_text()`, vendored
#' into the script after the Analysis package withdrew its version. It is here
#' rather than in a panel because three panels call it.
#'
#' @param pathway A name in MOLECULAR_SIGNATURES, e.g. "WP_VEGFA_VEGFR2_SIGNALING".
#' @returns The signature's members, deduplicated.
cmeans_pathway_features <- function(pathway) {
  signatures <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES
  hit <- unlist(lapply(signatures, function(collection) collection[[pathway]]),
                use.names = FALSE)
  if (length(hit) == 0) {
    stop("no molecular signature named '", pathway, "'.",
         "\n  The cmeans selection in config/highlights.json names it.", call. = FALSE)
  }
  sort(unique(hit))
}

# ---- palettes --------------------------------------------------------------

#' Tissue colours for the tissues actually present, keyed by name.
cmeans_tissue_colors <- function(tissues) {
  hits <- MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[tolower(tissues)]
  if (anyNA(hits)) {
    stop("tissue not in HUMAN_TISSUE_COLORS: ",
         paste(tissues[is.na(hits)], collapse = ", "), call. = FALSE)
  }
  stats::setNames(unname(hits), tissues)
}

# ---- how clusters are grouped and coloured ---------------------------------

# The pathway row colours of FIG4D and ED5B: the legacy's ten-colour list,
# recycled down the rows in plotting order.
CMEANS_SET_PALETTE <- c(
  "#CC0000", "#006400", "#E65C00", "#0000B2", "#B8860B",
  "#008080", "#8B008B", "#556B2F", "#800000", "#4B0082"
)

# Each cluster is grouped by WHERE its trajectory peaks and in which direction,
# and every group gets one RColorBrewer ramp so clusters that behave alike are
# shades of one hue. The mapping is editorial and is stated here rather than
# derived.
CMEANS_GROUP_PALETTES <- c(
  "P3.5/4H_Negative" = "Reds",
  "P3.5/4H_Positive" = "Greens",
  "P15-45M_Positive" = "Blues",
  "P24H_Positive"    = "Purples",
  "P15-45M_Negative" = "Oranges",
  "P24H_Negative"    = "YlOrBr"
)

#' Label each cluster by its peak timepoint and direction, and colour it.
#'
#' Adds peak_timepoint, peak_sign, peak_group and cluster_color to the centroid
#' frame. One definition, so every panel groups and colours a cluster the same
#' way and the blood timepoint merge cannot differ between them.
#'
#' The peak is taken over the largest absolute centroid across both arms, so a
#' cluster is grouped by its strongest response whichever arm produced it.
#'
#' @param centroids Output of cmeans_centroids().
cmeans_cluster_groups <- function(centroids) {
  peak <- do.call(rbind, lapply(split(centroids,
                                      list(centroids$tissue, centroids$cluster),
                                      drop = TRUE), function(rows) {
    at <- which.max(abs(rows$value))
    data.frame(
      tissue = rows$tissue[1],
      cluster = rows$cluster[1],
      peak_timepoint = rows$timepoint_short[at],
      peak_sign = if (rows$value[at] >= 0) "Positive" else "Negative",
      stringsAsFactors = FALSE
    )
  }))
  # Here is where blood's P10M folds into P15-45M: a blood trajectory peaking at
  # either is the same response, and separating them split one cluster's palette
  # group across two headers.
  peak$peak_timepoint <- cmeans_merge_blood_timepoints(peak$tissue,
                                                       peak$peak_timepoint)
  peak$peak_group <- paste(peak$peak_timepoint, peak$peak_sign, sep = "_")

  unmapped <- setdiff(unique(peak$peak_group), names(CMEANS_GROUP_PALETTES))
  if (length(unmapped) > 0) {
    stop("cluster peak group with no palette: ",
         paste(sort(unmapped), collapse = ", "),
         "\n  Add it to CMEANS_GROUP_PALETTES. The legacy warned and fell back ",
         "to grey, which silently made two groups look like one.", call. = FALSE)
  }

  # One ramp per tissue x group, clusters ordered so the shade is stable.
  peak$cluster_color <- NA_character_
  for (key in unique(paste(peak$tissue, peak$peak_group))) {
    rows <- which(paste(peak$tissue, peak$peak_group) == key)
    rows <- rows[order(as.integer(peak$cluster[rows]))]
    palette_name <- CMEANS_GROUP_PALETTES[[peak$peak_group[rows[1]]]]
    n <- length(rows)
    shades <- if (n > 9) {
      grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(9, palette_name)))(n)
    } else {
      rev(RColorBrewer::brewer.pal(max(3, n), palette_name))[seq_len(n)]
    }
    peak$cluster_color[rows] <- shades
  }

  merge(centroids, peak, by = c("tissue", "cluster"), sort = FALSE)
}

#' The frame every c-means panel plots from.
#'
#' cmeans_centroids() tidied, grouped and coloured, with the timepoint axis
#' ordered. One call replaces the ~120 lines each legacy chunk opened with.
#'
#' Memoised, and memoised HERE rather than in a figure script: FIG4C, FIG4D and
#' FIG4E plot from it, ED5B and ED5C read their colours off it, and
#' cmeans_cluster_row_annotation() calls it once per heatmap the three
#' cross-tissue panels draw. It takes no arguments, so every caller wants the
#' same frame and one FCM_CLUSTERS pass serves a whole figure.
cmeans_plot_frame <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    frame <- cmeans_cluster_groups(cmeans_centroids())
    present <- intersect(unname(CMEANS_TIMEPOINTS), unique(frame$timepoint_short))
    frame$timepoint_short <- factor(frame$timepoint_short, levels = present)
    frame$arm <- factor(frame$arm, levels = c("EE", "RE"))
    cache <<- frame[order(frame$tissue, as.integer(frame$cluster), frame$arm,
                          frame$timepoint_short), ]
    cache
  }
})

# ---- hard assignment and enrichment ----------------------------------------

#' Features assigned to a cluster above the membership threshold.
#'
#' Fuzzy c-means gives every feature a membership in every cluster. The panels
#' that count features per cluster need a hard assignment, which is the cluster
#' with the highest membership PROVIDED that membership clears
#' CMEANS_CLUSTER_THRESHOLD; below it the feature is assigned to none.
#'
#' This is not the same as fclust's own `$cluster`, which assigns every feature
#' to its maximum regardless of how weak that maximum is. cmeans_assignments()
#' returns that; this returns the thresholded version FCM_ORA was built on.
#'
#' Not memoised: it takes a threshold, and a memoised accessor keyed on nothing
#' would hand a caller that passed one the frame built for another.
#'
#' @param threshold Membership a feature must exceed. Defaults to the config.
#' @returns A data.frame: tissue, assay, feature_id, cluster (character, NA
#'   below threshold), membership.
cmeans_hard_assignments <- function(threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- as.numeric(cmeans_setting("CLUSTER_THRESHOLD"))
  }
  if (!is.finite(threshold)) {
    stop("CMEANS_CLUSTER_THRESHOLD is not a number", call. = FALSE)
  }
  clusters <- MotrpacHumanPreSuspensionAnalysis::FCM_CLUSTERS

  rows <- lapply(names(clusters), function(tissue) {
    membership <- clusters[[tissue]]$membership
    if (is.null(membership)) {
      stop("FCM_CLUSTERS$", tissue, " has no membership matrix", call. = FALSE)
    }
    top <- max.col(membership, ties.method = "first")
    best <- membership[cbind(seq_len(nrow(membership)), top)]
    labels <- rownames(membership)

    data.frame(
      tissue = tissue,
      assay = sub(" .*$", "", labels),
      feature_id = sub("^[^ ]+ ", "", labels),
      cluster = ifelse(best > threshold, colnames(membership)[top], NA_character_),
      membership = best,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' The significant cluster enrichments.
#'
#' FCM_ORA at FDR < 0.05, which is what "this cluster is enriched for that
#' pathway" means everywhere in these panels.
cmeans_ora_hits <- function(fdr = 0.05) {
  ora <- cmeans_ora()
  hits <- ora[!is.na(ora$adj_p_value) & ora$adj_p_value < fdr, , drop = FALSE]
  hits$set <- as.character(hits$set)
  hits$cluster <- as.character(hits$cluster)
  hits
}

#' Which tissue x cluster combinations are enriched for one signature.
#'
#' @param pathway A signature name or set_id.
#' @param ome The assay the enrichment was computed on.
#' @returns A data.frame: tissue, cluster, group_id ("<tissue>-<cluster>").
cmeans_enriched_clusters <- function(pathway, ome) {
  hits <- cmeans_ora_hits()
  keep <- (hits$set == pathway | as.character(hits$set_id) == pathway) &
    as.character(hits$assay) == ome
  found <- hits[keep, c("tissue", "cluster"), drop = FALSE]
  if (nrow(found) == 0) {
    stop("no cluster is enriched for '", pathway, "' on ", ome,
         "\n  config/highlights.json (cmeans) and CMEANS_OME in figure_4/cmeans.env name them.",
         call. = FALSE)
  }
  # as.numeric first: FCM_ORA stores the cluster as "01" in places and "1" in
  # others, and the group id has to match the centroid frame's rownames.
  found$cluster <- as.character(as.numeric(found$cluster))
  found <- unique(found)
  found$group_id <- paste0(found$tissue, "-", found$cluster)
  found
}

# ---- label placement -------------------------------------------------------

#' One timepoint per cluster, spread along the x axis.
#'
#' Greedy: rank every (cluster, timepoint) pair by how far that cluster travels
#' from zero at that timepoint, take the highest, then retire both that cluster
#' and that timepoint and repeat. A cluster is labelled at its own peak unless a
#' stronger cluster claimed that timepoint first.
#'
#' @param rows Centroid rows for one panel.
#' @param fallback When TRUE, a cluster left over once the timepoints run out is
#'   labelled at its own peak anyway, sharing a timepoint. FIG4C does this so
#'   every cluster is named; FIG4E does not, and reports what it dropped.
#' @param score Column ranked on. FIG4C ranks on distance from zero in either
#'   direction, FIG4E on the signed height.
cmeans_label_positions <- function(rows, fallback = TRUE, score = c("abs", "signed")) {
  score <- match.arg(score)
  rows$.score <- if (score == "abs") abs(rows$value) else rows$value
  ranked <- stats::aggregate(.score ~ cluster + timepoint_short, data = rows,
                             FUN = max)
  # as.character on both: these arrive as factors, and c(character(0), <factor>)
  # yields the level INDEX, which makes the retire list match nothing.
  ranked$cluster <- as.character(ranked$cluster)
  ranked$timepoint_short <- as.character(ranked$timepoint_short)
  ranked <- ranked[order(-ranked$.score), , drop = FALSE]

  all_clusters <- unique(ranked$cluster)
  taken_cluster <- character(0)
  taken_timepoint <- character(0)
  chosen <- list()
  remaining <- ranked
  while (nrow(remaining) > 0) {
    winner <- remaining[1, , drop = FALSE]
    chosen[[length(chosen) + 1L]] <- winner[, c("cluster", "timepoint_short")]
    taken_cluster <- c(taken_cluster, winner$cluster)
    taken_timepoint <- c(taken_timepoint, winner$timepoint_short)
    remaining <- remaining[!remaining$cluster %in% taken_cluster &
                             !remaining$timepoint_short %in% taken_timepoint, ,
                           drop = FALSE]
  }
  placed <- do.call(rbind, chosen)

  if (fallback) {
    left <- setdiff(all_clusters, placed$cluster)
    if (length(left) > 0) {
      extra <- do.call(rbind, lapply(left, function(cl) {
        best <- ranked[ranked$cluster == cl, , drop = FALSE][1, , drop = FALSE]
        best[, c("cluster", "timepoint_short")]
      }))
      placed <- rbind(placed, extra)
    }
  }
  placed
}

#' Apply the hand-set label moves from figure_4/cmeans.env.
#'
#' Both entries exist because the greedy result still collides at this size, and
#' both move the label only. Recorded as configuration rather than as an `if`
#' inside the panel, which is where the legacy kept them.
cmeans_apply_label_overrides <- function(placed, tissue) {
  raw <- Sys.getenv("CMEANS_LABEL_OVERRIDES", unset = "")
  if (!nzchar(raw)) return(placed)

  for (rule in trimws(strsplit(raw, ",", fixed = TRUE)[[1]])) {
    parts <- strsplit(rule, ":", fixed = TRUE)[[1]]
    if (length(parts) < 3 || parts[1] != tissue) next
    cluster <- parts[2]
    if (!cluster %in% placed$cluster) next

    target <- if (parts[3] == "match") {
      other <- placed$timepoint_short[placed$cluster == parts[4]]
      if (length(other) > 0) other[1] else parts[5]
    } else {
      parts[3]
    }
    placed$timepoint_short[placed$cluster == cluster] <- target
  }
  placed
}

# ---- the cross-tissue feature heatmaps -------------------------------------

#' Feature ids of a signature, restricted to one ome.
cmeans_pathway_feature_ids <- function(pathway, ome) {
  map <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
  key <- if (ome == "metab") "refmet_name" else "gene_symbol"
  members <- cmeans_pathway_features(pathway)
  ids <- unique(as.character(
    map$feature_id[as.character(map[[key]]) %in% members &
                     as.character(map$assay) == ome]
  ))
  if (length(ids) == 0) {
    stop("'", pathway, "' maps to no ", ome, " feature", call. = FALSE)
  }
  ids
}

#' The features one heatmap draws, and the cluster each falls in per tissue.
#'
#' A heatmap selects the pathway features that land in named clusters of a
#' PRIMARY tissue, optionally narrowed to those that also land in named clusters
#' of a SECOND tissue. The second is what makes the panel cross-tissue: it says
#' these features move together in muscle and in adipose, not merely that they
#' move in muscle.
#'
#' @param pathway,ome The signature and the assay it was enriched on.
#' @param tissue,clusters The primary tissue and its clusters.
#' @param with_tissue,with_clusters Optional second tissue to intersect with.
#' @returns A list: `feature_ids`, and `annotation` — one row per feature with a
#'   column per tissue (all three) carrying the cluster that tissue assigns it,
#'   NA where the feature is unassigned or unmeasured there.
cmeans_heatmap_features <- function(pathway, ome, tissue, clusters,
                                    with_tissue = NULL, with_clusters = NULL) {
  pathway_ids <- cmeans_pathway_feature_ids(pathway, ome)

  assignments <- cmeans_hard_assignments()
  assignments <- assignments[assignments$assay == ome &
                               assignments$feature_id %in% pathway_ids, ,
                             drop = FALSE]

  in_primary <- assignments$feature_id[assignments$tissue == tissue &
                                         !is.na(assignments$cluster) &
                                         assignments$cluster %in% clusters]
  keep <- unique(in_primary)

  if (!is.null(with_tissue)) {
    in_second <- assignments$feature_id[assignments$tissue == with_tissue &
                                          !is.na(assignments$cluster) &
                                          assignments$cluster %in% with_clusters]
    keep <- intersect(keep, unique(in_second))
  }

  if (length(keep) == 0) {
    stop("no ", ome, " feature of '", pathway, "' falls in ", tissue,
         " cluster(s) ", paste(clusters, collapse = ", "),
         if (!is.null(with_tissue)) {
           paste0(" and ", with_tissue, " cluster(s) ",
                  paste(with_clusters, collapse = ", "))
         } else "",
         call. = FALSE)
  }

  # One column per tissue, whichever the panel selected on: the legacy strip
  # names every tissue's cluster for every feature.
  annotation <- data.frame(feature_id = sort(keep), stringsAsFactors = FALSE)
  for (tis in CMEANS_TISSUES) {
    rows <- assignments[assignments$tissue == tis, , drop = FALSE]
    annotation[[tis]] <- rows$cluster[match(annotation$feature_id,
                                            rows$feature_id)]
  }

  # plot_feature_heatmap() relabels every row to its gene symbol where it has
  # one (falling back to the feature id), so the annotation has to carry that
  # same label or it cannot be matched to the drawn rows at all.
  map <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE
  map <- map[as.character(map$assay) == ome, , drop = FALSE]
  symbol <- as.character(map$gene_symbol)[match(annotation$feature_id,
                                                as.character(map$feature_id))]
  annotation$label <- ifelse(is.na(symbol), annotation$feature_id, symbol)

  clashes <- unique(annotation$label[duplicated(annotation$label)])
  if (length(clashes) > 0) {
    stop("two features share a row label, so the cluster annotation cannot be ",
         "matched to the heatmap rows: ", paste(clashes, collapse = ", "),
         call. = FALSE)
  }

  list(feature_ids = annotation$feature_id, annotation = annotation)
}

#' Row annotation showing which cluster each feature falls in, per tissue.
#'
#' Three strips, adipose, blood, muscle, coloured with the same per-cluster
#' palette FIG4C and FIG4D use, so a block of one colour down this annotation is
#' the same cluster the sparklines name.
#'
#' `feature_order` must be the row order of the matrix the annotation is joined
#' to, NOT the order the features were selected in: ComplexHeatmap pairs an
#' annotation to heatmap rows positionally, so a mismatch silently labels every
#' row with another row's cluster.
cmeans_cluster_row_annotation <- function(annotation, feature_order) {
  # Matched on `label`, not feature_id: the heatmap draws gene symbols.
  missing <- setdiff(feature_order, annotation$label)
  if (length(missing) > 0) {
    stop("no cluster assignment for heatmap row(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  frame <- annotation[match(feature_order, annotation$label), , drop = FALSE]

  # Every tissue, a tissue selected on one cluster included: that strip is one
  # colour, and the legacy draws it.
  tissues <- intersect(CMEANS_TISSUES, names(frame))
  palette <- unique(cmeans_plot_frame()[, c("tissue", "cluster", "cluster_color")])
  colours <- lapply(tissues, function(tis) {
    per_tissue <- palette[palette$tissue == tis, , drop = FALSE]
    stats::setNames(per_tissue$cluster_color, as.character(per_tissue$cluster))
  })
  names(colours) <- tissues

  strip <- as.data.frame(lapply(frame[tissues], as.character),
                         stringsAsFactors = FALSE)
  names(strip) <- tissues

  # ComplexHeatmap defaults for size, NA colour and legend text, as the legacy
  # rowAnnotation(df = , col = ) call left them.
  ComplexHeatmap::rowAnnotation(df = strip, col = colours)
}

#' plot_feature_heatmap() with a cluster annotation down the right-hand side,
#' rows restricted to the features measured in every tissue and clustered.
#'
#' Draws on the open device, so call it inside export_panel().
cmeans_feature_heatmap <- function(annotation, ...) {
  heatmap <- MotrpacHumanPreSuspensionAnalysis::plot_feature_heatmap(
    ...,
    multi_tissue_clust_rows = TRUE,
    right_annotation = function(row_labels) {
      cmeans_cluster_row_annotation(annotation, row_labels)
    },
    return_drawing = TRUE
  )
  heatmap$draw()
  invisible(NULL)
}
