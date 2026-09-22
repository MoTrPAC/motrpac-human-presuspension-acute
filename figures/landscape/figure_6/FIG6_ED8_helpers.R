# FIG6_ED8.R — what the Figure 6 and Extended Data 8 panels share, and what
# tables/ST6.R reads.
#
# One file because both figures and the table script read it:
#   exercise_contrast_matrices()  FIG6C, FIG6E, ED8A, ED8C, ED8E
#   logfc_star_heatmap()          the same five
#   tfeb_chip_evidence()          FIG6D, ED8A, ED8E
#   scion_tfeb_targets()          FIG6D, FIG6E, ED8A
#   scion_edges()                 FIG6E, ST6a
#
# Sourced after lib/panel_export.R and lib/highlights.R, which define
# panel_source() and highlight_feature(). Nothing here attaches a package:
# panel_init() attaches the Analysis package and its objects are reached by name
# through the search path.
#
# exercise_contrast_matrices() turns a DA table into logFC and adj_p matrices,
# one column per exercise-with-controls contrast, EE then RE, timepoints
# ascending; logfc_star_heatmap() draws them as the blue-white-red heatmap every
# heatmap panel of these two figures is, with a star where adj_p < 0.05.

EXERCISE_GROUP_OF <- c("EE-CON" = "ADUEndur", "RE-CON" = "ADUResist")

#' Muscle transcript DA over the exercise-with-controls contrasts, with symbols.
#'
#' Lazy and memoised here rather than in either figure script: FIG6D, FIG6E,
#' ED8A and ED8E all want the same frame with no arguments, and one run of
#' either script asks for it twice.
muscle_transcript_exercise_da <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    da <- MUSCLE_TRNSCRPT_DA[
      MUSCLE_TRNSCRPT_DA$contrast_type == "exercise_with_controls" &
        MUSCLE_TRNSCRPT_DA$contrast_category %in% names(EXERCISE_GROUP_OF), ]
    symbols <- HUMAN_FEATURE_TO_GENE[HUMAN_FEATURE_TO_GENE$assay == "transcript-rna-seq",
                                     c("feature_id", "gene_symbol")]
    symbols <- symbols[!duplicated(symbols$feature_id), ]
    da <- merge(as.data.frame(da), symbols, by = "feature_id", all.x = TRUE)
    da$gene_symbol <- as.character(da$gene_symbol)
    cache <<- da
    cache
  }
})

#' Every muscle transcript that carries a gene symbol.
muscle_transcript_symbols <- function() {
  da <- muscle_transcript_exercise_da()
  out <- unique(da[!is.na(da$gene_symbol), c("feature_id", "gene_symbol")])
  out$feature_id <- as.character(out$feature_id)
  out$gene_symbol <- as.character(out$gene_symbol)
  out[order(out$feature_id), ]
}

#' The exercise-with-controls contrasts of one DA table as two matrices.
#'
#' @param da a DA table filtered to contrast_type exercise_with_controls and
#'   contrast_category EE-CON / RE-CON. The row key is `key`.
#' @param key the column that names a row: "feature_id" or "gene_symbol".
#' @param tissue a label carried into the annotation, for a heatmap that binds
#'   several tissues side by side.
#' @return list(fc, adj_p, annotation). Columns are EE then RE, timepoints in
#'   the order names(HUMAN_ACUTE_TIMEPOINT_COLORS) gives them.
exercise_contrast_matrices <- function(da, key = "gene_symbol", tissue = NULL) {
  da <- as.data.frame(da)
  da <- da[!is.na(da[[key]]), ]
  timepoints <- names(HUMAN_ACUTE_TIMEPOINT_COLORS)
  da$group <- EXERCISE_GROUP_OF[as.character(da$contrast_category)]
  da <- da[order(match(da$group, EXERCISE_GROUP_OF),
                 match(da$Timepoint, timepoints)), ]
  columns <- unique(da[, c("group", "Timepoint")])
  column_id <- paste(columns$group, columns$Timepoint)
  da$column <- paste(da$group, da$Timepoint)
  rows <- unique(da[[key]])
  cell <- function(value) {
    m <- matrix(NA_real_, nrow = length(rows), ncol = length(column_id),
                dimnames = list(rows, column_id))
    m[cbind(match(da[[key]], rows), match(da$column, column_id))] <- da[[value]]
    m
  }
  annotation <- data.frame(group = as.character(columns$group),
                           time = as.character(columns$Timepoint),
                           stringsAsFactors = FALSE)
  if (!is.null(tissue)) annotation$tissue <- tissue
  list(fc = cell("logFC"), adj_p = cell("adj_p_value"), annotation = annotation)
}

#' One transcript per gene symbol.
#'
#' A symbol carried by more than one transcript (a PAR_Y copy, a retired id)
#' cannot be one heatmap row. The PAR_Y copies go first; if a symbol is still
#' duplicated, the transcript with the smallest adjusted p over its contrasts
#' stays, and the choice is reported.
one_transcript_per_symbol <- function(da, panel) {
  da <- da[!grepl("PAR_Y", da$feature_id, fixed = TRUE), ]
  best <- stats::aggregate(adj_p_value ~ feature_id + gene_symbol, data = da,
                           FUN = min, na.action = stats::na.pass)
  best <- best[order(best$gene_symbol, best$adj_p_value, best$feature_id), ]
  dropped <- best[duplicated(best$gene_symbol), ]
  if (nrow(dropped) > 0) {
    message(sprintf(
      "        %s: %d symbol(s) carried by more than one transcript; keeping the most responsive of each (%s)",
      panel, length(unique(dropped$gene_symbol)),
      paste(unique(dropped$gene_symbol), collapse = ", ")))
  }
  keep <- best$feature_id[!duplicated(best$gene_symbol)]
  da[da$feature_id %in% keep, ]
}

#' White on a dark cell, black on a light one.
star_color <- function(hex) {
  rgb_value <- grDevices::col2rgb(hex)
  luminance <- (0.299 * rgb_value[1, ] + 0.587 * rgb_value[2, ] +
                  0.114 * rgb_value[3, ]) / 255
  ifelse(luminance < 0.5, "white", "black")
}

#' The logFC heatmap these figures draw: rows clustered, columns in contrast
#' order under a modality and timepoint annotation, a star where the adjusted
#' p clears `alpha`.
#'
#' @param matrices the list exercise_contrast_matrices() returns.
#' @param breaks three values for the blue-white-red ramp.
#' @param show_row_names FALSE for a heatmap too tall to label.
#' @param cluster_rows FALSE keeps the rows in the order given.
#' @param row_fontsize row label size in points.
logfc_star_heatmap <- function(matrices, breaks, show_row_names = TRUE,
                               cluster_rows = TRUE, row_fontsize = 8,
                               alpha = 0.05) {
  fc <- matrices$fc
  adj_p <- matrices$adj_p
  annotation <- matrices$annotation
  # Factors with levels in study order, so the legend reads EE then RE and the
  # timepoints chronologically rather than alphabetically; palettes cut to the
  # levels drawn.
  in_order <- function(values, palette, levels = names(palette)[names(palette) %in% values]) {
    list(values = factor(values, levels = levels), colors = palette[levels])
  }
  # The palette names resistance first; the columns run EE then RE.
  modality <- in_order(annotation$group, HUMAN_EXERCISE_GROUP_COLORS,
                       levels = unique(annotation$group))
  timepoint <- in_order(annotation$time, HUMAN_ACUTE_TIMEPOINT_COLORS)
  annotation_columns <- list(Modality = modality$values, Timepoint = timepoint$values)
  annotation_colors <- list(Modality = modality$colors, Timepoint = timepoint$colors)
  if (!is.null(annotation$tissue)) {
    tissue <- in_order(annotation$tissue, HUMAN_TISSUE_COLORS)
    annotation_columns <- c(list(Tissue = tissue$values), annotation_columns)
    annotation_colors <- c(list(Tissue = tissue$colors), annotation_colors)
  }
  top <- do.call(ComplexHeatmap::HeatmapAnnotation,
                 c(annotation_columns, list(col = annotation_colors)))
  cell <- grid::unit(3, "mm")
  col_fun <- circlize::colorRamp2(breaks = breaks, colors = c("blue", "white", "red"))
  ComplexHeatmap::Heatmap(
    fc,
    col = col_fun,
    na_col = "grey90",
    heatmap_legend_param = list(title = "logFC", legend_direction = "vertical",
                                legend_width = grid::unit(50, "mm")),
    width = ncol(fc) * cell, height = nrow(fc) * cell,
    top_annotation = top,
    cluster_columns = FALSE, cluster_rows = cluster_rows,
    show_column_names = FALSE, show_row_names = show_row_names,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x = x, y = y, width, height, gp = grid::gpar(col = "#555555"))
      if (!is.na(adj_p[i, j]) && adj_p[i, j] < alpha && !is.na(fc[i, j])) {
        star <- grid::textGrob("*")
        star_w <- grid::convertWidth(grid::grobWidth(star), "mm")
        star_h <- grid::convertHeight(grid::grobHeight(star), "mm")
        grid::grid.text("*", x, y - star_h * 0.5 + star_w * 0.4,
                        gp = grid::gpar(col = star_color(col_fun(fc[i, j]))))
      }
    },
    column_title_rot = 90, row_title_rot = 0,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = row_fontsize)
  )
}

#' What export_panel() is handed for a heatmap: draw() on the open device, the
#' legends merged to the right, no new page.
draw_heatmap <- function(heatmap, legend_side = "right", annotation_legend_list = NULL) {
  function() {
    ComplexHeatmap::draw(
      heatmap, merge_legends = TRUE, heatmap_legend_side = legend_side,
      annotation_legend_list = annotation_legend_list,
      legend_gap = grid::unit(10, "mm"), newpage = FALSE)
  }
}

# ---- TFEB ChIP evidence ----------------------------------------------------

#' The gene symbols each published TFEB ChIP source calls bound.
#'
#' Gambardella et al. 2020: a peak annotated Promoter (<=1kb).
#' MSigDB TFEB_TARGET_GENES: membership.
#' Settembre et al. 2013: TFEB and PPARGC1A, by ChIP-qPCR.
#' Doronzo et al. 2019: a peak annotatr places in a promoter or 1-5 kb upstream,
#' against hg19_basicgenes.
#'
#' Lazy and memoised: building the hg19_basicgenes annotation is the expensive
#' step of these two figures, and ED8A and ED8E ask for it in the same run.
tfeb_chip_sources <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    gambardella <- openxlsx::read.xlsx(
      panel_source("figure_6", "Gambardella-et-al-2020-table-s2.xlsx"), sheet = 2)
    gambardella_bound <- unique(gambardella$SYMBOL[
      grepl("Promoter", gambardella$annotation, fixed = TRUE) &
        grepl("<=1", gambardella$annotation, fixed = TRUE)])

    msigdb_lines <- readLines(panel_source("figure_6", "TFEB_TARGET_GENES.v2024.1.Hs.tsv"), warn = FALSE)
    symbol_line <- msigdb_lines[startsWith(msigdb_lines, "GENE_SYMBOLS\t")]
    if (length(symbol_line) != 1) {
      stop("TFEB_TARGET_GENES.v2024.1.Hs.tsv has no single GENE_SYMBOLS row", call. = FALSE)
    }
    msigdb_bound <- strsplit(sub("^GENE_SYMBOLS\t", "", symbol_line), ",", fixed = TRUE)[[1]]

    regions <- annotatr::read_regions(
      panel_source("figure_6", "GSM2354032_TFEB_e5_peaks.bed"), genome = "hg19")
    annotations <- annotatr::build_annotations(genome = "hg19", annotations = "hg19_basicgenes")
    annotated <- as.data.frame(annotatr::annotate_regions(
      regions, annotations, ignore.strand = TRUE, quiet = TRUE))
    doronzo_bound <- unique(annotated$annot.symbol[
      grepl("promoter", annotated$annot.type, fixed = TRUE) |
        grepl("1to5kb", annotated$annot.type, fixed = TRUE)])
    doronzo_bound <- doronzo_bound[!is.na(doronzo_bound)]

    cache <<- list(Gambardella_et_al_2020 = gambardella_bound,
                   MSigDb_TFEB_TARGET_GENES = msigdb_bound,
                   Settembre_et_al_2013 = c("TFEB", "PPARGC1A"),
                   Doronzo_et_al_2018 = doronzo_bound)
    cache
  }
})

#' Y/N per gene symbol for each ChIP source, and `any` across them.
tfeb_chip_evidence <- function(gene_symbols, sources = tfeb_chip_sources()) {
  out <- data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE)
  for (source_name in names(sources)) {
    out[[source_name]] <- ifelse(gene_symbols %in% sources[[source_name]], "Y", "N")
  }
  out$any <- ifelse(rowSums(out[, names(sources), drop = FALSE] == "Y") > 0, "Y", "N")
  message(sprintf("        TFEB ChIP evidence over %d symbols: %s; any %d",
                  nrow(out),
                  paste(sprintf("%s %d", names(sources),
                                colSums(out[, names(sources), drop = FALSE] == "Y")),
                        collapse = ", "),
                  sum(out$any == "Y")))
  out
}

# ---- the SC-ION network -----------------------------------------------------
#
# The network inference, the permutation trim at edge weight 0.1 and the
# Cytoscape merge behind these two tables are not reproduced here. Both are
# vendored under figure_6/sources/ and read as given; see
# docs/external_dependencies.md.

#' The SC-ION-predicted TFEB targets, as the PR recorded them.
scion_tfeb_targets <- function() {
  targets <- readxl::read_excel(panel_source("figure_6", "TFEB_targets_v2.xlsx"))
  as.character(targets[[1]])
}

#' The merged EE+RE network's edges, with the regulator and target parsed out
#' of the edge name.
#'
#' Not memoised: FIG6E is the only panel that reads it, and tables/ST6.R reads
#' it in a run of its own.
scion_edges <- function() {
  path <- panel_source("figure_6", "SCION_supplemental_table_edges.csv",
                       env_var = "FIG6_SCION_EDGES_CSV")
  # Comma-delimited, quoted: a Cytoscape export rather than a freeze file.
  edges <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("", "NA"))
  required <- c("name", "regulator_id", "site", "target_id", "modality",
                "tissue", "cluster", "weight")
  missing_columns <- setdiff(required, names(edges))
  if (length(missing_columns) > 0) {
    stop("the SC-ION edge table is missing column(s): ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  edges$regulator <- sub(" \\(regulates\\) .*$", "", edges$name)
  edges$target <- sub("^.* \\(regulates\\) ", "", edges$name)
  edges
}

#' The targets one regulator points at: outgoing edges only. A first-neighbour
#' definition would also return the regulator's own upstream inputs.
scion_targets_of <- function(regulator, edges) {
  unique(edges$target[edges$regulator == regulator])
}
