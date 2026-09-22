#!/usr/bin/env Rscript
# Supplementary Figure 1 — skeletal muscle ATAC-seq QC
#
# Panels:  SF1A  nuclei isolation batch, isolation-to-tagmentation lag and
#                sequencing batch, per participant and timepoint
#          SF1B  fraction of reads in peaks per timepoint and exercise group
#          SF1C  correlation of PC1-PC5 with participant, sample-prep and QC
#                covariates, unadjusted and M0-adjusted
#          SF1D  M0-adjusted PC1 against nuclei isolation batch
#          SF1E  M0-adjusted PC2 by exercise group, frozen against fresh nuclei
#          SF1F  M0-adjusted PC1 against PC2, 95% ellipse per group and timepoint
#          SF1G  significant peaks per exercise arm and timepoint under the two
#                contrast types, and their overlap
#
# Seven panels over three objects: load_qc()'s muscle ATAC sample sheet, the
# TMM-normalised voom log-CPM of the raw peak counts, and the PCA of that matrix
# unadjusted and after the M0 covariates are regressed out with the exercise
# terms protected. They are built once per run by the accessors below.
#
# Needs consortium data access. The muscle epigen-atac-seq QC object is not in
# MotrpacHumanPreSuspensionData; load_qc(epigen = TRUE) downloads it into
# EPIGEN_QC_DIR on first use and reuses it after that.
#
# SF1C-SF1F additionally need the raw peak counts, which no data package ships.
# They are resolved by listing QUANT_BUCKET, so those four panels need gsutil on
# PATH, read access to the consortium bucket, and QUANT_BUCKET set — it is
# declared in config/landscape.env, which lib/panel_export.R loads on source,
# and the read below stops rather than guessing a default if it is unset.
#
#   Rscript figures/landscape/SF1.R          every panel
#   Rscript figures/landscape/SF1.R SF1B     one panel

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))

# ---- the ATAC QC objects ---------------------------------------------------
#
# Read by all seven panels of this figure and by nothing else, so they live here
# rather than under helpers/.
#
# The qc-norm matrix, sample sheet and DA table come through the data packages'
# loaders from the epigenomics staging prefix. The raw counts are the one file
# those loaders do not cover; they are resolved by listing the bucket, never by
# composing a versioned name.
#
# Everything here is namespace-qualified. epigen_qc_dir() is lib/panel_export.R's.

ATAC_TISSUE <- "muscle"
ATAC_TISSUE_CODE <- "t06-muscle"
ATAC_OME <- "epigen-atac-seq"
ATAC_SIGNIFICANCE <- 0.05

# The covariates of the M0 model, and the terms left in the data when it is
# regressed out.
ATAC_M0_COVARIATES <- c("FRiP", "BMI", "Sex", "calculatedAge", "codedsiteid",
                        "Timepoint", "randomGroupCode")
ATAC_PROTECT_EXERCISE <- c("randomGroupCode", "Timepoint")
ATAC_MODEL_LABELS <- c(no_adj = "Voom Counts", M0 = "~M0")

# Timepoint shapes, in the muscle timepoint order.
ATAC_TIMEPOINT_SHAPES <- c(15, 17, 18, 16)

# SF1D fits its trend lines to nuclei isolation batches below this.
ATAC_TREND_MAX_BATCH <- 50

# The exercise arms and contrast types SF1G counts peaks for.
ATAC_EXERCISE_GROUPS <- c("ADUEndur", "ADUResist")
ATAC_CONTRAST_TYPES <- c("exercise_with_controls", "exercise_no_controls")

# The QC metrics SF1C correlates with the PCs; one representative per cluster of
# correlated pipeline metrics.
ATAC_QC_METRICS <- c(
  "FRiP", "align.samstat.total_reads", "align.dup.pct_duplicate_reads",
  "align.nodup_samstat.total_reads", "align_enrich.tss_enrich.tss_enrich",
  "align.frag_len_stat.frac_reads_in_nfr", "peak_stat.peak_region_size.idr_opt.mean",
  "align.samstat.pct_properly_paired_reads", "peak_enrich.frac_reads_in_annot.fri_blacklist",
  "replication.num_peaks.num_peaks", "align.dup.paired_optical_duplicate_reads",
  "align.frac_mito.frac_mito_reads", "peak_stat.peak_region_size.mean",
  "peak_stat.peak_region_size.overlap_opt.mean", "peak_stat.peak_region_size.max_size"
)

# ---- where the files are ---------------------------------------------------

atac_env_or_stop <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) stop(name, " is not set — it is defined in config/landscape.env",
                       call. = FALSE)
  v
}

#' Exactly one bucket object matching a glob.
atac_bucket_file <- function(pattern) {
  gsutil <- Sys.getenv("GSUTIL", unset = "gsutil")
  hits <- suppressWarnings(
    system2(gsutil, c("ls", shQuote(pattern)), stdout = TRUE, stderr = FALSE)
  )
  hits <- hits[nzchar(hits)]
  if (length(hits) == 0) {
    stop("no bucket object matches ", pattern, call. = FALSE)
  }
  if (length(hits) > 1) {
    stop("several bucket objects match ", pattern, ": ",
         paste(basename(hits), collapse = ", "), " — pin one", call. = FALSE)
  }
  hits[[1]]
}

atac_read_bucket_file <- function(path) {
  MotrpacBicQC::dl_read_gcp(path, sep = "\t", tmpdir = epigen_qc_dir(),
                            gsutil_path = Sys.getenv("GSUTIL", unset = "gsutil"),
                            check_first = TRUE)
}

#' The gs:// path of the raw peak counts.
atac_raw_counts_path <- function() {
  atac_bucket_file(file.path(
    atac_env_or_stop("QUANT_BUCKET"), "*", "epigenomics", ATAC_TISSUE_CODE, ATAC_OME,
    paste0("*", ATAC_OME, "_counts_*.txt.gz")
  ))
}

# ---- building the shared objects -------------------------------------------

#' The muscle ATAC qc-norm matrix and sample sheet, from the staging prefix.
atac_qc_objects <- function() {
  qc <- MotrpacHumanPreSuspensionData::load_qc(
    selected_tissues = ATAC_TISSUE,
    selected_omes = ATAC_OME,
    epigen = TRUE,
    repo_local_dir = epigen_qc_dir(),
    verbose = FALSE
  )
  objects <- qc[[ATAC_TISSUE]][[ATAC_OME]]
  if (is.null(objects) || is.null(objects$qc_norm) || is.null(objects$sample_metadata)) {
    stop("load_qc(epigen = TRUE) returned no ", ATAC_TISSUE, " ", ATAC_OME,
         " objects — check gsutil access to the staging prefix", call. = FALSE)
  }
  objects
}

#' The muscle ATAC acute differential analysis, with the contrast columns
#' load_differential_analysis() merges in from CONTRAST_CONVERTER.
atac_differential_analysis <- function() {
  da <- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
    selected_omes = ATAC_OME,
    selected_tissues = ATAC_TISSUE,
    epigen = TRUE,
    repo_local_dir = epigen_qc_dir(),
    verbose = FALSE
  )
  table <- da[[ATAC_TISSUE]][[ATAC_OME]]
  if (is.null(table)) {
    stop("load_differential_analysis(epigen = TRUE) returned no ", ATAC_TISSUE, " ",
         ATAC_OME, " table — check gsutil access to the staging prefix", call. = FALSE)
  }
  as.data.frame(table)
}

#' One 0/1 column per level of each factor column, named <column>___<level>.
atac_indicator_columns <- function(sheet, columns) {
  indicators <- lapply(columns, function(column) {
    values <- droplevels(as.factor(sheet[[column]]))
    design <- stats::model.matrix(~ 0 + values)
    colnames(design) <- paste0(column, "___", levels(values))
    design
  })
  as.data.frame(do.call(cbind, indicators), check.names = FALSE)
}

#' load_qc()'s muscle ATAC sample sheet with the sample-prep covariates added,
#' one row per library, ordered by group and nuclei isolation batch.
#'
#' @param qc The list atac_qc_objects() returns.
atac_sample_metadata <- function(qc = atac_qc_objects()) {
  sheet <- qc$sample_metadata
  sheet$vialLabel <- as.character(sheet$vialLabel)
  stopifnot(setequal(sheet$vialLabel, colnames(qc$qc_norm)))

  for (column in c("codedsiteid", "Seq_flowcell_lane", "Seq_flowcell_run", "randomGroupCode")) {
    sheet[[column]] <- factor(sheet[[column]])
  }
  sheet$Seq_batch_short <- factor(paste0(sheet$Seq_flowcell_lane, "__", sheet$Seq_flowcell_run))
  sheet$FRiP <- sheet$peak_enrich.frac_reads_in_peaks.macs2.frip
  for (column in c("Nuclei_extr_date", "Tagmentation_date", "PCR_date")) {
    sheet[[column]] <- as.Date(sheet[[column]], format = "%m/%d/%y")
  }
  sheet$nucleiIsolation2Tagmentation <- as.numeric(sheet$Tagmentation_date - sheet$Nuclei_extr_date)
  sheet$Tagmentation2PCR <- as.numeric(sheet$PCR_date - sheet$Tagmentation_date)
  sheet$nuclei2atac_lag <- sheet$nucleiIsolation2Tagmentation > 0
  sheet$nuclei_isolation_batch <- sheet$Sample_batch
  sheet$library <- sheet$vialLabel
  sheet <- droplevels(sheet)

  # Participants ordered by their earliest isolation batch within group.
  earliest_batch <- stats::ave(sheet$nuclei_isolation_batch,
                               sheet$randomGroupCode, sheet$pid, FUN = min)
  sheet <- sheet[order(sheet$randomGroupCode, earliest_batch, sheet$pid, sheet$Timepoint), , drop = FALSE]
  rownames(sheet) <- NULL

  indicators <- cbind(
    data.frame(Sex = as.numeric(sheet$Sex), nuclei2atac_lag = as.numeric(sheet$nuclei2atac_lag)),
    atac_indicator_columns(sheet, c("randomGroupCode", "Timepoint", "codedsiteid", "Seq_batch_short"))
  )
  colnames(indicators) <- paste0("binarized___", colnames(indicators))
  cbind(sheet, indicators)
}

#' voom log-CPM, samples in sample-sheet order.
#'
#' @param sample_sheet From atac_sample_metadata().
#' @param qc The list atac_qc_objects() returns; its qc-norm rows decide which
#'   peaks are kept.
atac_log_cpm <- function(sample_sheet, qc = atac_qc_objects()) {
  peaks_to_keep <- gsub("-", ":", rownames(qc$qc_norm), fixed = TRUE)

  counts <- atac_read_bucket_file(atac_raw_counts_path())
  counts$feature_id <- paste(counts$chrom, counts$start, counts$end, sep = ":")
  counts <- counts[counts$feature_id %in% peaks_to_keep, , drop = FALSE]
  count_matrix <- as.matrix(counts[, sample_sheet$vialLabel, drop = FALSE])
  rownames(count_matrix) <- counts$feature_id

  log_cpm <- limma::voom(count_matrix)$E
  stopifnot(identical(colnames(log_cpm), sample_sheet$vialLabel))
  log_cpm
}

#' Regress every column of covariate_df that is not protected out of a features
#' x samples matrix. The protected terms form the design, so their effects stay.
atac_adjust_covariates <- function(log_counts, covariate_df, protect) {
  stopifnot(identical(colnames(log_counts), rownames(covariate_df)))
  covariate_df <- droplevels(covariate_df)
  protect <- intersect(protect, colnames(covariate_df))
  adjust <- setdiff(colnames(covariate_df), protect)
  if (length(adjust) == 0) return(log_counts)
  design <- if (length(protect) > 0) {
    stats::model.matrix(stats::reformulate(protect), data = covariate_df)
  } else {
    matrix(1, nrow = nrow(covariate_df), ncol = 1)
  }
  covariates <- stats::model.matrix(stats::reformulate(adjust), data = covariate_df)[, -1, drop = FALSE]
  limma::removeBatchEffect(log_counts, covariates = covariates, design = design)
}

#' PCA of the log-CPM after covariate adjustment.
#'
#' @param covariates Sample-sheet columns to put in the model; character(0) for
#'   no adjustment. Columns constant across the samples are dropped.
#' @param protect The subset of covariates whose effects are kept.
#' @return scores (the PCs joined to the sample sheet), variance (axis labels
#'   keyed by PC), in_model and protect.
atac_pca <- function(log_cpm, sample_sheet, covariates, protect) {
  covariate_df <- sample_sheet[, intersect(covariates, colnames(sample_sheet)), drop = FALSE]
  rownames(covariate_df) <- sample_sheet$library
  is_constant <- vapply(covariate_df, function(column) length(unique(column)) == 1, logical(1))
  covariate_df <- droplevels(covariate_df[, !is_constant, drop = FALSE])
  protect <- intersect(protect, colnames(covariate_df))

  adjusted <- atac_adjust_covariates(log_cpm[, rownames(covariate_df), drop = FALSE],
                                     covariate_df, protect)
  pca <- stats::prcomp(scale(t(adjusted)), rank. = 10)

  pc_names <- colnames(pca$x)
  explained <- scales::percent(pca$sdev^2 / sum(pca$sdev^2), 0.01)[seq_along(pc_names)]
  variance <- stats::setNames(paste0(pc_names, " (", explained, ")"), pc_names)

  scores <- data.frame(library = rownames(pca$x), pca$x, check.names = FALSE,
                       stringsAsFactors = FALSE)
  scores <- merge(scores, sample_sheet, by = "library", sort = FALSE)

  list(scores = scores, variance = variance,
       in_model = colnames(covariate_df), protect = protect)
}

#' The two all-sample PCAs the panels read: unadjusted, and M0-adjusted with the
#' exercise terms protected.
atac_pca_models <- function(log_cpm, sample_sheet) {
  list(
    no_adj = atac_pca(log_cpm, sample_sheet, character(0), ATAC_PROTECT_EXERCISE),
    M0 = atac_pca(log_cpm, sample_sheet, ATAC_M0_COVARIATES, ATAC_PROTECT_EXERCISE)
  )
}

# ---- SF1C's covariate rows -------------------------------------------------

#' The sample-sheet columns SF1C correlates with the PCs: category, the column
#' holding the value, and the label it is drawn with.
atac_covariate_columns <- function(sample_sheet) {
  selected <- rbind(
    data.frame(category = "Biological",
               col_name = c("randomGroupCode", "Timepoint", "Sex", "calculatedAge", "BMI")),
    data.frame(category = "Experimental",
               col_name = c("codedsiteid", "Seq_batch_short", "nuclei_isolation_batch", "nuclei2atac_lag")),
    data.frame(category = "QC Metrics",
               col_name = intersect(ATAC_QC_METRICS, colnames(sample_sheet)))
  )

  binarized_names <- grep("^binarized___", colnames(sample_sheet), value = TRUE)
  parts <- strsplit(sub("^binarized___", "", binarized_names), "___", fixed = TRUE)
  binarized <- data.frame(
    full_col_name = binarized_names,
    col_name = vapply(parts, `[`, character(1), 1),
    short_name = vapply(parts, function(x) if (length(x) > 1) x[2] else x[1], character(1)),
    stringsAsFactors = FALSE
  )
  is_coded <- binarized$col_name %in% c("codedsiteid", "Seq_batch_short")
  binarized$short_name[is_coded] <- paste0(binarized$col_name[is_coded], "_", binarized$short_name[is_coded])

  rows <- merge(selected, binarized, by = "col_name", all.x = TRUE, sort = FALSE)
  rows$short_name[is.na(rows$short_name)] <- rows$col_name[is.na(rows$short_name)]
  rows$full_col_name[is.na(rows$full_col_name)] <- rows$col_name[is.na(rows$full_col_name)]
  rows <- rows[order(match(rows$category, unique(selected$category)),
                     match(rows$col_name, selected$col_name)), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

#' Pearson correlation and its p-value for every column pair, pairwise complete.
atac_pairwise_correlation <- function(x) {
  x <- as.matrix(x)
  r <- stats::cor(x, use = "pairwise.complete.obs")
  n <- crossprod(!is.na(x))
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  p <- 2 * stats::pt(-abs(t_stat), n - 2)
  diag(p) <- NA
  list(r = r, P = p)
}

# ---- labels, colours, theme ------------------------------------------------

#' Manuscript labels for group codes, timepoints and contrasts.
atac_pretty_labels <- function(x, short = FALSE) {
  labels <- as.character(x)
  labels <- stringr::str_remove_all(labels, "pre_and_|(?<!^)ADU.*_(?=pre)")
  labels <- stringr::str_replace_all(labels, c(
    "15_30_45_min" = "15-45 min", "3.5_4_hr" = "3.5-4 hr", "24_hr" = "24 hr",
    "vs" = "/", "_" = " ", "ADUEndur" = "EE", "ADUResist" = "RE", "ADUControl" = "CON"
  ))
  labels <- stringr::str_remove_all(labels, " (?=hr|min)")
  labels <- stringr::str_replace_all(labels, c(
    "(?<=^)EE(?= )" = "EE:", "(?<=^)RE(?= )" = "RE:", "(?<=^)Ctrl(?= )" = "Ctrl:"
  ))
  if (short) {
    labels <- stringr::str_replace_all(labels, c("pre exercise" = "Pre", "post " = "P"))
  }
  labels
}

#' Exercise-group colours keyed by both the code and the short label.
atac_group_colors <- function() {
  colors <- MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS
  c(colors, stats::setNames(colors, atac_pretty_labels(names(colors))))
}

atac_set2_colors <- function() RColorBrewer::brewer.pal(8, "Set2")

atac_theme <- function() {
  ggplot2::theme_classic(base_size = 6) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "snow2", linetype = "blank"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.line = ggplot2::element_line(linewidth = ggplot2::rel(0.5), lineend = "square"),
      axis.ticks = ggplot2::element_line(linewidth = ggplot2::rel(0.5), color = "black"),
      legend.margin = ggplot2::margin(0, 0, 0, 0, "cm")
    )
}

#' The group colour scale and legend key size the scatter panels share.
atac_group_color_scale <- function() {
  list(
    ggplot2::scale_color_manual(values = atac_group_colors(),
                                labels = atac_pretty_labels, name = "Modality"),
    ggplot2::theme(legend.key.size = grid::unit(2, "mm"))
  )
}

atac_heatmap_legend_params <- function() {
  list(title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
       labels_gp = grid::gpar(fontsize = 6),
       grid_width = grid::unit(0.3, "cm"), grid_height = grid::unit(0.3, "cm"))
}

#' cell_fun for ComplexHeatmap::Heatmap: writes label_mat[i, j] into every cell
#' whose p_mat[i, j] is below threshold.
atac_cell_labels <- function(p_mat, label_mat, threshold = ATAC_SIGNIFICANCE, fontsize = 3) {
  p_mat <- as.matrix(p_mat)
  label_mat <- as.matrix(label_mat)
  stopifnot(identical(dim(p_mat), dim(label_mat)))
  function(j, i, x, y, width, height, fill) {
    if (is.na(p_mat[i, j]) || p_mat[i, j] >= threshold || is.na(label_mat[i, j])) {
      return(invisible(NULL))
    }
    grid::grid.text(sprintf("%.1f", label_mat[i, j]), x, y,
                    gp = grid::gpar(fontsize = fontsize))
  }
}

# ---- shared ----------------------------------------------------------------

# A shared load, not seven separate ones: every panel that reads one of these
# reads it with the same arguments. SF1A and SF1B take the sample sheet alone;
# SF1C reads both PCA models; SF1D, SF1E and SF1F read the M0 model, which is
# the same atac_pca() call with the same covariates and the same protected
# terms, so it is the same object rather than a per-panel fit. SF1G shares
# nothing with the rest — the DA table is its own and stays inside it.
#
# Lazy and memoised rather than loaded at the top of the script: the QC object
# is about half a gigabyte and the raw counts are a bucket download, and a run
# that selects SF1A, SF1B or SF1G alone should touch neither.

atac_qc <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- atac_qc_objects()
    }
    cache
  }
})

sample_sheet <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- atac_sample_metadata(atac_qc())
    }
    cache
  }
})

# The voom step, which is where the raw peak counts are downloaded.
voom_log_cpm <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- atac_log_cpm(sample_sheet(), atac_qc())
    }
    cache
  }
})

pca_models <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- atac_pca_models(voom_log_cpm(), sample_sheet())
    }
    cache
  }
})

# ---- SF1A — sample prep per participant and timepoint ----------------------

# Three heatmaps side by side, one row per participant split by exercise group,
# one column per muscle timepoint. A white cell is a timepoint with no library.
# The middle heatmap prints the number of days wherever it is above zero.

sf1a <- function() {
  panel_init("SF1A")

  # ---- data ----

  sheet <- sample_sheet()

  # One participant per row, one timepoint per column.
  by_participant <- function(value_column) {
    wide <- sheet %>%
      dplyr::select(pid, randomGroupCode, Sex, Timepoint, dplyr::all_of(value_column)) %>%
      tidyr::pivot_wider(id_cols = c(pid, randomGroupCode, Sex), names_from = Timepoint,
                         values_from = dplyr::all_of(value_column), names_sort = TRUE)
    values <- as.matrix(wide[, setdiff(colnames(wide), c("pid", "randomGroupCode", "Sex")), drop = FALSE])
    rownames(values) <- wide$pid
    list(frame = wide, matrix = values)
  }

  isolation_batch <- by_participant("nuclei_isolation_batch")
  isolation_lag <- by_participant("nucleiIsolation2Tagmentation")
  sequencing_batch <- by_participant("Seq_batch_short")

  flowcell_levels <- unique(stats::na.omit(as.vector(sequencing_batch$matrix)))
  if (length(flowcell_levels) > length(atac_set2_colors())) {
    stop(length(flowcell_levels), " sequencing batches, more than the Set2 palette holds",
         call. = FALSE)
  }

  # ---- plot ----

  legend_params <- atac_heatmap_legend_params()
  name_gp <- grid::gpar(fontsize = 6)
  title_gp <- grid::gpar(fontsize = 7)
  row_height <- grid::unit(2.5, "mm")
  timepoint_labels <- atac_pretty_labels(colnames(isolation_batch$matrix))

  hm_isolation_batch <- ComplexHeatmap::Heatmap(
    isolation_batch$matrix, name = "isolation_batch",
    na_col = "white", row_split = isolation_batch$frame$randomGroupCode,
    cluster_rows = FALSE, cluster_columns = FALSE,
    col = rev(grDevices::rainbow(10)), column_labels = timepoint_labels,
    column_title = "Nuclei\nIsolation\nBatch", column_title_gp = title_gp,
    row_names_gp = name_gp, column_names_gp = name_gp, row_title_gp = title_gp,
    heatmap_legend_param = c(
      list(title = "Nuclei Isolation\nBatch (ordered)",
           at = round(seq(min(isolation_batch$matrix, na.rm = TRUE),
                          max(isolation_batch$matrix, na.rm = TRUE), length.out = 4)),
           direction = "horizontal", legend_width = grid::unit(1.75, "cm"),
           grid_height = grid::unit(0.25, "cm")),
      legend_params[c("title_gp", "labels_gp")]
    ),
    height = nrow(isolation_batch$matrix) * row_height,
    width = ncol(isolation_batch$matrix) * grid::unit(2.5, "mm"),
    left_annotation = ComplexHeatmap::rowAnnotation(
      df = data.frame(Group = isolation_batch$frame$randomGroupCode, Sex = isolation_batch$frame$Sex),
      col = list(Sex = MotrpacHumanPreSuspensionAnalysis::HUMAN_SEX_COLORS,
                 Group = MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS),
      simple_anno_size = grid::unit(3, "mm"), annotation_name_gp = name_gp,
      annotation_legend_param = legend_params
    )
  )

  hm_isolation_lag <- ComplexHeatmap::Heatmap(
    isolation_lag$matrix, name = "isolation_lag",
    na_col = "white", row_split = isolation_lag$frame$randomGroupCode,
    cluster_rows = FALSE, cluster_columns = FALSE,
    col = circlize::colorRamp2(c(0, 1), c("#CBFBF1", "#36BBA7")), column_labels = timepoint_labels,
    column_title = "ATAC From\nFrozen Nuclei", column_title_gp = title_gp,
    row_names_gp = name_gp, column_names_gp = name_gp, row_title_gp = title_gp,
    heatmap_legend_param = c(
      list(title = "Days Between\nNuclei Isolation\n& Tagmentation", at = c(1, 0),
           labels = c("≥1", "0"), color_bar = "discrete"),
      legend_params
    ),
    height = nrow(isolation_lag$matrix) * row_height,
    width = ncol(isolation_lag$matrix) * grid::unit(3.5, "mm"),
    layer_fun = function(j, i, x, y, width, height, fill) {
      days <- ComplexHeatmap::pindex(isolation_lag$matrix, i, j)
      show <- !is.na(days) & days > 0
      if (any(show)) {
        grid::grid.text(sprintf("%.0f", days[show]), x[show], y[show], gp = grid::gpar(fontsize = 5))
      }
    }
  )

  hm_sequencing_batch <- ComplexHeatmap::Heatmap(
    sequencing_batch$matrix, name = "sequencing_batch",
    na_col = "white", row_split = sequencing_batch$frame$randomGroupCode,
    cluster_rows = FALSE, cluster_columns = FALSE,
    col = stats::setNames(atac_set2_colors()[seq_along(flowcell_levels)], flowcell_levels),
    column_labels = timepoint_labels,
    column_title = "Sequence\nBatch", column_title_gp = title_gp,
    row_names_gp = name_gp, column_names_gp = name_gp, row_title_gp = title_gp,
    heatmap_legend_param = c(list(title = "Sequence\nBatch"), legend_params),
    height = nrow(sequencing_batch$matrix) * row_height,
    width = ncol(sequencing_batch$matrix) * grid::unit(2.5, "mm")
  )

  export_panel(function() {
    ComplexHeatmap::draw(
      hm_isolation_batch + hm_isolation_lag + hm_sequencing_batch,
      heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
      legend_grid_width = grid::unit(1, "cm"), legend_grouping = "original"
    )
  }, "SF1A")
}

# ---- SF1B — fraction of reads in peaks per timepoint and group -------------

# One box per timepoint over the libraries, with every library jittered on top.

sf1b <- function() {
  panel_init("SF1B")

  # ---- data ----

  sheet <- sample_sheet()

  # ---- plot ----

  p <- ggplot(sheet, aes(x = Timepoint, y = FRiP, color = randomGroupCode)) +
    atac_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ randomGroupCode) +
    xlab(NULL) +
    scale_x_discrete(labels = function(x) atac_pretty_labels(x, short = TRUE)) +
    geom_boxplot(outliers = FALSE, show.legend = FALSE) +
    geom_jitter(aes(shape = Timepoint), size = 1.5) +
    atac_group_color_scale() +
    scale_shape_manual(values = ATAC_TIMEPOINT_SHAPES) +
    theme(legend.position = "none")

  export_panel(p, "SF1B")
}

# ---- SF1C — PC correlation with covariates, before and after M0 adjustment --

# Rows are the covariates, factors expanded to one 0/1 row per level; columns
# are PC1-PC5 of the unadjusted PCA and of the PCA after the M0 covariates are
# regressed out with the exercise terms protected. A cell prints its Pearson r
# where p < 0.05.

sf1c <- function() {
  panel_init("SF1C")

  # ---- data ----

  sheet <- sample_sheet()
  models <- pca_models()
  covariate_rows <- atac_covariate_columns(sheet)

  # PC1-PC5 of each model, each column named "<PC (variance)>___<model label>".
  pc_columns <- Map(function(model, label) {
    scores <- model$scores[, c("library", names(model$variance)[1:5]), drop = FALSE]
    colnames(scores)[-1] <- paste0(model$variance[1:5], "___", label)
    scores
  }, models, ATAC_MODEL_LABELS[names(models)])
  pc_wide <- Reduce(function(a, b) merge(a, b, by = "library"), pc_columns)

  correlation_input <- merge(sheet[, c("library", covariate_rows$full_col_name), drop = FALSE],
                             pc_wide, by = "library")
  correlation_input$library <- NULL

  row_vars <- covariate_rows$full_col_name
  col_vars <- setdiff(colnames(correlation_input), row_vars)
  column_model <- sub("^.*___", "", col_vars)
  column_pc <- sub("___.*$", "", col_vars)

  correlation <- atac_pairwise_correlation(correlation_input)
  r <- correlation$r[row_vars, col_vars, drop = FALSE]
  p_value <- correlation$P[row_vars, col_vars, drop = FALSE]

  # ---- plot ----

  legend_params <- atac_heatmap_legend_params()
  adjustment_levels <- unique(column_model)
  annotation_df <- data.frame(
    `Correlation Inclusion` = "All Samples",
    samples = "All",
    adjustment = column_model,
    check.names = FALSE
  )

  hm <- ComplexHeatmap::Heatmap(
    r, name = "r",
    col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                               c("#444EB8", "#A2A7DB", "white", "#D38A90", "#B8444E")),
    na_col = "gray90",
    width = ncol(r) * grid::unit(2.5, "mm"), height = nrow(r) * grid::unit(2.5, "mm"),
    row_names_gp = grid::gpar(fontsize = 6), column_names_gp = grid::gpar(fontsize = 6),
    heatmap_legend_param = c(list(title = "r"), legend_params[c("title_gp", "labels_gp", "grid_width")]),
    column_title_gp = grid::gpar(fontsize = 8), row_title_gp = grid::gpar(fontsize = 8),
    column_title_side = "bottom",
    column_split = factor(gsub(" ", "\n", column_model), levels = unique(gsub(" ", "\n", column_model))),
    row_split = factor(covariate_rows$category, levels = unique(covariate_rows$category)),
    column_labels = column_pc,
    row_labels = covariate_rows$short_name,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    cell_fun = atac_cell_labels(p_value, r),
    row_names_max_width = grid::unit(2, "cm"),
    top_annotation = ComplexHeatmap::columnAnnotation(
      df = annotation_df,
      col = list(
        adjustment = stats::setNames(c("#E78AC3", "#A6D854", "#FFD92F", "#6699CC")[seq_along(adjustment_levels)],
                                     adjustment_levels),
        `Correlation Inclusion` = c("All Samples" = "#999933", "Within Group" = "#882255"),
        samples = c("All" = "#E5C494", MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS)
      ),
      annotation_name_gp = grid::gpar(fontsize = 6),
      annotation_legend_param = c(legend_params, list(nrow = 2)),
      simple_anno_size = grid::unit(3, "mm"),
      show_annotation_name = FALSE
    )
  )

  export_panel(function() {
    ComplexHeatmap::draw(hm, heatmap_legend_side = "right", annotation_legend_side = "top",
                         align_heatmap_legend = "heatmap_top")
  }, "SF1C")
}

# ---- SF1D — M0-adjusted PC1 against nuclei isolation batch -----------------

# One point per library. The trend lines and Spearman correlations are fitted
# over batches below ATAC_TREND_MAX_BATCH, pooled (tan) and per exercise group.

sf1d <- function() {
  panel_init("SF1D")

  # ---- data ----

  m0 <- pca_models()$M0

  scores <- m0$scores %>%
    dplyr::mutate(
      Timepoint = forcats::fct_relabel(Timepoint, function(x) atac_pretty_labels(x, short = TRUE)),
      batch_index = as.numeric(factor(nuclei_isolation_batch))
    )
  early_batches <- scores %>% dplyr::filter(nuclei_isolation_batch < ATAC_TREND_MAX_BATCH)

  # ---- plot ----

  pooled_color <- "#E5C494"

  p <- ggplot(scores, aes(x = factor(nuclei_isolation_batch), y = PC1)) +
    atac_theme() +
    labs(x = "Nuclei Isolation Batch", y = paste("~M0", m0$variance[["PC1"]])) +
    geom_point(aes(color = randomGroupCode, shape = Timepoint), size = 2, alpha = 0.5, stroke = 0) +
    atac_group_color_scale() +
    scale_shape_manual(values = ATAC_TIMEPOINT_SHAPES) +
    theme(legend.position = "none") +
    geom_smooth(data = early_batches, aes(x = batch_index), method = "lm", se = FALSE,
                color = pooled_color) +
    geom_smooth(data = early_batches, aes(x = batch_index, color = randomGroupCode),
                method = "lm", se = FALSE) +
    ggpubr::stat_cor(aes(x = batch_index), method = "spearman", cor.coef.name = "rho",
                     color = pooled_color, label.y.npc = 0.27, label.x.npc = 0.72,
                     p.digits = 2, p.accuracy = 1e-10) +
    ggpubr::stat_cor(aes(x = batch_index, color = randomGroupCode), method = "spearman",
                     cor.coef.name = "rho", label.y.npc = 0.2, label.x.npc = 0.72,
                     p.accuracy = 0.001)

  export_panel(p, "SF1D")
}

# ---- SF1E — M0-adjusted PC2, frozen against fresh nuclei -------------------

# Left: one box per group and nuclei state. Right: the same, split by timepoint.

sf1e <- function() {
  panel_init("SF1E")

  # ---- data ----

  m0 <- pca_models()$M0

  # ---- plot ----

  base <- ggplot(m0$scores, aes(x = randomGroupCode, y = PC2, color = nuclei2atac_lag, shape = Timepoint)) +
    atac_theme()

  finish <- function(p) {
    p +
      geom_point(position = position_jitterdodge(), size = 1.25, stroke = 0) +
      guides(color = guide_legend(title = "Use of Frozen Nuclei", override.aes = list(size = 2)),
             shape = "none") +
      theme(legend.position = "bottom", legend.margin = margin(),
            legend.key.width = unit(2, "mm"), legend.key.height = unit(1, "mm")) +
      labs(x = NULL, y = paste("~M0", m0$variance[["PC2"]])) +
      scale_shape_manual(values = ATAC_TIMEPOINT_SHAPES) +
      scale_x_discrete(labels = atac_pretty_labels)
  }

  pooled <- finish(base + geom_boxplot(aes(shape = NULL), outliers = FALSE, show.legend = FALSE, linewidth = 0.25))
  by_timepoint <- finish(base + geom_boxplot(outliers = FALSE, show.legend = FALSE, linewidth = 0.25))

  p <- ggpubr::ggarrange(pooled, by_timepoint, nrow = 1, common.legend = TRUE,
                         widths = c(0.7, 1), legend = "bottom")

  export_panel(p, "SF1E")
}

# ---- SF1F — M0-adjusted PC1 against PC2 ------------------------------------

# With a 95% ellipse per exercise group and timepoint.

sf1f <- function() {
  panel_init("SF1F")

  # ---- data ----

  m0 <- pca_models()$M0

  scores <- m0$scores %>%
    dplyr::mutate(Timepoint = forcats::fct_relabel(Timepoint, function(x) atac_pretty_labels(x, short = TRUE)))

  # ---- plot ----

  p <- ggplot(scores, aes(x = PC1, y = PC2, color = randomGroupCode, shape = Timepoint,
                          linetype = Timepoint, linewidth = Timepoint)) +
    atac_theme() +
    geom_point(size = 1.5, stroke = 0, alpha = 0.8) +
    stat_ellipse(alpha = 0.8, segments = 100) +
    theme(aspect.ratio = 1) +
    labs(x = paste("~M0", m0$variance[["PC1"]]), y = paste("~M0", m0$variance[["PC2"]])) +
    scale_linewidth_manual(values = rep(0.25, length(ATAC_TIMEPOINT_SHAPES))) +
    scale_linetype_manual(values = c("solid", "1342", "dashed", "dotted")) +
    atac_group_color_scale() +
    scale_shape_manual(values = ATAC_TIMEPOINT_SHAPES) +
    theme(legend.position = "right", legend.key.width = unit(5, "mm")) +
    guides(linetype = guide_legend(order = 1, override.aes = list(alpha = 1)),
           linewidth = guide_legend(order = 1),
           shape = guide_legend(order = 2, override.aes = list(size = 2)),
           color = guide_legend(order = 3, override.aes = list(linetype = 0, size = 2, alpha = 1)))

  export_panel(p, "SF1F")
}

# ---- SF1G — significant peaks per exercise arm and timepoint ---------------

# Left: one tile per contrast type with the number of peaks at adjusted
# p < 0.05. Right: those peaks split into difference-in-changes only, common,
# and within-group only. A cell with no significant peak reads 0.
#
# The DA table is read here and nowhere else in this figure, so it is not one of
# the shared accessors above.

sf1g <- function() {
  panel_init("SF1G")

  # ---- data ----

  da <- atac_differential_analysis()

  timepoint_levels <- levels(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER$Timepoint)
  post_timepoints <- setdiff(intersect(timepoint_levels, unique(as.character(da$Timepoint))), "pre_exercise")

  significant <- da %>%
    dplyr::filter(adj_p_value < ATAC_SIGNIFICANCE, contrast_type %in% ATAC_CONTRAST_TYPES) %>%
    dplyr::transmute(
      contrast_type = factor(contrast_type, levels = ATAC_CONTRAST_TYPES),
      randomGroupCode = factor(as.character(randomGroupCode), levels = ATAC_EXERCISE_GROUPS),
      Timepoint = factor(as.character(Timepoint), levels = post_timepoints),
      feature_id = as.character(feature_id)
    )
  stopifnot(!anyNA(significant$randomGroupCode), !anyNA(significant$Timepoint))

  peak_counts <- significant %>%
    dplyr::count(randomGroupCode, Timepoint, contrast_type, name = "n") %>%
    tidyr::complete(randomGroupCode, Timepoint, contrast_type, fill = list(n = 0))

  with_controls <- function(df) df$feature_id[df$contrast_type == "exercise_with_controls"]
  no_controls <- function(df) df$feature_id[df$contrast_type == "exercise_no_controls"]

  peak_overlap <- significant %>%
    dplyr::group_by(randomGroupCode, Timepoint) %>%
    dplyr::group_modify(function(df, key) data.frame(
      dd_only = length(setdiff(with_controls(df), no_controls(df))),
      common = length(intersect(with_controls(df), no_controls(df))),
      within_group_only = length(setdiff(no_controls(df), with_controls(df)))
    )) %>%
    dplyr::ungroup() %>%
    tidyr::complete(randomGroupCode, Timepoint,
                    fill = list(dd_only = 0, common = 0, within_group_only = 0)) %>%
    dplyr::arrange(randomGroupCode, Timepoint)

  message(sprintf(
    "        %.1f%% of difference-in-changes peaks are also significant within group",
    100 * sum(peak_overlap$common) / (sum(peak_overlap$common) + sum(peak_overlap$dd_only))
  ))

  # ---- plot ----

  timepoint_axis <- scale_y_discrete(labels = function(x) atac_pretty_labels(x, short = TRUE),
                                     expand = expansion(0))

  p_counts <- ggplot(peak_counts, aes(y = forcats::fct_rev(Timepoint), x = contrast_type, fill = n, label = n)) +
    atac_theme() +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(fontface = "bold") +
    facet_wrap(~ randomGroupCode) +
    scale_fill_gradientn(colours = c("#EFEDF5", "#BCBDDC", "#756BB1")) +
    scale_x_discrete(labels = c("difference-in-\nchanges", "within-\ngroup"), expand = expansion(0)) +
    timepoint_axis +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    ggtitle("Count of Significant Peaks")

  overlap_positions <- c(dd_only = 0.125, common = 0.5, within_group_only = 0.875)
  p_overlap <- peak_overlap %>%
    tidyr::pivot_longer(c(dd_only, common, within_group_only), names_to = "category", values_to = "n") %>%
    dplyr::mutate(x_val = overlap_positions[category],
                  category = factor(category, levels = names(overlap_positions))) %>%
    ggplot(aes(x = x_val, y = forcats::fct_rev(Timepoint), label = n, fill = category)) +
    atac_theme() +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(fontface = "bold") +
    facet_wrap(~ randomGroupCode) +
    labs(y = NULL, x = NULL) +
    timepoint_axis +
    scale_x_continuous(breaks = unname(overlap_positions),
                       labels = c("difference-in-\nchanges only", "common", "within-group\nonly"),
                       expand = expansion(0)) +
    scale_fill_manual(values = atac_set2_colors()[-(1:3)]) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    ggtitle("Distinct & Overlapping Significant Peaks")

  p <- ggpubr::ggarrange(p_counts, p_overlap, nrow = 1, align = "h", widths = c(0.6, 1))

  export_panel(p, "SF1G")
}

run_panels(list(SF1A = sf1a, SF1B = sf1b, SF1C = sf1c, SF1D = sf1d,
                SF1E = sf1e, SF1F = sf1f, SF1G = sf1g))
