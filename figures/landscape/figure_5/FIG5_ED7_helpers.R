# FIG5_ED7.R — the HOMER transcription-factor enrichment, and the bubble heatmap
# it is drawn as, as Figure 5, Extended Data 7 and ST5a read them.
#
# One file because more than one figure reads each half:
#   the enrichment functions   every FIG5 panel, both ED7 panels, tables/ST5.R
#   tf_bubble_heatmap()        FIG5D, ED7A, ED7B
#
# Four steps, in order:
#
#   tf_comparisons()        the 22 tissue x modality x timepoint HOMER runs
#   tf_enrichment()         motif x comparison q-values
#   tf_annotation()         motif -> gene symbol -> transcript / prot-pr / prot-ol ids
#   tf_da_stats()           logFC and adj_p for a feature set over one tissue's grid
#
# Sourced after lib/panel_export.R, which defines panel_source() and
# landscape_root().
#
# Everything is namespace-qualified: a script sourcing this file inherits no
# attached packages from it.

# ---- the comparison grid ---------------------------------------------------

# The manuscript's short timepoint codes, keyed by the DA tables' Timepoint.
TF_TIMEPOINT_CODES <- c(
  "during_20_min"     = "D20M",
  "during_40_min"     = "D40M",
  "post_10_min"       = "P10M",
  "post_15_30_45_min" = "P15-45M",
  "post_3.5_4_hr"     = "P3.5/4H",
  "post_24_hr"        = "P24H"
)

TF_MODALITIES <- c("ADUEndur" = "EE", "ADUResist" = "RE")

TF_TISSUES <- c("muscle" = "Muscle", "adipose" = "Adipose", "blood" = "Blood")

# The slug HOMER's output files carry for each timepoint. The DEG lists were
# written under this naming and the results were never renamed, so it is the
# only link between a file on disk and the contrast it came from.
TF_HOMER_TIMEPOINT_SLUGS <- c(
  "during_20_min"     = "during20",
  "during_40_min"     = "during40",
  "post_10_min"       = "immpost",
  "post_15_30_45_min" = "early",
  "post_3.5_4_hr"     = "mid",
  "post_24_hr"        = "late"
)

TF_HOMER_MODALITY_SLUGS <- c("ADUEndur" = "endur", "ADUResist" = "resist")

#' The 22 HOMER runs, in manifest column order.
#'
#' Muscle then adipose (RE, then EE, over the three post timepoints), then blood
#' (RE over four, EE over six). Rows are ordered so that column i of every matrix
#' in this file is row i here.
#'
#' @returns A data.frame with one row per comparison: `comparison` (the key used
#'   as a matrix column name), `tissue`, `Tissue`, `randomGroupCode`, `Modality`,
#'   `Timepoint`, `tp_code` and `homer_file`.
tf_comparisons <- function() {
  grid <- function(tissue, modality, timepoints) {
    data.frame(
      tissue = tissue,
      randomGroupCode = modality,
      Timepoint = timepoints,
      stringsAsFactors = FALSE
    )
  }

  post <- c("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")
  blood_re <- c("post_10_min", post)
  blood_ee <- c("during_20_min", "during_40_min", "post_10_min", post)

  cmp <- rbind(
    grid("muscle",  "ADUResist", post),
    grid("muscle",  "ADUEndur",  post),
    grid("adipose", "ADUResist", post),
    grid("adipose", "ADUEndur",  post),
    grid("blood",   "ADUResist", blood_re),
    grid("blood",   "ADUEndur",  blood_ee)
  )

  cmp$Tissue <- unname(TF_TISSUES[cmp$tissue])
  cmp$Modality <- unname(TF_MODALITIES[cmp$randomGroupCode])
  cmp$tp_code <- unname(TF_TIMEPOINT_CODES[cmp$Timepoint])
  cmp$comparison <- paste(cmp$Tissue, cmp$Modality, cmp$tp_code, sep = "|")
  cmp$homer_file <- sprintf(
    "%s_%s_%s_sig_knownResults.txt",
    cmp$tissue,
    TF_HOMER_TIMEPOINT_SLUGS[cmp$Timepoint],
    TF_HOMER_MODALITY_SLUGS[cmp$randomGroupCode]
  )

  rownames(cmp) <- NULL
  cmp[, c("comparison", "tissue", "Tissue", "randomGroupCode", "Modality",
          "Timepoint", "tp_code", "homer_file")]
}

#' Tissue, modality and timepoint as a heatmap column annotation frame.
#'
#' Rownames are the comparison keys, so pheatmap and ComplexHeatmap both match a
#' reordered matrix by name rather than by position.
tf_column_annotation <- function(comparisons = tf_comparisons()) {
  df <- data.frame(
    Timepoint = factor(comparisons$tp_code, levels = unname(TF_TIMEPOINT_CODES)),
    Modality  = factor(comparisons$Modality, levels = c("EE", "RE")),
    Tissue    = factor(comparisons$Tissue, levels = c("Muscle", "Adipose", "Blood")),
    row.names = comparisons$comparison
  )
  df
}

#' Annotation colours, from the Analysis package rather than restated.
#'
#' The legacy hardcoded the modality and timepoint hexes in the two panels it
#' rewrote for ComplexHeatmap and read the package vectors in the rest, so one
#' figure carried two spellings of the same three scales.
tf_annotation_colors <- function() {
  tissue <- MotrpacHumanPreSuspensionAnalysis::HUMAN_TISSUE_COLORS[
    c("muscle", "adipose", "blood")]
  names(tissue) <- c("Muscle", "Adipose", "Blood")

  modality <- MotrpacHumanPreSuspensionAnalysis::HUMAN_EXERCISE_GROUP_COLORS[
    c("ADUEndur", "ADUResist")]
  names(modality) <- c("EE", "RE")

  timepoint <- MotrpacHumanPreSuspensionAnalysis::HUMAN_ACUTE_TIMEPOINT_COLORS[
    names(TF_TIMEPOINT_CODES)]
  names(timepoint) <- unname(TF_TIMEPOINT_CODES)

  list(Tissue = tissue, Modality = modality, Timepoint = timepoint)
}

# ---- the HOMER results -----------------------------------------------------

# HOMER reports q down to 1e-4 and writes anything smaller as 0. The legacy
# floored the zeros at 1e-5 for muscle and blood only, because adipose happens to
# have none; the floor is applied to every tissue here so that an adipose zero in
# a future run cannot become an infinite -log10.
TF_QVALUE_FLOOR <- 1e-05

#' Path to one vendored HOMER knownResults file.
tf_homer_path <- function(homer_file) {
  panel_source("figure_5",
               file.path("Precovid_DEG_HOMER_KnownTF_Results", homer_file))
}

#' One HOMER knownResults table, one row per motif.
#'
#' HOMER emits a motif more than once when several of its variants score; the
#' legacy kept the first, which is the best-scoring one because the file is
#' sorted by p-value. Kept.
tf_homer_results <- function(homer_file) {
  res <- utils::read.csv(tf_homer_path(homer_file), sep = "\t", check.names = FALSE)
  required <- c("Motif Name", "q-value (Benjamini)")
  missing <- setdiff(required, colnames(res))
  if (length(missing) > 0) {
    stop("HOMER file ", homer_file, " has no column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  res <- res[!duplicated(res[["Motif Name"]]), , drop = FALSE]
  rownames(res) <- res[["Motif Name"]]
  res
}

#' Every HOMER run, as a motif x comparison q-value matrix.
#'
#' @returns A list with `q` (motifs x 22, floored), `q_raw`, `motifs`,
#'   `comparisons` and `results`.
tf_enrichment <- function(comparisons = tf_comparisons()) {
  results <- lapply(comparisons$homer_file, tf_homer_results)
  names(results) <- comparisons$comparison

  motifs <- rownames(results[[1]])
  for (i in seq_along(results)) {
    if (!setequal(rownames(results[[i]]), motifs)) {
      stop("HOMER file ", comparisons$homer_file[i],
           " scores a different motif set than ", comparisons$homer_file[1],
           call. = FALSE)
    }
  }

  q_raw <- vapply(results,
                  function(res) res[motifs, "q-value (Benjamini)"],
                  numeric(length(motifs)))
  dimnames(q_raw) <- list(motifs, comparisons$comparison)
  q <- q_raw
  q[q == 0] <- TF_QVALUE_FLOOR

  list(
    q = q,
    # Unfloored, for ST5a. A panel needs the floor so a reported 0 does not
    # become an infinite -log10; a table reports what HOMER reported.
    q_raw = q_raw,
    motifs = motifs,
    comparisons = comparisons,
    # Kept in HOMER's own row order, which is by ascending p-value. Selecting the
    # strongest motifs of a comparison off the matrix instead would break ties by
    # motif name; off these it breaks them the way the published figure did.
    results = results
  )
}

#' The `n` strongest motifs of every comparison, unioned.
#'
#' Far fewer than 22n motifs: the same factors top several comparisons.
tf_top_motifs <- function(enrichment, n = 5) {
  per_comparison <- lapply(enrichment$results, function(res) {
    res[order(res[["q-value (Benjamini)"]]), "Motif Name"][seq_len(n)]
  })
  Reduce(union, per_comparison)
}

#' Motifs significantly enriched anywhere in one tissue.
tf_enriched_motifs <- function(enrichment, tissue, alpha = 0.05) {
  cols <- enrichment$comparisons$comparison[enrichment$comparisons$Tissue == tissue]
  if (length(cols) == 0) {
    stop("no comparison for tissue '", tissue, "'", call. = FALSE)
  }
  q <- enrichment$q[, cols, drop = FALSE]
  rownames(q)[apply(q, 1, min) < alpha]
}

# ---- the motif -> gene annotation ------------------------------------------

#' Motif name to human gene symbol, and from there to feature ids.
#'
#' The default reading of a HOMER motif name is everything before the first
#' parenthesis, uppercased: "Foxo3(Forkhead)/..." -> "FOXO3". That resolves 299
#' of the 471 motifs to a symbol this study measured. `tfproanno.RDS`, a curated
#' motif -> gene table carried over from PASS1B, overrides it for 402 motifs and
#' brings the total to 416: it is what maps "ETS:RUNX(ETS,Runt)/..." to RUNX1 and
#' "TATA-Box(TBP)/Promoter/Homer" to TBP.
#'
#' Its Ensembl and per-tissue protein-id columns are rat and are dropped; only
#' Gene.Name is read, and every id below is resolved against
#' HUMAN_FEATURE_TO_GENE.
#'
#' The legacy applied 15 hand corrections to that table by ROW NUMBER, against a
#' version of it with different contents: 13 landed on unrelated motifs (row 11
#' is Zic3, and was overwritten with "JUN") and 2 were past its 402 rows and
#' appended NA-named rows. They are not reproduced.
#'
#' @param motifs Motif names, in the order the returned frame should carry.
#' @returns A data.frame keyed by motif: `gene_symbol`, `ensembl_gene`,
#'   `prot_pr_id`, `prot_ol_id`. Unresolved ids are NA.
tf_annotation <- function(motifs) {
  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  feature_to_gene$gene_symbol <- as.character(feature_to_gene$gene_symbol)

  symbol <- toupper(gsub("\\(.*", "", motifs))

  curated <- readRDS(panel_source("figure_5", "tfproanno.RDS"))
  hit <- motifs %in% rownames(curated)
  symbol[hit] <- toupper(as.character(curated[motifs[hit], "Gene.Name"]))

  # First id per symbol per assay, as the legacy took. A symbol carrying several
  # features in one assay is a real case (isoforms, several Olink panels); which
  # one is chosen only affects the single-id lookups, and every count below goes
  # through the symbol.
  first_id_by_symbol <- function(assay) {
    rows <- feature_to_gene[feature_to_gene$assay == assay, ]
    rows <- rows[!duplicated(rows$gene_symbol), ]
    stats::setNames(as.character(rows$feature_id), rows$gene_symbol)
  }
  first_ensembl_by_symbol <- function() {
    rows <- feature_to_gene[!is.na(feature_to_gene$ensembl_gene) &
                              nzchar(as.character(feature_to_gene$ensembl_gene)), ]
    rows <- rows[!duplicated(rows$gene_symbol), ]
    stats::setNames(as.character(rows$ensembl_gene), rows$gene_symbol)
  }

  measured <- unique(feature_to_gene$gene_symbol)
  anno <- data.frame(
    gene_symbol  = ifelse(symbol %in% measured, symbol, NA_character_),
    ensembl_gene = unname(first_ensembl_by_symbol()[symbol]),
    prot_pr_id   = unname(first_id_by_symbol("prot-pr")[symbol]),
    prot_ol_id   = unname(first_id_by_symbol("prot-ol")[symbol]),
    row.names    = motifs,
    stringsAsFactors = FALSE
  )

  unresolved <- sum(is.na(anno$gene_symbol))
  message(sprintf(
    "        %d of %d motifs map to a gene symbol this study measured (%d do not)",
    nrow(anno) - unresolved, nrow(anno), unresolved))
  anno
}

#' The prot-ph feature ids belonging to each TF, as a symbol -> ids list.
tf_phosphosites <- function(symbols) {
  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  ph <- feature_to_gene[feature_to_gene$assay == "prot-ph", ]
  split(as.character(ph$feature_id), as.character(ph$gene_symbol))[symbols]
}

# ---- one parse, one annotation, per run ------------------------------------

# Every panel but FIG5A reads the 22 HOMER files through tf_enrichment(), and
# every panel from FIG5D on reads tf_annotation() over the full motif list. Both
# are called with the same arguments everywhere — FIG5, ED7 and tables/ST5.R
# alike — so both are memoised here rather than once per script. Lazy, so
# `Rscript FIG5.R FIG5A` still reads no HOMER file.
tf_enrichment_cached <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- tf_enrichment()
    }
    cache
  }
})

tf_annotation_cached <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- tf_annotation(tf_enrichment_cached()$motifs)
    }
    cache
  }
})

# ---- differential abundance over the comparison grid -----------------------

#' logFC and adj_p for a feature set, over one tissue's exercise-with-controls grid.
#'
#' Replaces roughly 500 lines in which each tissue x ome x contrast cell was
#' filled by a separate `DA[DA$feature_id %in% id & DA$contrast %in% "<the full
#' 160-character contrast string>", "logFC"][[1]]`. The columns are read off
#' `contrast_type`, `randomGroupCode` and `Timepoint`, which the DA tables have
#' carried since 2.0.0, so a tissue gaining or losing a timepoint changes the
#' matrix rather than silently missing a column.
#'
#' @param da A per-tissue DA table.
#' @param feature_ids Features to keep. NULL keeps every feature in `da`.
#' @param comparisons The rows of tf_comparisons() for this tissue.
#' @returns A list of two matrices, `logFC` and `adj_p`, features x comparisons.
tf_da_stats <- function(da, feature_ids, comparisons) {
  da <- as.data.frame(da)
  da <- da[da$contrast_type == "exercise_with_controls", ]
  da$feature_id <- as.character(da$feature_id)
  if (!is.null(feature_ids)) {
    da <- da[da$feature_id %in% feature_ids, ]
  }

  key <- paste(da$randomGroupCode, da$Timepoint)
  want <- paste(comparisons$randomGroupCode, comparisons$Timepoint)
  absent <- setdiff(want, unique(key))
  if (length(absent) > 0) {
    stop("this DA table has no exercise-with-controls contrast for: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  da <- da[key %in% want, ]
  da$comparison <- comparisons$comparison[match(paste(da$randomGroupCode,
                                                      da$Timepoint), want)]

  features <- unique(da$feature_id)
  fill <- function(value_col) {
    m <- matrix(NA_real_, nrow = length(features), ncol = nrow(comparisons),
                dimnames = list(features, comparisons$comparison))
    m[cbind(match(da$feature_id, features), match(da$comparison, colnames(m)))] <-
      da[[value_col]]
    m
  }

  logFC <- fill("logFC")
  adj_p <- fill("adj_p_value")
  if (anyNA(logFC) || anyNA(adj_p)) {
    stop(sum(is.na(logFC) | is.na(adj_p)),
         " feature x contrast cell(s) are absent from the DA table; the grid is ",
         "not rectangular", call. = FALSE)
  }
  list(logFC = logFC, adj_p = adj_p)
}

#' Features significant in any exercise-with-controls contrast of one DA table.
tf_significant_features <- function(da, alpha = 0.05) {
  da <- as.data.frame(da)
  da <- da[da$contrast_type == "exercise_with_controls" & da$adj_p_value < alpha, ]
  unique(as.character(da$feature_id))
}

#' Every feature id in a DA table's exercise-with-controls contrasts.
tf_measured_features <- function(da) {
  da <- as.data.frame(da)
  unique(as.character(da$feature_id[da$contrast_type == "exercise_with_controls"]))
}

#' The smallest adj p among a motif's features, per comparison.
#'
#' A motif resolves to one gene symbol and a symbol to one or more features of an
#' assay - one transcript, usually several phosphosites. The strongest of them
#' stands for the factor, which is what the legacy did for phosphorylation. For
#' transcripts it took `grep(ensembl, feature_id)[1]` instead, an unanchored
#' substring match that would also select a longer id; every symbol resolves to
#' exactly one transcript feature, so min() over the same set agrees with it
#' today and cannot pick the wrong one tomorrow.
#'
#' @param annotation tf_annotation() over the motif list, in row order.
#' @param assay An assay id in HUMAN_FEATURE_TO_GENE.
#' @param da The tissue's DA table for that assay.
#' @param comparisons The rows of tf_comparisons() for that tissue.
#' @returns A motifs x comparisons matrix. NA where the factor has no feature.
tf_min_adj_p_by_motif <- function(annotation, assay, da, comparisons) {
  feature_to_gene <- as.data.frame(
    MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  rows <- feature_to_gene[feature_to_gene$assay == assay, ]
  symbol_of <- stats::setNames(as.character(rows$gene_symbol),
                               as.character(rows$feature_id))

  ids <- intersect(names(symbol_of), tf_measured_features(da))
  stats <- tf_da_stats(da, ids, comparisons)
  by_symbol <- split(rownames(stats$adj_p), symbol_of[rownames(stats$adj_p)])

  motifs <- rownames(annotation)
  out <- matrix(NA_real_, nrow = length(motifs), ncol = nrow(comparisons),
                dimnames = list(motifs, comparisons$comparison))
  for (i in seq_along(motifs)) {
    symbol <- annotation$gene_symbol[i]
    if (is.na(symbol) || is.null(by_symbol[[symbol]])) next
    out[i, ] <- apply(stats$adj_p[by_symbol[[symbol]], , drop = FALSE], 2, min)
  }
  out
}

# ---- shared labelling ------------------------------------------------------

#' "SYMBOL S663 S665" out of a prot-ph feature id and its gene symbol.
#'
#' The legacy stripped every lowercase s, t and y from the whole label with
#' gsub("s|t|y", ""), which is right for the residue markers and would silently
#' eat a letter from any symbol that carried one. Same rule as FIG2E.
tf_phosphosite_label <- function(feature_id, gene_symbol) {
  sites <- sub("^.*_", "", feature_id)
  sites <- gsub("([A-Z][0-9]+)[a-z]", "\\1 ", sites)
  trimws(paste(gene_symbol, trimws(sites)))
}

#' Motif name trimmed to what the figures label rows with: "Foxo3", "ETS:RUNX".
tf_motif_label <- function(motifs) {
  gsub("\\(.*", "", motifs)
}

# ============================================================================
# the bubble heatmap
# ============================================================================
#
# A grid of empty cells with a filled circle in each: area carries -log10(adj_p),
# colour carries log2FC. FIG5D, ED7A and ED7B are drawn this way. The legacy drew
# it with the bubbleHeatmap package and then replaced that with ComplexHeatmap so
# the panel would not have to be composited against the rest of the figure in
# Illustrator. The ComplexHeatmap version is what is ported; bubbleHeatmap is not
# a dependency of this repo.

TF_BUBBLE_CELL_MM <- 6
TF_BUBBLE_MAX_MM <- 5
# Where the colour ramp saturates. The legacy used 1; at 0.5 the same responses
# spread across more of the ramp, and it is the bound the PLIER bubble panels
# use. FIG5D, ED7A and ED7B all read it.
TF_BUBBLE_LOGFC_LIM <- 0.5

#' The log2FC colour ramp both the cells and the legend read.
tf_bubble_color_fun <- function(limit = TF_BUBBLE_LOGFC_LIM) {
  circlize::colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))
}

#' A bubble heatmap over two aligned matrices.
#'
#' @param logFC,adj_p Same dimnames, same order. Rows are drawn top to bottom in
#'   the order given; neither axis is clustered.
#' @param top_annotation A ComplexHeatmap column annotation, or NULL.
#' @param outline_alpha Outline the bubbles whose adjusted p clears this, so the
#'   panel marks significance rather than leaving it to be read off the area. NA
#'   draws no outline.
tf_bubble_heatmap <- function(logFC, adj_p, top_annotation = NULL,
                              outline_alpha = NA) {
  stopifnot(identical(dimnames(logFC), dimnames(adj_p)))
  col_fun <- tf_bubble_color_fun()
  size <- pmin(-log10(adj_p), TF_BUBBLE_MAX_MM)
  cell <- grid::unit(TF_BUBBLE_CELL_MM, "mm")

  ComplexHeatmap::Heatmap(
    matrix = matrix(NA_real_, nrow = nrow(logFC), ncol = ncol(logFC),
                    dimnames = dimnames(logFC)),
    width = cell * ncol(logFC),
    height = cell * nrow(logFC),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(fill = NA, col = "black"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    show_heatmap_legend = FALSE,
    top_annotation = top_annotation,
    cell_fun = function(j, i, x, y, width, height, fill) {
      outlined <- !is.na(outline_alpha) && !is.na(adj_p[i, j]) &&
        adj_p[i, j] < outline_alpha
      grid::grid.circle(
        x = x, y = y,
        r = grid::unit(size[i, j], "mm") / 2,
        gp = grid::gpar(fill = col_fun(logFC[i, j]),
                        col = if (outlined) "black" else NA,
                        lwd = 1)
      )
    }
  )
}

#' The legends a bubble heatmap needs: one for colour, one for area, and — when
#' the heatmap outlines its significant bubbles — one saying so.
tf_bubble_legends <- function(outline_alpha = NA) {
  legends <- list(
    ComplexHeatmap::Legend(
      title = "log2FC",
      col_fun = tf_bubble_color_fun(),
      at = c(-TF_BUBBLE_LOGFC_LIM, 0, TF_BUBBLE_LOGFC_LIM)
    ),
    ComplexHeatmap::Legend(
      title = "-log10(adj p)",
      type = "points",
      pch = 16,
      size = grid::unit(0:TF_BUBBLE_MAX_MM, "mm"),
      labels = c(as.character(0:(TF_BUBBLE_MAX_MM - 1)),
                 paste0(TF_BUBBLE_MAX_MM, "+"))
    )
  )
  if (!is.na(outline_alpha)) {
    legends[[length(legends) + 1L]] <- ComplexHeatmap::Legend(
      title = "",
      type = "points",
      pch = 1,
      size = grid::unit(3, "mm"),
      labels = paste("adj_p <", outline_alpha),
      legend_gp = grid::gpar(col = "black")
    )
  }
  legends
}
