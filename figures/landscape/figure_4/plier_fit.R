#!/usr/bin/env Rscript
# The two cross-tissue PLIER fits — the input FIG4G, FIG4H and all ten ED6
# panels read.
#
# One PLIER per ome over one matrix holding every tissue at once:
#
#   rna     14,428 genes measured in adipose, blood and muscle, each tissue
#           z-scored on its own, against the nine gene-set databases in
#           MOLECULAR_SIGNATURES (14,200 sets, 8,482 with at least 10 genes)
#   metab   the metabolites shared by all three tissues, each tissue's matrix
#           z-scored on its own, against the RefMet classes (162 sets)
#
# and, off each fit, the LV x comparison response grid every bubble panel draws.
#
# A fit rather than panel code because twelve panels read these fits and the
# provenance of a figure's input is part of the provenance of the figure, plus
# one reason this fit has on its own: a PLIER fit's latent variables are
# NUMBERED, and every panel selects on those numbers (figure_4/plier_lvs.env).
# Refitting per panel would let two panels draw two different LV 33s. The legacy
# wrote both fits to .rds files in the analyst's working directory
# (`crosstissuernacombofin_plierresult_081026.rds`, `metab_plier_081426.rds`)
# and neither file is in any repo.
#
# Needs consortium data access: it reads individual-level QC matrices out of
# MotrpacHumanPreSuspensionData and the molecular signatures out of the Analysis
# package. THE RNA ARM IS EXPENSIVE — a bit over two hours and ~8.5 GB of
# memory. The metabolomics arm is minutes. Both are idempotent: a fit already on
# disk is reused, never refitted, unless PLIER_FORCE=TRUE.
#
# Parameters: figure_4/plier.env, plus PLIER_DIR from config/landscape.env.
# The hand-set latent-variable selections are figure_4/plier_lvs.env's and are not
# read here except by the report below.
#
#   Rscript figures/landscape/analysis/04_plier.R
#   PLIER_FORCE=TRUE Rscript figures/landscape/analysis/04_plier.R
#
# The `report` sub-command prints what is in an EXISTING fit, one row per latent
# variable, and fits nothing. It exists because of the problem
# figure_4/plier_lvs.env's header describes: the LV numbers those twelve panels
# select on were read off a fit that is not in any repo, and a refit renumbers
# them. The question "which latent variables should this figure draw now" has to
# be answerable without opening the fit by hand.
#
#   Rscript figures/landscape/analysis/04_plier.R report rna   [n_pathways]
#   Rscript figures/landscape/analysis/04_plier.R report metab [n_pathways]

suppressPackageStartupMessages({
  library(dplyr)
  library(PLIER)
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R for landscape_root(), and for the config/landscape.env read it
# performs on load — PLIER_DIR comes from there.
source(file.path(here, "..", "lib", "panel_export.R"))

# ---- figure_4/plier.env -------------------------------------------------------

# The same line parser load_landscape_env() uses, over this fit's own file.
# Values already exported WIN, which is how a one-off override works. Anything
# still carrying a substitution after ${LANDSCAPE_ROOT} is resolved cannot be
# evaluated here, and setting the literal is worse than leaving the variable
# unset for the code default to handle.
load_plier_env <- function() {
  root <- landscape_root()
  path <- file.path(root, "figure_4", "plier.env")
  if (!file.exists(path)) {
    stop("figure_4/plier.env is missing from ", root, call. = FALSE)
  }
  for (line in readLines(path, warn = FALSE)) {
    line <- trimws(line)
    if (!startsWith(line, "export PLIER_")) next
    body <- sub("^export ", "", line)
    name <- sub("=.*$", "", body)
    if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", name)) next
    if (nzchar(Sys.getenv(name, unset = ""))) next          # already set: leave it
    value <- sub("^[^=]*=", "", body)
    value <- sub('^"', "", sub('"$', "", value))
    # "${NAME:-default}" -> default
    value <- sub("^\\$\\{[A-Za-z0-9_]+:-", "", value)
    value <- sub("\\}$", "", value)
    value <- gsub("\\$\\{?LANDSCAPE_ROOT\\}?", root, value)
    if (grepl("[$`]", value)) next
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(NULL)
}

# Before plier_helpers.R, not after: that file reads figure_4/plier_lvs.env on
# load with the same already-set-wins rule, so loading this one first is what
# keeps the two files' names from being able to collide.
load_plier_env()

source(file.path(here, "plier_helpers.R"))

# Progress and outcome go to the console. These rows were once written to a TSV
# that a shell stage folded into its own accounting; there is no stage here, and
# a run's report is what the reader sees on screen.
record <- function(status, check, detail) {
  message(sprintf("[%s] %s — %s", status, check, detail))
}

# plier_env_int() and plier_env_num() name figure_4/plier_lvs.env when a variable
# is unset, because lib/plier_helpers.R is shared with the panels and that is
# the file they read. The five below are figure_4/plier.env's.
ARMS <- trimws(strsplit(Sys.getenv("PLIER_ARMS", "rna,metab"), ",", fixed = TRUE)[[1]])
FORCE <- toupper(Sys.getenv("PLIER_FORCE", "FALSE")) == "TRUE"
K_HALF <- plier_env_int("PLIER_K_HALF", "figure_4/plier.env")
FRAC <- plier_env_num("PLIER_FRAC", "figure_4/plier.env")
SEED <- plier_env_int("PLIER_SEED", "figure_4/plier.env")
METAB_NUM_PC <- toupper(Sys.getenv("PLIER_METAB_USE_NUM_PC", "TRUE")) == "TRUE"
OUT_DIR <- plier_dir()
MANIFEST <- Sys.getenv("PLIER_MANIFEST_TSV", file.path(OUT_DIR, "plier_manifest.tsv"))
VERSION <- Sys.getenv("PLIER_VERSION", "2.0")

# ---- the report sub-command -------------------------------------------------

# Out of band: what is in a PLIER fit, one row per latent variable. No panel
# reads it and nothing below it runs when it is asked for.
#
# For each latent variable it prints:
#   cor_AB / cor_AM / cor_BM   its agreement between each pair of tissues — high
#                              in all three is what makes an LV "cross-tissue"
#   n_sig                      response cells with adjusted p < 0.05, of 22
#   markers                    features above the arm's affinity cutoff
#   top pathways               what the prior says it is
#
# The selected LVs are marked with *.
lv_report <- function(arm, n_pathways = 3) {
  arm <- plier_check_arm(arm)

  fit <- plier_fit(arm)
  input <- plier_input(arm)
  response <- plier_response(arm)
  selected <- plier_lvs(arm)
  sd_cutoff <- plier_marker_sd(arm)

  B <- fit$B
  rownames(B) <- plier_lv_names(nrow(B))
  correlation <- plier_tissue_correlation(B, input$metadata)

  scored <- if (arm == "rna") scale(fit$Z) else t(scale(t(fit$Z)))
  marker_counts <- colSums(scored > sd_cutoff)

  significant <- rowSums(response$neglog10_padj > -log10(0.05), na.rm = TRUE)

  cat(sprintf("%s fit: %d latent variables, %d features x %d samples\n",
              arm, ncol(fit$Z), nrow(input$matrix), ncol(input$matrix)))
  cat(sprintf("selected in figure_4/plier_lvs.env: %s\n\n",
              paste0("LV", selected, collapse = ", ")))
  cat(sprintf("%-8s %6s %6s %6s %6s %8s  %s\n",
              "LV", "cor_AB", "cor_AM", "cor_BM", "n_sig", "markers", "top pathways"))

  order_by <- order(-correlation$Blood_v_Muscle)
  for (i in order_by) {
    lv <- rownames(B)[i]
    top <- rownames(fit$Uauc)[order(-fit$Uauc[, i])][seq_len(n_pathways)]
    top <- top[fit$Uauc[top, i] > 0]
    cat(sprintf("%-8s %6.2f %6.2f %6.2f %6d %8d  %s\n",
                paste0(lv, if (i %in% selected) "*" else ""),
                correlation$Adipose_v_Blood[i], correlation$Adipose_v_Muscle[i],
                correlation$Blood_v_Muscle[i], significant[i], marker_counts[i],
                paste(plier_pathway_label(top), collapse = "; ")))
  }
  invisible(NULL)
}

# Dispatched here, above everything the fit defines, so the report path
# provably refits nothing.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && identical(args[1], "report")) {
  lv_report(
    arm = if (length(args) >= 2L) args[2] else "rna",
    n_pathways = if (length(args) >= 3L) as.integer(args[3]) else 3
  )
  quit(save = "no", status = 0)
}
if (length(args) > 0L) {
  stop("unknown argument: ", paste(args, collapse = " "),
       "\n  This script takes no arguments to fit, or `report <arm> [n_pathways]`.",
       call. = FALSE)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(MANIFEST), showWarnings = FALSE, recursive = TRUE)

message(sprintf("Fitting into %s — version: v%s", OUT_DIR, VERSION))
message(sprintf("  arms=%s  k=2x%d  frac=%s  seed=%d",
                paste(ARMS, collapse = ","), K_HALF, format(FRAC), SEED))

# ---- the nine gene-set databases, as one prior matrix -----------------------

# The legacy built one 0/1 matrix per database over the union of all nine gene
# lists and cbind()ed them, in 90 lines of nine near-identical blocks. Same
# matrix, one loop.
GENE_SET_DATABASES <- c("BIOCARTA", "KEGG_MEDICUS", "PID", "REACTOME", "WP",
                        "GOBP", "GOCC", "GOMF", "MITOCARTA")

signature_matrix <- function(databases) {
  sets <- unlist(
    lapply(databases, function(db) MOLECULAR_SIGNATURES[[db]]),
    recursive = FALSE
  )
  members <- sort(unique(unlist(sets, use.names = FALSE)))
  out <- matrix(0L, nrow = length(members), ncol = length(sets),
                dimnames = list(members, names(sets)))
  for (j in seq_along(sets)) {
    out[sets[[j]], j] <- 1L
  }
  out
}

#' Lift a feature x set matrix onto the rows of a data matrix.
#'
#' The legacy did this in a 15,000-iteration loop, each iteration scanning the
#' whole 200,000-row feature-to-gene table for one feature. Same result from one
#' match().
#'
#' Sets with fewer than `min_genes` of this matrix's features are dropped here
#' rather than inside PLIER. It is the same filter — PLIER applies `minGenes`
#' over the genes the data and the prior share, which is every row of this matrix
#' by construction — and the fit is identical either way. Doing it first is what
#' keeps the RNA arm inside a few GB: 14,200 gene sets over 14,428 genes is a
#' 1.6 GB matrix once PLIER makes it double, and 5,718 of those sets are about to
#' be discarded.
prior_for <- function(feature_keys, signatures, min_genes = 10) {
  out <- matrix(0L, nrow = length(feature_keys), ncol = ncol(signatures),
                dimnames = list(names(feature_keys), colnames(signatures)))
  hit <- match(feature_keys, rownames(signatures))
  found <- !is.na(hit)
  out[found, ] <- signatures[hit[found], ]
  n_sets_before <- ncol(out)
  out <- out[, colSums(out) >= min_genes, drop = FALSE]
  list(prior = out, n_mapped = sum(found),
       n_sets_before = n_sets_before, n_sets = ncol(out))
}

# ---- the RNA arm's input ----------------------------------------------------

# The three transcriptomics matrices, each restricted to the adult baseline
# visit, z-scored per tissue over the genes all three measure, and concatenated.
rna_input <- function() {
  # Muscle, adipose, blood is the order the legacy cbind()ed the three matrices
  # in, and so the column order its fit was produced in. fit_arm() fits on that
  # order; see the note there.
  qc <- list(
    Muscle  = MotrpacHumanPreSuspensionData::MUSCLE_TRNSCRPT_QC,
    Adipose = MotrpacHumanPreSuspensionData::ADIPOSE_TRNSCRPT_QC,
    Blood   = list(
      qc_norm = cbind(MotrpacHumanPreSuspensionData::blood_transcript_1$qc_norm,
                      MotrpacHumanPreSuspensionData::blood_transcript_2$qc_norm),
      sample_metadata = MotrpacHumanPreSuspensionData::blood_transcript_1$sample_metadata
    )
  )

  norms <- list()
  metas <- list()
  for (tissue in names(qc)) {
    meta <- qc[[tissue]]$sample_metadata
    rownames(meta) <- meta$vialLabel
    meta <- meta[meta$visitcode == "ADU_BAS", ]
    mat <- qc[[tissue]]$qc_norm[, rownames(meta), drop = FALSE]
    prefix <- paste0(tolower(tissue), "_")
    colnames(mat) <- paste0(prefix, colnames(mat))
    rownames(meta) <- paste0(prefix, rownames(meta))
    norms[[tissue]] <- mat
    metas[[tissue]] <- meta[, c("Tissue", "Sex", "Timepoint", "randomGroupCode",
                                "calculatedAge", "BMI", "pid", "visitcode")]
  }

  shared <- Reduce(intersect, lapply(norms, rownames))
  z <- do.call(cbind, unname(lapply(norms, function(m) {
    t(scale(t(m[shared, , drop = FALSE])))
  })))
  metadata <- do.call(rbind, unname(metas))

  metadata$Tissue <- plier_relabel(metadata$Tissue, PLIER_QC_TISSUE_LABELS, "tissue")
  metadata$Timepoint <- plier_relabel(metadata$Timepoint, PLIER_TIMEPOINT_LABELS, "timepoint")
  metadata$Modality <- plier_relabel(metadata$randomGroupCode, PLIER_MODALITY_LABELS, "modality")
  metadata$randomGroupCode <- NULL

  # The order the matrix was assembled in, before the sort below puts it in
  # display order. PLIER is fitted on the first and every panel draws the second.
  fit_order <- colnames(z)

  metadata <- metadata[order(
    match(metadata$Tissue, PLIER_TISSUE_ORDER),
    match(metadata$Modality, PLIER_SAMPLE_ORDER_MODALITY),
    match(metadata$Timepoint, PLIER_TIMEPOINT_ORDER),
    metadata$Sex
  ), ]
  z <- z[, rownames(metadata), drop = FALSE]

  # The prior is keyed on gene symbol, the matrix on versioned Ensembl id.
  f2g <- as.data.frame(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE)
  f2g <- f2g[f2g$assay == "transcript-rna-seq", ]
  symbol <- setNames(as.character(f2g$gene_symbol), as.character(f2g$feature_id))
  keys <- symbol[rownames(z)]
  names(keys) <- rownames(z)

  list(matrix = z, metadata = metadata, prior_keys = keys, fit_order = fit_order)
}

# ---- the metabolomics arm's input -------------------------------------------

# Each tissue's metabolomics platforms are combined and z-scored on their own,
# then the tissues are joined on the metabolites all three measure. Per-tissue
# rather than one load_qc() over everything, so that a tissue measuring a
# timepoint the others do not costs nothing.
metab_input <- function() {
  per_tissue <- list()
  for (tissue in MotrpacHumanPreSuspensionAnalysis::tissue_available_list()) {
    qc <- MotrpacHumanPreSuspensionData::load_qc(selected_tissues = tissue,
                                                 selected_omes = "metab")
    combined <- MotrpacHumanPreSuspensionData::combine_qc_matrixes(qc)
    rownames(combined) <- sub("\\..*", "", rownames(combined))
    colnames(combined) <- paste(colnames(combined), tissue, sep = "-")
    per_tissue[[tissue]] <- t(scale(t(combined)))
  }

  shared <- Reduce(intersect, lapply(per_tissue, rownames))
  z <- do.call(cbind, lapply(per_tissue, function(m) m[shared, , drop = FALSE]))

  pheno <- MotrpacHumanPreSuspensionData::load_pheno()$data
  pid_tp <- sub("-.*", "", colnames(z))
  meta <- pheno |>
    mutate(pid_tp = paste(.data$pid, .data$Timepoint, sep = "..")) |>
    filter(.data$pid_tp %in% pid_tp) |>
    distinct(.data$pid_tp, .keep_all = TRUE) |>
    select("pid", "Timepoint", "Sex", "randomGroupCode", "pid_tp") |>
    as.data.frame()
  rownames(meta) <- meta$pid_tp
  metadata <- meta[pid_tp, ]
  metadata$Tissue <- plier_relabel(sub(".*-", "", colnames(z)),
                                   PLIER_TISSUE_LABELS, "tissue")
  rownames(metadata) <- colnames(z)
  metadata$pid_tp <- NULL

  metadata$Timepoint <- plier_relabel(metadata$Timepoint, PLIER_TIMEPOINT_LABELS, "timepoint")
  metadata$Modality <- plier_relabel(metadata$randomGroupCode, PLIER_MODALITY_LABELS, "modality")
  metadata$randomGroupCode <- NULL

  # The order the matrix was assembled in, before the sort below puts it in
  # display order. PLIER is fitted on the first and every panel draws the second.
  fit_order <- colnames(z)

  metadata <- metadata[order(
    match(metadata$Tissue, PLIER_TISSUE_ORDER),
    match(metadata$Modality, PLIER_SAMPLE_ORDER_MODALITY),
    match(metadata$Timepoint, PLIER_TIMEPOINT_ORDER),
    metadata$Sex
  ), ]
  z <- z[, rownames(metadata), drop = FALSE]

  keys <- rownames(z)
  names(keys) <- rownames(z)
  list(matrix = z, metadata = metadata, prior_keys = keys, fit_order = fit_order)
}

# ---- one arm ----------------------------------------------------------------

fit_arm <- function(arm) {
  paths <- c(plier = plier_path(arm, "plier"),
             input = plier_path(arm, "input"),
             response = plier_path(arm, "response"))

  if (!FORCE && all(file.exists(paths))) {
    record("SKIP", paste0("plier:", arm, ":reuse"),
           paste0("already fitted: ", basename(paths[["plier"]]),
                  " — set PLIER_FORCE=TRUE to refit"))
    return(invisible(NULL))
  }

  # A fit on disk with no response grid beside it: rebuild the grid off the fit
  # rather than refitting. The grid is the cheap half and it moves on its own —
  # it is Wilcoxon tests through ggpubr, and compare_means() changed how it
  # adjusts p at ggpubr 1.0.0 (config/required_packages.tsv). Deleting the
  # response .rds is how to redo it after such a change; the RNA arm is two
  # hours and none of it is the grid.
  if (!FORCE && all(file.exists(paths[c("plier", "input")]))) {
    message(sprintf("  %s — reusing the fit; rebuilding the response grid", arm))
    fit <- readRDS(paths[["plier"]])
    input <- readRDS(paths[["input"]])
    response <- plier_lv_response_stats(fit, input$metadata)
    saveRDS(response, paths[["response"]])
    record("PASS", paste0("plier:", arm, ":response"),
           sprintf("%d LVs x %d comparisons rebuilt from %s (ggpubr %s)",
                   nrow(response$effect), nrow(response$comparisons),
                   basename(paths[["plier"]]),
                   as.character(utils::packageVersion("ggpubr"))))
    return(invisible(NULL))
  }

  message(sprintf("  %s — assembling the input matrix", arm))
  input <- if (arm == "rna") rna_input() else metab_input()
  message(sprintf("        %d features x %d samples (%s)",
                  nrow(input$matrix), ncol(input$matrix),
                  paste(sprintf("%s %d", names(table(input$metadata$Tissue)),
                                as.integer(table(input$metadata$Tissue))),
                        collapse = ", ")))

  message(sprintf("  %s — assembling the prior", arm))
  signatures <- if (arm == "rna") {
    signature_matrix(GENE_SET_DATABASES)
  } else {
    signature_matrix("REFMET")
  }
  prior <- prior_for(input$prior_keys, signatures)
  message(sprintf("        %d of %d sets keep at least 10 of these features; %d of %d features map into at least one",
                  prior$n_sets, prior$n_sets_before, prior$n_mapped, nrow(input$matrix)))
  if (prior$n_mapped == 0) {
    stop("no feature of the ", arm, " matrix maps into the prior", call. = FALSE)
  }
  # The signature matrix is the largest object here after the prior itself and
  # nothing below reads it. PLIER is about to allocate several copies of the
  # prior, so it is dropped rather than held for the duration of the fit.
  rm(signatures)
  gc(verbose = FALSE)

  # PLIER sees the matrix in assembly order, not in the display order the
  # metadata carries. The fit is not invariant to column order: refitting the
  # sorted matrix converges elsewhere, and the latent variables
  # figure_4/plier_lvs.env names are renumbered — seven metabolomics LVs whose
  # markers then do not overlap the legacy's at all. The legacy sorted only after
  # fitting, so assembly order is the order both published fits were produced in.
  fit_matrix <- as.matrix(input$matrix[, input$fit_order, drop = FALSE])

  k_half <- K_HALF
  if (arm == "metab" && METAB_NUM_PC) {
    estimated <- PLIER::num.pc(as.data.frame(fit_matrix), seed = SEED)
    k_half <- min(estimated, K_HALF)
    message(sprintf("        num.pc() = %d, capped at %d -> k = %d",
                    estimated, K_HALF, 2 * k_half))
  }

  message(sprintf("  %s — PLIER(k = %d, frac = %s, seed = %d); this is the slow part",
                  arm, 2 * k_half, format(FRAC), SEED))
  started <- Sys.time()
  fit <- PLIER::PLIER(fit_matrix, prior$prior,
                      k = 2 * k_half, trace = FALSE, frac = FRAC,
                      scale = TRUE, seed = SEED)
  elapsed <- difftime(Sys.time(), started, units = "mins")
  message(sprintf("        %d latent variables in %.1f min",
                  ncol(fit$Z), as.numeric(elapsed)))

  # Back to display order: B and residual are the only sample-columned pieces of
  # a fit. plier_lv_response_stats() asserts colnames(B) == rownames(metadata),
  # plier_tissue_correlation() splits B by tissue positionally, and ED6A plots
  # sample index. The legacy subscripted B by column name and never reordered it.
  fit$B <- fit$B[, colnames(input$matrix), drop = FALSE]
  fit$residual <- fit$residual[, colnames(input$matrix), drop = FALSE]

  saveRDS(fit, paths[["plier"]])
  saveRDS(input[c("matrix", "metadata")], paths[["input"]])

  message(sprintf("  %s — the LV x comparison response grid", arm))
  response <- plier_lv_response_stats(fit, input$metadata)
  saveRDS(response, paths[["response"]])

  record("PASS", paste0("plier:", arm),
         sprintf("%d LVs over %d features x %d samples, %d prior sets, %d comparisons, %.1f min",
                 ncol(fit$Z), nrow(input$matrix), ncol(input$matrix),
                 prior$n_sets, nrow(response$comparisons),
                 as.numeric(elapsed)))

  manifest_row <- data.frame(
    arm = arm,
    version = VERSION,
    file = basename(paths[["plier"]]),
    n_features = nrow(input$matrix),
    n_samples = ncol(input$matrix),
    n_prior_sets = prior$n_sets,
    k = ncol(fit$Z),
    frac = FRAC,
    seed = SEED,
    fitted_minutes = round(as.numeric(elapsed), 1),
    fitted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )
  write.table(
    manifest_row, MANIFEST, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = !file.exists(MANIFEST), append = file.exists(MANIFEST)
  )
  invisible(NULL)
}

# ---- the arms ---------------------------------------------------------------

# An arm that fails costs its own fit and not the other's, and the run exits
# non-zero, which is how a non-interactive caller learns about a failure it did
# not read. The shell stage this replaces did the same through its check counts.
failed <- character(0)
for (arm in ARMS) {
  arm <- plier_check_arm(arm)
  tryCatch(fit_arm(arm), error = function(e) {
    record("FAIL", paste0("plier:", arm), conditionMessage(e))
    failed <<- c(failed, arm)
  })
}

if (length(failed) > 0) {
  record("FAIL", "plier:fit",
         paste0(length(failed), " arm(s) failed: ", paste(failed, collapse = ", ")))
  quit(save = "no", status = 1)
}
