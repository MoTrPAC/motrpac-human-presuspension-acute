#!/usr/bin/env Rscript
# Supplementary Table 2 — the standalone sub-tables
#
# Sub-tables:  ST2a  covariates included in the primary differential analysis
#              ST2d  DA percent overlap
#
# The rest of ST2 is a panel's own numbers and is written by that panel's
# figure script, which config/table_map.json records:
#   ST2b  total features detected            FIG1.R  (FIG1C)
#   ST2c  number of DA features              FIG2.R  (FIG2A)
#   ST2e  cross-tissue ORA                   FIG2.R  (FIG2B)
#   ST2f  phospho ORA (PTMsigDB)             ED2.R   (ED2E)
#   ST2g  sex-specific DA estimate correlation   ED3.R (ED3A)
#   ST2h  sex-specific differences enrichment    ED3.R (ED3A)
#
# ST2d needs consortium data access. load_differential_analysis(epigen = TRUE)
# downloads the epigenomics DA tables into EPIGEN_QC_DIR on first use, then
# reuses the cache. ST2a reads a stored object and needs neither.
#
#   Rscript figures/landscape/tables/ST2.R          every sub-table here
#   Rscript figures/landscape/tables/ST2.R ST2a     one of them

suppressPackageStartupMessages({
  library(dplyr)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R for epigen_qc_dir(); table_export.R after it.
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))

# ---- shared ----------------------------------------------------------------

TISSUES <- c("adipose", "blood", "muscle")

# The two platforms whose feature ids are replaced by their UniProt accession.
UNIPROT_KEYED <- c("prot-pr", "prot-ol")

# The label each DA assay is reported under, in the order the rows are written.
# Named locally rather than taken from ASSAY_DISPLAY_ORDER in
# lib/supp_table_helpers.R: this table pools prot-pr and prot-ol into one "prot"
# label, which is not an assay that order knows.
ASSAY_GROUPS <- list(
  "metab" = "metab",
  "transcript-rna-seq" = "transcript-rna-seq",
  "prot" = UNIPROT_KEYED,
  "prot-ph" = "prot-ph",
  "epigen-methylcap-seq" = "epigen-methylcap-seq",
  "epigen-atac-seq" = "epigen-atac-seq"
)

ASSAY_LABEL <- stats::setNames(
  rep(names(ASSAY_GROUPS), lengths(ASSAY_GROUPS)),
  unlist(ASSAY_GROUPS, use.names = FALSE)
)

assay_label <- function(assay) {
  out <- unname(ASSAY_LABEL[as.character(assay)])
  unmapped <- unique(as.character(assay)[is.na(out)])
  if (length(unmapped) > 0) {
    stop("ST2d: assay(s) with no label: ", paste(sort(unmapped), collapse = ", "),
         "\n  Add them to ASSAY_GROUPS in this script.", call. = FALSE)
  }
  out
}

# Only ST2d reads this, so it stays lazy: `Rscript tables/ST2.R ST2a` must not
# download the epigenomics DA tables to write a stored object.
all_da <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
        selected_omes = "all",
        selected_tissues = "all",
        epigen = TRUE,
        repo_local_dir = epigen_qc_dir(),
        single_matrix = TRUE
      )
    }
    cache
  }
})

# ---- ST2a — covariates included in the primary differential analysis -------

# Straight off COVARIATES_FILE. The object already is this table: one row per
# ome x covariate, carrying whether the covariate is technical or a design
# parameter, whether it is coded as a factor or numerically, and which tissue it
# applies to. Nothing is counted, joined or reshaped here, so there is nothing
# to get wrong beyond the column contract, which is checked.
#
# Row order is the object's own. It is a stored data frame, so that order is
# stable across runs, and it groups the rows by covariate, which is how the
# table reads.
st2a <- function() {
  spec <- table_init("ST2a")

  out <- COVARIATES_FILE

  missing <- setdiff(spec$columns, names(out))
  if (length(missing) > 0) {
    stop("ST2a: COVARIATES_FILE has no column(s) ",
         paste(missing, collapse = ", "),
         "\n  It carries: ", paste(names(out), collapse = ", "),
         "\n  Reconcile `columns` for ST2a in config/table_map.json with the object.",
         call. = FALSE)
  }

  # Selecting rather than reordering in place: the manifest decides which
  # columns ship and in what order, so a column added to the object upstream
  # does not silently appear in the table.
  out <- out[, spec$columns, drop = FALSE]

  # Factors would export as their labels anyway, but as.character makes that
  # explicit rather than incidental, and keeps export_table()'s cell check from
  # having to reason about levels.
  out[] <- lapply(out, as.character)

  export_table(out, "ST2a")
}

# ---- ST2d — DA percent overlap ---------------------------------------------

# For each assay and each ORDERED pair of tissues: how many of the features DA
# in tissue_query are also DA in tissue_compare. Directional, so the two
# directions of a pair share n_overlap and differ in n_query and percent.
#
# DA is FDR < 0.05 in any exercise-with-controls contrast, either exercise arm,
# any timepoint.
#
# prot-pr and prot-ol are reported together as "prot" AND keyed on UniProt, so a
# protein measured by mass spec in adipose or muscle and on Olink in blood is
# one feature. Olink ids are OID accessions and global-proteomics ids are
# UniProt accessions, so without this the proteomics overlap would be zero by
# construction rather than by measurement. Every other assay is keyed on
# feature_id as stored.
#
# The grid is complete: every assay against every ordered tissue pair, whether
# or not that assay was run in that tissue. An assay with no DA feature in
# tissue_query is n_query = 0 and a blank percent.
st2d <- function() {
  spec <- table_init("ST2d")

  # as.data.frame() at the boundary: load_differential_analysis() returns a
  # data.table, and [.data.table evaluates its i expression inside the table, so
  # a base subset by a variable that shares a column's name matches every row
  # instead of the intended ones.
  significant <- all_da() %>%
    filter(.data$contrast_type == "exercise_with_controls",
           .data$adj_p_value < 0.05) %>%
    distinct(.data$assay, .data$feature_id, .data$tissue) %>%
    mutate(raw_assay = as.character(.data$assay),
           assay = assay_label(.data$assay),
           feature_id = as.character(.data$feature_id),
           tissue = as.character(.data$tissue)) %>%
    as.data.frame()

  unplaced_tissue <- setdiff(unique(significant$tissue), TISSUES)
  if (length(unplaced_tissue) > 0) {
    stop("ST2d: tissue(s) the grid does not carry: ",
         paste(sort(unplaced_tissue), collapse = ", "),
         "\n  Add them to TISSUES in this script.", call. = FALSE)
  }

  # ---- UniProt keying for the two proteomics platforms ---------------------

  bridge <- MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE %>%
    as.data.frame() %>%
    filter(.data$assay %in% UNIPROT_KEYED) %>%
    distinct(.data$assay, .data$feature_id, .data$uniprot) %>%
    transmute(raw_assay = as.character(.data$assay),
              feature_id = as.character(.data$feature_id),
              mapped_uniprot = as.character(.data$uniprot))

  # A one-to-many mapping would count one protein once per accession.
  ambiguous <- bridge %>%
    count(.data$raw_assay, .data$feature_id, name = "n_uniprot") %>%
    filter(.data$n_uniprot > 1)
  if (nrow(ambiguous) > 0) {
    stop("ST2d: ", nrow(ambiguous), " proteomics feature(s) map to more than one ",
         "UniProt accession, e.g. ", ambiguous$feature_id[1],
         "\n  The overlap would count them once per accession.", call. = FALSE)
  }

  significant <- left_join(significant, bridge, by = c("raw_assay", "feature_id"))

  unmapped_prot <- significant[significant$raw_assay %in% UNIPROT_KEYED &
                                 is.na(significant$mapped_uniprot), ]
  if (nrow(unmapped_prot) > 0) {
    stop("ST2d: ", nrow(unmapped_prot), " proteomics feature(s) have no UniProt ",
         "accession in HUMAN_FEATURE_TO_GENE, e.g. ",
         unmapped_prot$raw_assay[1], " / ", unmapped_prot$feature_id[1],
         "\n  They cannot be matched across platforms.", call. = FALSE)
  }

  # Distinct again: the two platforms collapse onto one identifier here, so a
  # protein measured on both in one tissue is one feature rather than two.
  significant <- significant %>%
    mutate(feature_id = ifelse(.data$raw_assay %in% UNIPROT_KEYED,
                               .data$mapped_uniprot, .data$feature_id)) %>%
    distinct(.data$assay, .data$feature_id, .data$tissue)

  # ---- the DA feature set of every assay x tissue --------------------------

  da_features <- lapply(names(ASSAY_GROUPS), function(label) {
    rows <- significant$assay == label
    fid <- significant$feature_id[rows]
    tis <- significant$tissue[rows]
    stats::setNames(lapply(TISSUES, function(t) unique(fid[tis == t])), TISSUES)
  })
  names(da_features) <- names(ASSAY_GROUPS)

  # ---- the grid ------------------------------------------------------------

  # Assay outer, then tissue_query, then tissue_compare, each in the order
  # stated above. A tissue is never compared against itself.
  grid <- do.call(rbind, lapply(names(ASSAY_GROUPS), function(label) {
    do.call(rbind, lapply(TISSUES, function(query) {
      data.frame(
        assay = label,
        tissue_query = query,
        tissue_compare = setdiff(TISSUES, query),
        stringsAsFactors = FALSE
      )
    }))
  }))

  grid$n_query <- vapply(seq_len(nrow(grid)), function(i) {
    length(da_features[[grid$assay[i]]][[grid$tissue_query[i]]])
  }, integer(1))

  grid$n_overlap <- vapply(seq_len(nrow(grid)), function(i) {
    sets <- da_features[[grid$assay[i]]]
    length(intersect(sets[[grid$tissue_query[i]]], sets[[grid$tissue_compare[i]]]))
  }, integer(1))

  # Undefined rather than zero where the query tissue has no DA feature at all,
  # which covers both "not assayed here" and "assayed, nothing significant".
  grid$percent_overlap <- ifelse(grid$n_query > 0,
                                 100 * grid$n_overlap / grid$n_query,
                                 NA_real_)

  # ---- invariants ----------------------------------------------------------

  # An intersection is symmetric, so the two directions of a pair must report
  # the same n_overlap even though their percentages differ. This is what
  # catches a set built for the wrong assay or the wrong tissue; the percentages
  # alone would look plausible.
  key <- paste(grid$assay, grid$tissue_query, grid$tissue_compare)
  mirror <- match(paste(grid$assay, grid$tissue_compare, grid$tissue_query), key)
  asymmetric <- which(grid$n_overlap != grid$n_overlap[mirror])
  if (length(asymmetric) > 0) {
    i <- asymmetric[1]
    stop("ST2d: n_overlap is not symmetric for ", grid$assay[i], " ",
         grid$tissue_query[i], "/", grid$tissue_compare[i], ": ",
         grid$n_overlap[i], " vs ", grid$n_overlap[mirror[i]], call. = FALSE)
  }

  too_many <- which(grid$n_overlap > grid$n_query)
  if (length(too_many) > 0) {
    i <- too_many[1]
    stop("ST2d: n_overlap exceeds n_query for ", grid$assay[i], " ",
         grid$tissue_query[i], " -> ", grid$tissue_compare[i], ": ",
         grid$n_overlap[i], " > ", grid$n_query[i], call. = FALSE)
  }

  export_table(grid[, spec$columns, drop = FALSE], "ST2d")
}

run_tables(list(ST2a = st2a, ST2d = st2d))
