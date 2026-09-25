# highlights.R — reader for config/highlights.json: the curated pathway
# selections, single-feature plots, gene lists, named features and clinical x
# omics features the panel and table scripts draw.
#
# Sourced after panel_export.R, which defines landscape_root().

.highlights_cache <- new.env(parent = emptyenv())

highlights <- function() {
  if (is.null(.highlights_cache$doc)) {
    path <- file.path(landscape_root(), "config", "highlights.json")
    if (!file.exists(path)) {
      stop("config/highlights.json not found at ", path, call. = FALSE)
    }
    .highlights_cache$doc <- jsonlite::read_json(path, simplifyVector = FALSE)
  }
  .highlights_cache$doc
}

.highlights_entry <- function(section, key) {
  entries <- highlights()[[section]]
  if (!key %in% names(entries)) {
    stop("config/highlights.json has no ", section, " entry '", key, "'. Known: ",
         paste(names(entries), collapse = ", "), call. = FALSE)
  }
  entries[[key]]
}

#' One pathway selection as a data frame, one row per set, in file order.
#'
#' Columns are the fields the entries carry (set, set_id, label).
highlight_pathways <- function(key) {
  entry <- .highlights_entry("pathways", key)
  fields <- unique(unlist(lapply(entry$sets, names)))
  out <- as.data.frame(
    lapply(stats::setNames(fields, fields), function(field) {
      vapply(entry$sets, function(s) {
        if (is.null(s[[field]])) NA_character_ else as.character(s[[field]])
      }, character(1))
    }),
    stringsAsFactors = FALSE
  )
  if (!entry$match_on %in% names(out) || anyNA(out[[entry$match_on]])) {
    stop("pathway selection '", key, "' has a set with no ", entry$match_on,
         call. = FALSE)
  }
  out
}

#' One single-feature plot's feature: feature, tissues, omes, and title if any.
highlight_feature <- function(panel, name) {
  plots <- .highlights_entry("features", panel)
  if (!name %in% names(plots)) {
    stop("config/highlights.json has no feature '", name, "' for ", panel,
         ". Known: ", paste(names(plots), collapse = ", "), call. = FALSE)
  }
  entry <- plots[[name]]
  list(
    feature = entry$feature,
    tissues = unlist(entry$tissues),
    omes = unlist(entry$omes),
    title = entry$title
  )
}

#' The named features for ST1e: feature_id, tissue, assay, notes.
highlight_named_features <- function() {
  features <- highlights()$named_features$features
  field <- function(name) {
    vapply(features, function(f) {
      if (is.null(f[[name]])) "" else as.character(f[[name]])
    }, character(1))
  }
  data.frame(feature_id = field("feature_id"), tissue = field("tissue"),
             assay = field("assay"), notes = field("notes"),
             stringsAsFactors = FALSE)
}

#' One hand-picked gene list, as a character vector in file order.
highlight_gene_list <- function(key) {
  entry <- .highlights_entry("gene_lists", key)
  genes <- unlist(entry$genes)
  if (length(genes) == 0) {
    stop("gene list '", key, "' in config/highlights.json is empty", call. = FALSE)
  }
  as.character(genes)
}

#' The features one clinical x omics selection names, one row each, in file
#' order: trait, tissue, assay, feature_id, gene, timepoint, notes.
highlight_clinical_omics_features <- function(table) {
  entry <- .highlights_entry("clinical_omics_features", table)
  fields <- c("trait", "tissue", "assay", "feature_id", "gene", "timepoint", "notes")
  out <- as.data.frame(
    lapply(stats::setNames(fields, fields), function(field) {
      vapply(entry$features, function(f) {
        if (is.null(f[[field]])) NA_character_ else as.character(f[[field]])
      }, character(1))
    }),
    stringsAsFactors = FALSE
  )
  keys <- c("trait", "tissue", "assay", "feature_id", "timepoint")
  if (anyNA(out[keys])) {
    stop("clinical x omics selection '", table, "' has a feature missing one of ",
         paste(keys, collapse = ", "), call. = FALSE)
  }
  out
}
