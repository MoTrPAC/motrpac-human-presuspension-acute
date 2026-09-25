# table_export.R — the export contract every supplementary-table script uses.
#
# A table script never names its own output file. It calls table_init() to learn
# which sub-table it is building, and export_table() to write it. Both read
# config/table_map.json, so renumbering a sub-table (ST2d -> ST2c, say) is an
# edit to that file alone and the exported TSV name follows.
#
# This is the twin of panel_export.R, and deliberately a separate file. A
# sub-table has no seed, no page size and no graphics device; a panel has no
# column schema. One helper covering both would be half-inapplicable in every
# call.
#
# What the write settings are for: these files are opened in Excel, and read
# back by scripts that expect the same conventions as the rest of MoTrPAC.
#   sep = "\t"        most MoTrPAC files are tab-delimited, and a pathway or
#                     participant label containing a comma is ordinary.
#   na = ""           an empty cell reads as missing in every spreadsheet
#                     program. The literal string "NA" reads as text, and in a
#                     numeric column as a value.
#   quote = FALSE     with a tab separator nothing needs quoting, and unquoted
#                     text survives the round trip through Excel unchanged.
#   row.names = FALSE a leading unnamed column shifts every header one to the
#                     left on read-back, which is the single most common way a
#                     delimited table is silently corrupted.
#
# Usage in a standalone table script (tables/ST1.R and its siblings):
#
#   here <- dirname(sub("^--file=", "",
#                       grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
#   source(file.path(here, "..", "lib", "table_export.R"))
#   spec <- table_init("ST1c")             # announces the output path
#   ...build df...
#   export_table(df, "ST1c")               # writes it, records the export
#
# One script per supplementary table NUMBER, so a table script usually writes
# several sub-tables off one data load. run_tables() runs them the way
# run_panels() runs a figure's panels, and takes the same selection:
#
#   run_tables(list(ST3a = st3a, ST3b = st3b, ST3g = st3g))
#
#   Rscript figures/landscape/tables/ST3.R            # every sub-table
#   Rscript figures/landscape/tables/ST3.R ST3g       # one of them
#
# Usage in a FIGURE script whose panel also writes the table its own numbers
# are — nineteen of the twenty-nine sub-tables built here are written this way,
# from the object the panel was about to plot. No table_init(): panel_init() has
# already attached the data packages and announced the build, and export_table()
# resolves the manifest row itself.
#
#   source(file.path(here, "..", "lib", "panel_export.R"))
#   source(file.path(here, "..", "lib", "table_export.R"))
#   spec <- panel_init("FIG2A")
#   ...build all_da...
#   export_table(all_da_by_grp_tp, "ST2c")
#   export_panel(p, "FIG2A")
#
# Both helpers are sourced into one script now, so anything defined in both has
# to mean the same thing in both. `%||%` and installed_version() are identical
# by inspection. attach_data_packages() was NOT — the two differed in which
# manifest field they named in the error — so this file's is
# attach_table_data_packages() and nothing here shadows panel_export.R.
# landscape_root() and its two helpers are defined in both, identically, so a
# table script can source this file alone and source order changes nothing.

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- the landscape root -----------------------------------------------------

# Identical to panel_export.R's, so that either file can be sourced alone and
# neither shadows the other with different behaviour.

#' Absolute path to figures/landscape/.
#'
#' Resolved from LANDSCAPE_ROOT, else the directory of the running script walked
#' up, else the working directory. A directory is the root when it holds
#' lib/panel_export.R and config/panel_map.json.
landscape_root <- function() {
  from_env <- Sys.getenv("LANDSCAPE_ROOT", unset = "")
  if (nzchar(from_env)) {
    if (!is_landscape_root(from_env)) {
      stop(
        "LANDSCAPE_ROOT does not look like figures/landscape/: ", from_env,
        "\n  No lib/panel_export.R and config/panel_map.json there.",
        call. = FALSE
      )
    }
    return(normalizePath(from_env, mustWork = TRUE))
  }
  from_script <- landscape_root_from_script()
  if (!is.null(from_script)) {
    return(from_script)
  }
  if (is_landscape_root(getwd())) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  stop(
    "cannot tell where figures/landscape/ is.",
    "\n  Run the script with Rscript, or set LANDSCAPE_ROOT to the",
    " figures/landscape directory.",
    call. = FALSE
  )
}

#' TRUE when a directory holds the two files that make it the root.
is_landscape_root <- function(dir) {
  file.exists(file.path(dir, "lib", "panel_export.R")) &&
    file.exists(file.path(dir, "config", "panel_map.json"))
}

#' The root, walked up from the directory of the running script, or NULL.
landscape_root_from_script <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg) == 0) {
    return(NULL)
  }
  path <- sub("^--file=", "", arg[length(arg)])
  if (!nzchar(path)) {
    return(NULL)
  }
  dir <- tryCatch(normalizePath(dirname(path), mustWork = TRUE),
                  error = function(e) NULL)
  while (!is.null(dir)) {
    if (is_landscape_root(dir)) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) {
      return(NULL)
    }
    dir <- parent
  }
  NULL
}

.table_map_cache <- new.env(parent = emptyenv())

# ---- the manifest ----------------------------------------------------------

#' Read config/table_map.json, once per session.
table_map <- function() {
  if (!is.null(.table_map_cache$map)) {
    return(.table_map_cache$map)
  }
  path <- file.path(landscape_root(), "config", "table_map.json")
  if (!file.exists(path)) {
    stop("table manifest not found: ", path, call. = FALSE)
  }
  map <- jsonlite::read_json(path, simplifyVector = FALSE)
  .table_map_cache$map <- map
  map
}

#' Every sub-table id in the manifest, in manifest order.
#'
#' Includes the sub-tables this repo does not build. `built_here = TRUE` is what
#' marks the ones a script here writes; see table_ids_built().
table_ids <- function() {
  vapply(table_map()$tables, function(t) t$table, character(1))
}

#' The sub-table ids this repo builds, in manifest order.
table_ids_built <- function() {
  tabs <- table_map()$tables
  vapply(Filter(function(t) isTRUE(t$built_here), tabs),
         function(t) t$table, character(1))
}

#' The manifest entry for one sub-table, with supp-table fields resolved onto it.
#'
#' @param table Sub-table id, e.g. "ST1c". Must appear in config/table_map.json.
table_spec <- function(table) {
  map <- table_map()
  ids <- table_ids()
  if (!table %in% ids) {
    stop(
      "unknown sub-table '", table, "'. Known sub-tables: ",
      paste(ids, collapse = ", "), call. = FALSE
    )
  }
  spec <- map$tables[[match(table, ids)]]

  supp_ids <- vapply(map$supp_tables, function(s) s$supp_id, character(1))
  supp <- map$supp_tables[[match(spec$supp_id, supp_ids)]]
  spec$supp_dir <- supp$supp_dir
  spec$supp_title <- supp$title

  spec$columns <- unlist(spec$columns %||% list())
  spec
}

#' Absolute path of the TSV a sub-table is required to write.
table_output_path <- function(table) {
  spec <- table_spec(table)
  if (is.null(spec$output)) {
    stop("sub-table '", table, "' declares no output: it is not built here ",
         "(", spec$reason_not_built %||% "built_here = false", ")", call. = FALSE)
  }
  file.path(landscape_root(), spec$output)
}

# ---- start of a table script -----------------------------------------------

#' Attach the data packages a sub-table declares.
#'
#' Named apart from panel_export.R's so a panel script that sources both keeps
#' the panel version for its own init; see the header note.
attach_table_data_packages <- function(spec) {
  for (pkg in unlist(spec$data_packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("sub-table '", spec$table, "' needs ", pkg,
           ", which is not installed", call. = FALSE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, warn.conflicts = FALSE)
    )
  }
  invisible(NULL)
}

#' Announce the sub-table being built and attach what it reads.
#'
#' No seed, unlike panel_init(): nothing in the table arm is stochastic. A
#' sub-table that ever needs one should carry it in the manifest the way a panel
#' does, rather than calling set.seed() in the script.
table_init <- function(table) {
  spec <- table_spec(table)
  if (!isTRUE(spec$built_here)) {
    stop("sub-table '", table, "' is not built in this repo (",
         spec$reason_not_built %||% "built_here = false", ")", call. = FALSE)
  }
  attach_table_data_packages(spec)
  message(sprintf(
    "[%s] %s\n        %d columns  ->  %s",
    spec$table, spec$title, length(spec$columns), spec$output
  ))
  invisible(spec)
}

# ---- the export ------------------------------------------------------------

#' Write one sub-table to its declared path and record the export.
#'
#' @param df    A data frame. Its columns must match the manifest schema exactly,
#'              in order — see check_schema().
#' @param table Sub-table id, e.g. "ST1c".
export_table <- function(df, table) {
  spec <- table_spec(table)
  if (!is.data.frame(df)) {
    stop("export_table() needs a data frame, got ", class(df)[1], call. = FALSE)
  }
  check_schema(df, spec)
  check_cells(df, spec)

  out <- table_output_path(table)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

  defaults <- table_map()$defaults
  write.table(
    df, file = out,
    sep = defaults$sep %||% "\t",
    na = defaults$na %||% "",
    quote = isTRUE(defaults$quote),
    row.names = isTRUE(defaults$row_names),
    col.names = TRUE
  )

  message(sprintf("        wrote %d rows x %d cols  (%s)",
                  nrow(df), ncol(df), format_bytes(file.size(out))))
  record_table_export(spec, out, df)
  invisible(out)
}

#' The declared columns, in order, or a message saying exactly how they differ.
#'
#' A schema check rather than a column count: the failure this catches is an
#' upstream rename, where the count is unchanged and the header is wrong. Left
#' unchecked it ships a table whose header no longer describes its contents,
#' which is not visible in the file and not recoverable from it.
check_schema <- function(df, spec) {
  want <- spec$columns
  if (length(want) == 0) {
    stop("sub-table '", spec$table, "' declares no columns in table_map.json",
         call. = FALSE)
  }
  got <- names(df)
  if (identical(got, want)) {
    return(invisible(TRUE))
  }
  missing <- setdiff(want, got)
  extra <- setdiff(got, want)
  detail <- character(0)
  if (length(missing)) {
    detail <- c(detail, paste0("  missing: ", paste(missing, collapse = ", ")))
  }
  if (length(extra)) {
    detail <- c(detail, paste0("  unexpected: ", paste(extra, collapse = ", ")))
  }
  if (!length(detail)) {
    # Same names, wrong order. Worth its own message: the columns are all there,
    # so the diff above would be empty and the failure would read as a mystery.
    detail <- c(
      "  same columns, different order",
      paste0("  declared: ", paste(want, collapse = ", ")),
      paste0("  built:    ", paste(got, collapse = ", "))
    )
  }
  stop(
    "sub-table '", spec$table, "' does not match the schema in table_map.json\n",
    paste(detail, collapse = "\n"),
    "\n  Fix the script, or update `columns` for ", spec$table,
    " if the schema really changed.",
    call. = FALSE
  )
}

#' Refuse a cell that would not survive the round trip through a TSV.
#'
#' A tab or a newline inside a value silently splits a row, and the file stays
#' well-formed enough that nothing downstream notices. Checked on the way out so
#' the failure names the column and the value while the script that produced it
#' is still on the stack.
check_cells <- function(df, spec) {
  for (col in names(df)) {
    v <- df[[col]]
    if (!is.character(v) && !is.factor(v)) next
    v <- as.character(v)
    bad <- which(grepl("[\t\r\n]", v))
    if (length(bad)) {
      stop(
        "sub-table '", spec$table, "': column '", col, "' has a tab or newline ",
        "inside a value (row ", bad[1], "): ",
        encodeString(substr(v[bad[1]], 1, 60), quote = "'"),
        "\n  A tab-delimited file cannot carry it. Strip or replace it in the ",
        "script.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

format_bytes <- function(n) {
  if (is.na(n)) return("unknown size")
  units <- c("B", "KB", "MB", "GB")
  i <- max(1, min(length(units), floor(log(max(n, 1), 1024)) + 1))
  sprintf("%.1f %s", n / 1024^(i - 1), units[i])
}

installed_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)),
           error = function(e) NA_character_)
}

# ---- running a table script's sub-tables ------------------------------------

#' What a sub-table function wrote, for the run report.
#'
#' export_table() returns the path it wrote, so a function that ends with it
#' says where it went, and one that writes two sub-tables returns both if it
#' returns c(export_table(a, "ST6a"), export_table(b, "ST6b")). Anything else
#' falls back to what the manifest declares for the id.
table_result_path <- function(table, value) {
  if (is.character(value) && length(value) > 0 && all(nzchar(value))) {
    return(paste(value, collapse = "; "))
  }
  tryCatch(table_output_path(table), error = function(e) "")
}

#' Why a panel or sub-table failed, in one or two lines.
#'
#' The raw condition message is what R says went wrong; this says what to do
#' about it. Most failures here have one of a handful of causes, and the
#' difference between them is the difference between installing a package,
#' authenticating to a bucket, and reporting a real defect.
#'
#' Returns character(0) when nothing is recognised, so an unfamiliar error is
#' shown unadorned rather than guessed at.
#'
#' @param e A condition.
explain_failure <- function(e) {
  msg <- conditionMessage(e)
  hit <- function(p) grepl(p, msg, ignore.case = TRUE)

  if (hit("there is no package called")) {
    pkg <- sub(".*there is no package called [‘'\"]?([^’'\"]+).*", "\\1", msg)
    return(c(sprintf("why: the package '%s' is not installed.", pkg),
             "     config/panel_packages.txt lists every package the figures need."))
  }
  if (hit("could not find function")) {
    return(c("why: a function is not in scope, so a package is missing or a helper",
             "     was not sourced. config/panel_packages.txt lists the packages."))
  }
  if (hit("object '[A-Z][A-Z0-9_]+' not found")) {
    return(c("why: a data-package object is not on the search path.",
             "     Check MotrpacHumanPreSuspensionData and ...Analysis are installed",
             "     at the versions config/required_packages.tsv pins, and attached."))
  }
  if (hit("cannot open (the connection to|URL) 'https?://|Timeout of [0-9]+ seconds")) {
    return(c("why: an HTTPS download failed (the public epigenomics DA release on",
             "     CloudFront). No credentials are involved; check the network and",
             "     that the URL above resolves."))
  }
  if (hit("gsutil|CommandException|AccessDenied|ServiceException|401|403|credential|not authenticated")) {
    return(c("why: a consortium bucket read failed.",
             "     Needs gsutil on PATH, an authenticated account, and consortium",
             "     access. This is access, not a defect in the panel."))
  }
  if (hit("vendored input missing from the repo")) {
    return(c("why: a file that ships with this repo is not where it should be.",
             "     docs/external_dependencies.md says what it is and where from."))
  }
  if (hit("is missing from")) {
    return("why: an .env file the script reads is not where it looked.")
  }
  if (hit("cannot open|No such file or directory|does not exist")) {
    return("why: a path does not exist. The message above names it.")
  }
  if (hit("subscript out of bounds|undefined columns selected|arguments imply differing number")) {
    return(c("why: an input's shape is not what the script expects, which usually",
             "     means a data-package version other than the pinned one."))
  }
  if (hit("the manifest does not declare|has no column|Reconcile")) {
    return(c("why: the data carries something config/panel_map.json or",
             "     config/table_map.json does not declare. Reconcile the manifest."))
  }
  character(0)
}

#' Run the sub-tables a table script writes, in order, and report each one.
#'
#' The twin of run_panels(). One script per supplementary table number, so
#' tables/ST3.R writes ST3a through ST3h off one data load, and a selection
#' narrows the run:
#'
#'   Rscript figures/landscape/tables/ST3.R            # every sub-table
#'   Rscript figures/landscape/tables/ST3.R ST3g       # one of them
#'
#' Only for the standalone table scripts. The nineteen sub-tables a figure
#' script writes are written by export_table() inside the panel function that
#' computed them, and stay there.
#'
#' One entry may write more than one TSV, so an entry is a unit of work rather
#' than a promise of one file. Prefer one entry PER SUB-TABLE anyway, with the
#' shared work in a memoised accessor above them: every id the script writes is
#' then selectable and is named in the unknown-id error. ST1a and ST1b are two
#' entries over one memoised participant count for that reason, and nothing is
#' computed twice. Merge only where two ids are genuinely one write.
#'
#' @param tables Named list of zero-argument functions; names are sub-table ids,
#'   in write order.
#' @param selection Sub-table ids to run. Empty means all of them.
#' @return Invisibly, a data frame of table, status, output and message.
run_tables <- function(tables, selection = commandArgs(trailingOnly = TRUE)) {
  if (!is.list(tables) || length(tables) == 0) {
    stop("run_tables() needs a non-empty list of sub-table functions",
         call. = FALSE)
  }
  ids <- names(tables)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids) > 0) {
    stop("run_tables() needs every element named, once, with its sub-table id",
         call. = FALSE)
  }
  not_function <- ids[!vapply(tables, is.function, logical(1))]
  if (length(not_function) > 0) {
    stop("run_tables(): not a function: ", paste(not_function, collapse = ", "),
         call. = FALSE)
  }

  selection <- selection[nzchar(selection)]
  unknown <- setdiff(selection, ids)
  if (length(unknown) > 0) {
    stop("unknown sub-table ", paste(unknown, collapse = ", "),
         ". This script writes: ", paste(ids, collapse = ", "), call. = FALSE)
  }
  run <- if (length(selection) > 0) selection else ids

  rows <- lapply(run, function(table) {
    outcome <- tryCatch({
      value <- tables[[table]]()
      # A sub-table whose fit has not been run returns skip_if_no_fit()'s
      # sentinel. It wrote nothing, so reporting OK against a declared output
      # path would claim a file that is not there.
      if (inherits(value, "panel_skip")) {
        list(status = "SKIPPED", output = "", message = value$reason)
      } else {
        list(status = "OK", output = table_result_path(table, value), message = "")
      }
    }, error = function(e) {
      list(status = "FAIL", output = "", message = conditionMessage(e),
           why = explain_failure(e))
    })
    tag <- c(OK = " OK ", FAIL = "FAIL", SKIPPED = "SKIP")[[outcome$status]]
    detail <- if (nzchar(outcome$output)) outcome$output else outcome$message
    # One line, so the id and the raw message stay greppable and parseable; the
    # explanation follows indented and is not part of the outcome line.
    message(sprintf("[%s] %-8s %s", tag, table, detail))
    for (line in outcome$why) message("        ", line)
    data.frame(table = table, status = outcome$status, output = outcome$output,
               message = outcome$message, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  message(sprintf("        %d sub-table(s): %d ok, %d failed, %d skipped",
                  nrow(out), sum(out$status == "OK"),
                  sum(out$status == "FAIL"), sum(out$status == "SKIPPED")))
  if (any(out$status == "FAIL") && !interactive()) {
    quit(status = 1)
  }
  invisible(out)
}

#' Append this export to outputs/tables/export_manifest.tsv.
#'
#' Named apart from panel_export.R's record_export(), which takes a different
#' signature entirely — a panel script sourcing both would otherwise call
#' whichever was sourced last and fail on the argument list. See the header note.
#'
#' The row records the two data-package versions for the same reason the panel
#' manifest does: identical code over v1.3 objects and over v2.0 objects
#' produces two different tables, and nothing in the TSV would say which.
record_table_export <- function(spec, out, df) {
  manifest <- file.path(landscape_root(), table_map()$output$sidecar_manifest)
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)

  row <- data.frame(
    table = spec$table,
    supp_id = spec$supp_id,
    output = spec$output,
    bytes = file.size(out),
    n_row = nrow(df),
    n_col = ncol(df),
    data_pkg_version = installed_version("MotrpacHumanPreSuspensionData"),
    analysis_pkg_version = installed_version("MotrpacHumanPreSuspensionAnalysis"),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    exported_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  # Keep one row per sub-table: a re-run replaces its own row rather than stacking.
  if (file.exists(manifest)) {
    prior <- read.csv(manifest, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE)
    prior <- prior[prior$table != spec$table, , drop = FALSE]
    row <- rbind(prior, row)
    row <- row[order(match(row$table, table_ids())), , drop = FALSE]
  }
  write.table(row, file = manifest, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  invisible(manifest)
}
