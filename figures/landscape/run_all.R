#!/usr/bin/env Rscript
# run_all.R — build every figure and every supplementary table this repo owns.
#
#   Rscript figures/landscape/run_all.R              every figure and table
#   Rscript figures/landscape/run_all.R --fits       run the four fits first
#   Rscript figures/landscape/run_all.R FIG2 ED6     only those figures
#
# The script list comes from config/panel_map.json and config/table_map.json,
# so a figure added to a manifest is picked up here without editing this file.
#
# Each script runs in its own Rscript process. A figure that fails, or that grows
# to several GB loading its differential analysis, cannot affect the next one,
# and the memory is released between them. The cost is reloading packages per
# script, which is seconds against minutes of work.
#
# THE FOUR FITS ARE NOT RUN unless --fits is given. They are hours of compute and
# their output is cached, so the common case is to run them once by hand and then
# rebuild figures freely. Without them the sixteen fit-dependent panels report
# SKIPPED and name the fit to run, which is not a failure.
#
# extended_data_3/celltype_da_fit.R cannot complete without consortium access to
# the CIBERSORTx input, so --fits always reports it as failed in a public
# checkout. That is expected, not a defect in this runner.
#
# single_feature_plots.R is NOT run: all sixteen of its panels are drawn by the
# figure scripts that own them, so running it would rebuild the same PDFs and
# re-stamp the export manifest. Run it directly to rebuild those alone.

suppressPackageStartupMessages(library(jsonlite))

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "lib", "panel_export.R"))

root <- landscape_root()
args <- commandArgs(trailingOnly = TRUE)
with_fits <- "--fits" %in% args
selection <- setdiff(args, "--fits")

# ---- what to run, in order ---------------------------------------------------

# Fits first when asked: a figure that reads one is skipped until it exists.
FITS <- c("extended_data_3/sex_da_fit.R",
          "extended_data_3/celltype_da_fit.R",
          "figure_3/clinical_omics_fit.R",
          "figure_4/plier_fit.R")

pm <- jsonlite::read_json(file.path(root, "config", "panel_map.json"),
                          simplifyVector = FALSE)
tm <- jsonlite::read_json(file.path(root, "config", "table_map.json"),
                          simplifyVector = FALSE)

# Figure scripts in manifest order, each once. A figure's panels all name the
# same script, so unique() over the panel list is the figure list.
figure_scripts <- unique(vapply(pm$panels, function(p) p$script, character(1)))

# Only the standalone table scripts. The nineteen sub-tables a figure script
# writes come out of that figure's run and must not be run again here.
table_scripts <- unique(Filter(Negate(is.null),
  lapply(tm$tables, function(t) {
    if (is.null(t$script) || !isTRUE(t$built_here)) return(NULL)
    if (t$script %in% figure_scripts) return(NULL)
    t$script
  })))
table_scripts <- unlist(table_scripts)

scripts <- c(if (with_fits) FITS, figure_scripts, table_scripts)

# A bare argument selects by figure id, folder or path substring.
if (length(selection) > 0) {
  keep <- vapply(scripts, function(s) {
    any(vapply(selection, function(sel)
      grepl(sel, s, fixed = TRUE) ||
      identical(sub("\\.R$", "", basename(s)), sel), logical(1)))
  }, logical(1))
  unmatched <- selection[!vapply(selection, function(sel)
    any(grepl(sel, scripts, fixed = TRUE) |
        sub("\\.R$", "", basename(scripts)) == sel), logical(1))]
  if (length(unmatched) > 0) {
    stop("nothing matches: ", paste(unmatched, collapse = ", "),
         "\n  Known: ", paste(sub("\\.R$", "", basename(scripts)), collapse = ", "),
         call. = FALSE)
  }
  scripts <- scripts[keep]
}

missing <- scripts[!file.exists(file.path(root, scripts))]
if (length(missing) > 0) {
  stop("the manifest names a script that is not on disk: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

# ---- run ---------------------------------------------------------------------

# run_panels() and run_tables() both print "[ OK ] ID detail". Parsing their
# output rather than re-deriving it keeps one definition of what happened.
OUTCOME <- "^\\[( OK |FAIL|SKIP)\\] +(\\S+) *(.*)$"

rows <- list()
started <- Sys.time()

for (s in scripts) {
  message("\n==> ", s)
  out <- suppressWarnings(
    system2("Rscript", shQuote(file.path(root, s)),
            stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  hits <- grep(OUTCOME, out, value = TRUE)
  for (h in hits) {
    m <- regmatches(h, regexec(OUTCOME, h))[[1]]
    rows[[length(rows) + 1L]] <- data.frame(
      script = s,
      id = m[3],
      status = c(" OK " = "OK", "FAIL" = "FAIL", "SKIP" = "SKIPPED")[[m[2]]],
      detail = trimws(m[4]),
      stringsAsFactors = FALSE
    )
  }
  message(paste(hits, collapse = "\n"))

  # A script that dies before it reports anything leaves no rows of its own.
  # This is the case worth explaining: it means the failure was in the header,
  # not in a panel, so nothing built and the usual per-panel line never printed.
  if (status != 0L && length(hits) == 0L) {
    # In order of usefulness: the fit scripts report their own checks as
    # "[<fit>] FAIL <what>" and never raise, so an R-level Error is the SECOND
    # place to look, not the first. Falling straight to tail() picks up a
    # package startup warning and whatever passed last, which says nothing.
    err <- grep("\\bFAIL\\b", out, value = TRUE)
    if (length(err) == 0L) err <- grep("^(Error|Fatal|Calls:)", out, value = TRUE)
    if (length(err) == 0L) {
      err <- grep("^\\s*(Warning|package .* was built under)", out,
                  value = TRUE, invert = TRUE)
      err <- utils::tail(err[nzchar(trimws(err))], 3)
    }
    reason <- paste(utils::tail(err, 2), collapse = " | ")
    why <- explain_failure(simpleError(paste(out, collapse = "\n")))
    rows[[length(rows) + 1L]] <- data.frame(
      script = s, id = NA_character_, status = "FAIL",
      detail = paste0("script exited ", status, " before building anything: ", reason),
      stringsAsFactors = FALSE)
    message("[FAIL] ", s, " exited ", status, " before building anything")
    message("        ", reason)
    for (line in why) message("        ", line)
    message("        the failure is in the script's setup, not in a panel.")
  }
}

report <- if (length(rows) > 0) do.call(rbind, rows) else
  data.frame(script = character(), id = character(),
             status = character(), detail = character())

# ---- report -------------------------------------------------------------------

message("\n", strrep("-", 72))
for (s in scripts) {
  r <- report[report$script == s, , drop = FALSE]
  message(sprintf("  %-34s %d ok, %d failed, %d skipped", s,
                  sum(r$status == "OK"), sum(r$status == "FAIL"),
                  sum(r$status == "SKIPPED")))
}

failed <- report[report$status == "FAIL", , drop = FALSE]
if (nrow(failed) > 0) {
  message("\nFailed:")
  for (i in seq_len(nrow(failed))) {
    message(sprintf("  %-10s %-34s %s",
                    failed$id[i], failed$script[i],
                    substr(failed$detail[i], 1, 90)))
  }
}

message(sprintf("\n%d script(s) in %.1f min: %d ok, %d failed, %d skipped",
                length(scripts),
                as.numeric(difftime(Sys.time(), started, units = "mins")),
                sum(report$status == "OK"), nrow(failed),
                sum(report$status == "SKIPPED")))

out_dir <- file.path(root, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
report_path <- file.path(out_dir, "run_report.tsv")
utils::write.table(report, report_path, sep = "\t", row.names = FALSE,
                   quote = FALSE, na = "")
message("report: ", sub(paste0(root, "/"), "", report_path, fixed = TRUE))

if (nrow(failed) > 0 && !interactive()) quit(status = 1)
