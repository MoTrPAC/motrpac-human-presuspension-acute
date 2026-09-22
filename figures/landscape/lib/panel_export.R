# panel_export.R — the export contract every panel script in this repo uses.
#
# A panel script never names its own output file. It calls panel_init() to learn
# which panel it is building, and export_panel() to write it. Both read
# config/panel_map.json, so renumbering a panel (ED5A -> ED4B, say) is an edit to
# that file alone and the exported PDF name follows.
#
# What the device settings are for: these PDFs are opened in Adobe Illustrator.
#   onefile = TRUE       a second plot becomes a second PAGE rather than
#                        overwriting the file. This looks backwards and is not:
#                        with onefile = FALSE and a filename carrying no %d,
#                        base pdf() rewrites the same path on every new page, so
#                        a panel that draws twice silently ships only its last
#                        plot. Keeping the extra page turns silent loss into a
#                        defect that can be seen. Illustrator opens page 1 of a
#                        multi-page PDF and drops the rest, so more than one
#                        page is still a defect - it is just a detectable one
#                        now.
#   useDingbats = FALSE  point glyphs drawn as circles rather than as ZapfDingbats
#                        characters, which Illustrator renders as random letters.
#   family = "sans"      maps to Helvetica, a base-14 PDF font. It is referenced,
#                        not embedded, which is what keeps the text editable.
#                        Panels needing glyphs outside Latin-1 (Greek, arrows,
#                        typographic minus) must set device = "cairo_pdf" in
#                        panel_map.json; base pdf() cannot encode them.
#
# One script draws every panel of one figure. Each panel is a zero-argument
# function and run_panels() runs them, so the data load at the top of the script
# happens once for the whole figure and a single panel can still be rebuilt on
# its own.
#
# Usage in a figure script:
#
#   here <- dirname(sub("^--file=", "",
#                       grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
#   source(file.path(here, "..", "lib", "panel_export.R"))
#
#   ...load what the figure shares...
#
#   fig2a <- function() {
#     spec <- panel_init("FIG2A")            # sets the seed, returns the spec
#     ...build p...
#     export_panel(p, "FIG2A")               # ggplot / patchwork
#   }
#   fig2b <- function() {
#     panel_init("FIG2B")
#     export_panel(function() draw(ht), "FIG2B")  # anything drawn by side effect
#   }
#
#   run_panels(list(FIG2A = fig2a, FIG2B = fig2b))
#
#   Rscript figures/landscape/FIG2.R            # every panel
#   Rscript figures/landscape/FIG2.R FIG2B      # one panel

# ---- the pipe ---------------------------------------------------------------

# %>% is bound here, once, because every panel script and every figure helper
# sources this file before anything else. An operator has to be in scope to be
# parsed, so a helper that is otherwise fully namespace-qualified would
# otherwise have to attach dplyr just to use one. Binding the object adds
# nothing to the caller's search path.
#
# Panels that attach dplyr themselves get the identical object; a local binding
# and an attached one cannot disagree.
`%>%` = dplyr::`%>%`

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- the landscape root -----------------------------------------------------

#' Absolute path to figures/landscape/.
#'
#' Every path config/ carries is relative to it: the panel outputs, the vendored
#' sources/, the fits under outputs/. Resolved in this order:
#'
#'   LANDSCAPE_ROOT, if it is set;
#'   the directory of the running script, walking up — a figure script sits in
#'     the root, a tables/, analysis/ or helpers/ script one level below it;
#'   the working directory, if it is the root.
#'
#' A directory is the root when it holds lib/panel_export.R and
#' config/panel_map.json.
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
#'
#' NULL rather than an error: a session that has no --file (R -e, a console)
#' falls through to the working directory.
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

# ---- config/landscape.env ---------------------------------------------------

#' Read config/landscape.env into the environment, once per session.
#'
#' Shell syntax, because it is a list of settings and reads like one; parsed
#' here rather than run through a shell so an R session needs no subprocess.
#' Values already exported WIN, which is how a one-off override works:
#'
#'   RASTER_DPI=150 Rscript figures/landscape/FIG2.R FIG2A
#'
#' A missing file is not an error: everything the file declares has a default in
#' the function that reads it.
load_landscape_env <- function(root = NULL) {
  root <- root %||% tryCatch(landscape_root(), error = function(e) NULL)
  if (is.null(root)) {
    return(invisible(FALSE))
  }
  path <- file.path(root, "config", "landscape.env")
  if (!file.exists(path)) {
    return(invisible(FALSE))
  }
  for (line in readLines(path, warn = FALSE)) {
    line <- trimws(line)
    if (!startsWith(line, "export ")) next
    body <- sub("^export ", "", line)
    name <- sub("=.*$", "", body)
    if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", name)) next
    # LANDSCAPE_ROOT is the root this file was found under, so it is already
    # known and its own line is a shell self-locating expression this parser
    # cannot evaluate. Reading it back would overwrite a resolved path with an
    # unevaluated one.
    if (identical(name, "LANDSCAPE_ROOT")) next
    if (nzchar(Sys.getenv(name, unset = ""))) next          # already set: leave it
    value <- sub("^[^=]*=", "", body)
    value <- sub('^"', "", sub('"$', "", value))
    # "${NAME:-default}" -> default
    value <- sub("^\\$\\{[A-Za-z0-9_]+:-", "", value)
    value <- sub("\\}$", "", value)
    value <- gsub("\\$\\{?LANDSCAPE_ROOT\\}?", root, value)
    # This is not a shell. A value still carrying a substitution, a backtick or
    # a nested expansion cannot be resolved here, and setting the literal is
    # worse than leaving the variable unset for the code default to handle.
    if (grepl("[$`]", value)) next
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(TRUE)
}

load_landscape_env()

.panel_map_cache <- new.env(parent = emptyenv())

# ---- the manifest ----------------------------------------------------------

#' Read config/panel_map.json, once per session.
panel_map <- function() {
  if (!is.null(.panel_map_cache$map)) {
    return(.panel_map_cache$map)
  }
  path <- file.path(landscape_root(), "config", "panel_map.json")
  if (!file.exists(path)) {
    stop("panel manifest not found: ", path, call. = FALSE)
  }
  map <- jsonlite::read_json(path, simplifyVector = FALSE)
  .panel_map_cache$map <- map
  map
}

#' Every panel id in the manifest, in manifest order.
panel_ids <- function() {
  vapply(panel_map()$panels, function(p) p$panel, character(1))
}

#' The directory the ATAC and methylCap QC matrices are cached in.
#'
#' Neither data package ships them; load_qc(epigen = TRUE) downloads them on
#' first use and reuses the cache afterwards. Panels pass this as its
#' repo_local_dir. EPIGEN_QC_DIR overrides the default. The cache sits under
#' staging/ rather than outputs/: it is 7.7 GB of downloaded consortium data,
#' not something a run produces, and a clean of outputs/ must not discard it.
epigen_qc_dir <- function() {
  dir <- Sys.getenv("EPIGEN_QC_DIR", unset = "")
  if (!nzchar(dir)) {
    dir <- file.path(landscape_root(), "staging", "raw-files")
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dir)) {
    stop("EPIGEN_QC_DIR names a directory that cannot be created: ", dir,
         call. = FALSE)
  }
  dir
}

#' The rasterisation resolution set in config/landscape.env, or NULL.
#'
#' NULL rather than a default, so panel_spec() can fall through to the manifest
#' when the variable is not set.
env_raster_dpi <- function() {
  v <- Sys.getenv("RASTER_DPI", unset = "")
  if (!nzchar(v)) return(NULL)
  dpi <- suppressWarnings(as.integer(v))
  if (is.na(dpi) || dpi < 1) {
    stop("RASTER_DPI must be a positive integer, got '", v, "'", call. = FALSE)
  }
  dpi
}

#' The manifest entry for one panel, with figure-level fields resolved onto it.
#'
#' @param panel Panel id, e.g. "FIG2A". Must appear in config/panel_map.json.
panel_spec <- function(panel) {
  map <- panel_map()
  ids <- panel_ids()
  if (!panel %in% ids) {
    stop(
      "unknown panel '", panel, "'. Known panels: ", paste(ids, collapse = ", "),
      call. = FALSE
    )
  }
  spec <- map$panels[[match(panel, ids)]]

  figure_ids <- vapply(map$figures, function(f) f$figure_id, character(1))
  figure <- map$figures[[match(spec$figure_id, figure_ids)]]
  spec$figure_dir <- figure$figure_dir
  spec$figure_title <- figure$title

  spec$seed <- spec$seed %||% map$defaults$seed
  spec$device <- spec$device %||% map$defaults$device
  spec$font_family <- spec$font_family %||% map$defaults$font_family
  # RASTER_DPI wins over the manifest. Rasterisation resolution is a render
  # policy that applies to every panel that rasterises anything, not a fact
  # about one panel, so it is set once in config/landscape.env. The manifest
  # value is the fallback.
  spec$raster_dpi <- env_raster_dpi() %||% spec$raster_dpi %||% map$defaults$raster_dpi
  spec$width_in <- spec$width_in %||% 7
  spec$height_in <- spec$height_in %||% 5
  spec
}

#' Absolute path of the PDF a panel is required to write.
panel_output_path <- function(panel) {
  spec <- panel_spec(panel)
  file.path(landscape_root(), spec$output)
}

# ---- vendored inputs -------------------------------------------------------

#' Path to an input vendored under a figure's sources/ directory.
#'
#' The file ships with the repo, so a panel needs nothing supplied to build. An
#' environment variable may still override it, which is how you point a panel at
#' a newer copy without editing anything.
#'
#' @param figure_dir Figure directory, e.g. "figure_3".
#' @param file_name File name inside sources/<figure_dir>/.
#' @param env_var Optional variable that overrides the vendored copy.
#' The folder a figure's scripts and sources live in.
#'
#' A main figure and its Extended Data supplement share one folder, so this is
#' not figure_dir. The manifest declares it per figure.
figure_folder <- function(figure_dir) {
  for (f in panel_map()$figures) {
    if (identical(f$figure_dir, figure_dir)) return(f$folder)
  }
  stop("no figure with figure_dir '", figure_dir, "' in config/panel_map.json",
       call. = FALSE)
}

panel_source <- function(figure_dir, file_name, env_var = NULL) {
  if (!is.null(env_var)) {
    override <- Sys.getenv(env_var, unset = "")
    if (nzchar(override)) {
      if (!file.exists(override)) {
        stop(env_var, " names a file that is not there: ", override, call. = FALSE)
      }
      return(override)
    }
  }
  # Each figure pair's folder carries its own sources/, so the vendored input
  # sits beside the script that reads it. figure_dir is per-FIGURE and the
  # folder is per-PAIR, so the manifest's `folder` is what resolves it.
  folder <- figure_folder(figure_dir)
  path <- file.path(landscape_root(), folder, "sources", file_name)
  if (!file.exists(path)) {
    stop("vendored input missing from the repo: ",
         file.path(folder, "sources", file_name),
         if (!is.null(env_var)) paste0(" (or set ", env_var, ")") else "",
         call. = FALSE)
  }
  path
}

# ---- starting a panel ------------------------------------------------------

#' Attach the data packages a panel declares in the manifest.
#'
#' These two have to be ATTACHED, not merely namespace-qualified. Both resolve
#' their lazy-loaded data objects with `eval(parse(text = "MUSCLE_TRNSCRPT_DA"))`
#' and similar, which reaches the package's lazydata environment only through the
#' search path. `MotrpacHumanPreSuspensionAnalysis::load_differential_analysis()`
#' called on a package that is merely loaded fails with
#' *object 'ADIPOSE_METAB_DA' not found*, hundreds of lines into a build.
#'
#' Driving this from the manifest's `data_packages` rather than from a
#' `library()` line in each script means one declaration governs both the attach
#' and the record: the same field is what the export manifest reports the
#' version of.
attach_data_packages <- function(spec) {
  for (pkg in unlist(spec$data_packages)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("panel '", spec$panel, "' needs ", pkg, ", which is not installed",
           call. = FALSE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, warn.conflicts = FALSE)
    )
  }
  invisible(NULL)
}

#' Announce a panel and make its run reproducible.
#'
#' Sets the RNG seed from the manifest before any code that could consume random
#' numbers. Panels with an empty "stochastic" list are seeded anyway: a seed that
#' is set and unused costs nothing, and a dependency that starts sampling
#' internally should not silently make a panel irreproducible.
#'
#' Also attaches the panel's declared data packages — see attach_data_packages().
#'
#' @param panel Panel id, e.g. "FIG2A".
panel_init <- function(panel) {
  spec <- panel_spec(panel)
  set.seed(spec$seed)
  attach_data_packages(spec)
  message(sprintf(
    "[%s] %s\n        seed=%d  ->  %s",
    spec$panel, spec$title, spec$seed, spec$output
  ))
  if (length(spec$stochastic) > 0) {
    message(sprintf(
      "        stochastic steps under this seed: %s",
      paste(unlist(spec$stochastic), collapse = ", ")
    ))
  }
  invisible(spec)
}

# ---- exporting -------------------------------------------------------------

#' Write a panel to its manifest-declared PDF.
#'
#' @param x What to draw. A ggplot or patchwork object; a grid grob; a
#'   ComplexHeatmap object; or a zero-argument function whose side effect is the
#'   drawing (the escape hatch for base graphics and for anything that needs
#'   draw() called with arguments).
#' @param panel Panel id, e.g. "FIG2A".
#' @param width,height Inches. Default to the manifest values; pass explicitly
#'   only to override a size the manifest gets wrong, and update the manifest.
#' @param device One of "pdf" or "cairo_pdf". Defaults to the manifest value.
#'   "pdf" keeps text editable in Illustrator by referencing rather than
#'   embedding Helvetica; "cairo_pdf" embeds a font subset and is required only
#'   for glyphs base pdf() cannot encode.
export_panel <- function(x,
                         panel,
                         width = NULL,
                         height = NULL,
                         device = c("pdf", "cairo_pdf")) {
  spec <- panel_spec(panel)
  device <- if (missing(device)) spec$device else match.arg(device)
  device <- match.arg(device, c("pdf", "cairo_pdf"))
  width <- width %||% spec$width_in
  height <- height %||% spec$height_in

  out <- panel_output_path(panel)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

  # A previous good export is not this call's to destroy. If the drawing fails,
  # the file is restored rather than removed, so a transient failure does not
  # cost the artifact that was already on disk.
  prior <- if (file.exists(out)) readBin(out, "raw", file.size(out)) else NULL

  if (identical(device, "cairo_pdf")) {
    if (!capabilities("cairo")) {
      stop(
        "panel '", panel, "' asks for cairo_pdf but this R has no cairo support",
        call. = FALSE
      )
    }
    grDevices::cairo_pdf(
      filename = out, width = width, height = height,
      family = spec$font_family, onefile = TRUE
    )
  } else {
    grDevices::pdf(
      file = out, width = width, height = height,
      family = spec$font_family, onefile = TRUE, useDingbats = FALSE
    )
  }
  # Close the device this call opened, not whatever is current: drawing code
  # that opens a device of its own would otherwise leave ours open and this PDF
  # truncated, while the call reported success.
  own_dev <- grDevices::dev.cur()
  close_own <- function() {
    if (own_dev %in% grDevices::dev.list()) {
      try(grDevices::dev.off(own_dev), silent = TRUE)
    }
  }

  drawn <- FALSE
  on.exit({
    close_own()
    if (!drawn) {
      if (is.null(prior)) unlink(out) else writeBin(prior, out)
    }
  }, add = TRUE)

  # Base pdf() cannot encode anything outside Latin-1. It does not error: it
  # warns 'conversion failure ... in mbcsToSbcs' and silently drops the glyph,
  # or substitutes an ASCII lookalike. Promoting that warning to an error is the
  # difference between catching it here and shipping a legend that reads ">= 0.05"
  # or a count with its direction arrow missing.
  withCallingHandlers(
    draw_panel(x),
    warning = function(w) {
      if (grepl("mbcsToSbcs|conversion failure", conditionMessage(w))) {
        stop(
          "panel '", panel, "' draws a character base pdf() cannot encode: ",
          conditionMessage(w),
          "\n  Set \"device\": \"cairo_pdf\" for this panel in config/panel_map.json.",
          call. = FALSE
        )
      }
      invokeRestart("muffleWarning")
    }
  )
  drawn <- TRUE

  on.exit(NULL)
  close_own()

  if (!file.exists(out) || file.size(out) == 0) {
    stop("panel '", panel, "' produced no PDF at ", out, call. = FALSE)
  }
  # A device closed by the drawing code, or a second page, leaves a file that
  # looks plausible. Check the bytes before claiming the export.
  bytes <- readBin(out, "raw", file.size(out))
  if (length(grepRaw("%%EOF", bytes, all = TRUE, fixed = TRUE)) == 0) {
    stop("panel '", panel, "' wrote a truncated PDF at ", out,
         " — the device was closed by something other than export_panel()",
         call. = FALSE)
  }

  record_export(spec, out, width, height, device)
  message(sprintf(
    "        wrote %s (%.1f KB, %g x %g in, device=%s)",
    spec$output, file.size(out) / 1024, width, height, device
  ))
  invisible(out)
}

#' Installed version of a package, or NA when it is not installed.
installed_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)),
           error = function(e) NA_character_)
}

#' Draw whatever a panel handed to export_panel().
draw_panel <- function(x) {
  if (is.function(x)) {
    x()
  } else if (inherits(x, c("Heatmap", "HeatmapList", "HeatmapAnnotation"))) {
    ComplexHeatmap::draw(x)
  } else if (inherits(x, "gg") || inherits(x, "patchwork")) {
    print(x)
  } else if (inherits(x, c("grob", "gTree", "gList"))) {
    grid::grid.newpage()
    grid::grid.draw(x)
  } else if (inherits(x, "recordedplot")) {
    grDevices::replayPlot(x)
  } else {
    print(x)
  }
  invisible(NULL)
}

#' Append one row to outputs/panels/export_manifest.tsv.
#'
#' A row without a file means an export was claimed and lost, and a file without
#' a row means a panel wrote a PDF outside the export contract.
#'
#' The row carries the versions of the two data packages the panel was drawn
#' from. That is the fact a PDF cannot otherwise state and the one that decides
#' what the figure means: identical code over two versions of the objects
#' produces two different figures and nothing in the file would say so.
record_export <- function(spec, out, width, height, device) {
  manifest <- file.path(landscape_root(), panel_map()$output$sidecar_manifest)
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)

  row <- data.frame(
    panel = spec$panel,
    figure_id = spec$figure_id,
    output = spec$output,
    bytes = file.size(out),
    width_in = width,
    height_in = height,
    device = device,
    seed = spec$seed,
    data_pkg_version = installed_version("MotrpacHumanPreSuspensionData"),
    analysis_pkg_version = installed_version("MotrpacHumanPreSuspensionAnalysis"),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    exported_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    stringsAsFactors = FALSE
  )

  # Keep one row per panel: a re-run replaces its own row rather than stacking.
  if (file.exists(manifest)) {
    prior <- read.csv(manifest, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE)
    prior <- prior[prior$panel != spec$panel, , drop = FALSE]
    row <- rbind(prior, row)
    row <- row[order(match(row$panel, panel_ids())), , drop = FALSE]
  }
  write.table(row, file = manifest, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  invisible(manifest)
}

# ---- running a figure's panels ----------------------------------------------

#' A panel that cannot be run here, as a value a panel function returns.
panel_skip <- function(reason) {
  structure(list(reason = reason), class = "panel_skip")
}

is_panel_skip <- function(x) inherits(x, "panel_skip")

#' Skip a panel whose input is not in this repository.
#'
#' Three panels read individual-level inputs that are not released with this
#' repository. Missing files are ordinarily an error; for these three they are
#' the expected state for anyone without a data-access approval, so the panel
#' says so and the rest of the figure still builds.
#'
#' The panel calls this first and hands back what it returns:
#'
#'   ed1a <- function() {
#'     skip <- skip_if_gated(
#'       "ED1A",
#'       path = "figure_1/sources/1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS",
#'       env_var = "ED1A_PCA_RDS",
#'       input = "a genotype PCA carrying 238 MoTrPAC participant IDs and 32 genotype principal components per individual",
#'       access = "an approved MoTrPAC consortium data-access request"
#'     )
#'     if (!is.null(skip)) return(skip)
#'     ...
#'   }
#'
#' @param panel Panel id, for the message.
#' @param path Where the input would be: absolute, or relative to the landscape
#'   root.
#' @param input What the file is and what makes it individual-level.
#' @param access What it takes to obtain a copy.
#' @param env_var Variable naming a local copy; checked before `path`.
#' @return NULL when the input is there, otherwise a skip sentinel.
skip_if_gated <- function(panel, path, input, access, env_var = NULL) {
  candidates <- character(0)
  if (!is.null(env_var)) {
    from_env <- Sys.getenv(env_var, unset = "")
    if (nzchar(from_env)) {
      candidates <- c(candidates, from_env)
    }
  }
  full <- if (grepl("^(/|~|[A-Za-z]:)", path)) path else file.path(landscape_root(), path)
  candidates <- c(candidates, full)
  if (any(file.exists(candidates))) {
    return(invisible(NULL))
  }

  message(paste(
    sprintf("%s CANNOT BE RUN BY THE PUBLIC.", panel),
    sprintf("        It reads %s.", input),
    "        That is individual-level participant data. It is not in this",
    "        repository, it is not distributed with it, and it may not be",
    "        published. No copy of it will be found by a public checkout.",
    sprintf("        Running %s requires %s,", panel, access),
    sprintf("        and a local copy of the file at %s", full),
    if (!is.null(env_var)) sprintf("        (or %s pointing at one).", env_var) else "",
    "        Every other panel in this figure runs without it.",
    sep = "\n"
  ))
  panel_skip(sprintf("input not in this repository: %s", basename(path)))
}

#' Skip a panel whose fit has not been run yet.
#'
#' Sixteen panels read one of the four fits under analysis/. A pipeline stage
#' once produced those before any panel ran, which is why the manifest marks
#' them `must_exist`. Here the reader runs the fit themselves, and an unrun fit is an
#' ordinary state rather than a broken checkout, so the panel says which script
#' to run and the rest of the figure still builds.
#'
#' Manifest-driven: the panel repeats nothing. Every `must_exist` entry in the
#' panel's `external_inputs` is checked, and `produced_by` names the script.
#'
#'   fig4g <- function() {
#'     skip <- skip_if_no_fit("FIG4G")
#'     if (!is.null(skip)) return(skip)
#'     ...
#'   }
#'
#' @param id Panel or sub-table id.
#' @return NULL when every required fit is there, otherwise a skip sentinel.
skip_if_no_fit <- function(id) {
  # Panels and sub-tables both declare must_exist inputs, and both manifests
  # carry the same fields, so one function serves either. table_spec() lives in
  # table_export.R, which a figure script need not source.
  spec <- tryCatch(panel_spec(id), error = function(e) NULL)
  if (is.null(spec) && exists("table_spec", mode = "function")) {
    spec <- tryCatch(table_spec(id), error = function(e) NULL)
  }
  if (is.null(spec)) {
    stop(id, " is in neither config/panel_map.json nor config/table_map.json.",
         call. = FALSE)
  }
  inputs <- spec$external_inputs
  if (length(inputs) == 0) return(invisible(NULL))

  missing <- list()
  for (input in inputs) {
    if (!isTRUE(input$must_exist)) next
    var <- input$env_var
    if (is.null(var) || !nzchar(var)) next
    # override_env_var names a copy produced elsewhere. It satisfies the input
    # on its own: ED3E's is the documented way to draw it without the gated
    # deconvolution, and checking env_var alone would skip a panel that can
    # draw.
    override <- input$override_env_var
    if (!is.null(override) && nzchar(override)) {
      from_override <- Sys.getenv(override, unset = "")
      if (nzchar(from_override) && file.exists(from_override)) next
    }
    where <- Sys.getenv(var, unset = "")
    if (nzchar(where) && file.exists(where)) next
    missing[[length(missing) + 1L]] <- list(var = var, where = where,
                                            override = override,
                                            produced_by = input$produced_by)
  }
  if (length(missing) == 0) return(invisible(NULL))

  # produced_by opens with the script that writes it; the rest is prose. Matched
  # as <folder>/<name>.R rather than against a fixed directory, so moving a fit
  # between folders does not silently turn this into the whole paragraph.
  script <- function(text) {
    hit <- regmatches(text, regexpr("^[A-Za-z0-9_]+/[A-Za-z0-9_.]+\\.R", text))
    if (length(hit) == 1L) hit else text
  }
  for (m in missing) {
    message(paste(
      sprintf("%s needs a fit that has not been run.", id),
      sprintf("        %s is %s.", m$var,
              if (nzchar(m$where)) sprintf("set to %s, which is not there", m$where) else "unset"),
      sprintf("        Run:  Rscript figures/landscape/%s", script(m$produced_by)),
      if (!is.null(m$override) && nzchar(m$override))
        sprintf("        Or point %s at a copy produced elsewhere.", m$override) else "",
      "        The fits are expensive and are cached, so this is a one-time cost.",
      sep = "\n"
    ))
  }
  panel_skip(sprintf("fit not run: %s",
                     paste(vapply(missing, function(m) m$var, character(1)),
                           collapse = ", ")))
}

#' Where a panel wrote, for the run report.
#'
#' export_panel() returns the path, so a panel function that ends with it says
#' where it wrote; anything else falls back to what the manifest declares.
panel_result_path <- function(panel, value) {
  if (is.character(value) && length(value) == 1L && nzchar(value)) {
    return(value)
  }
  tryCatch(panel_output_path(panel), error = function(e) "")
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

#' Run the panels a figure script defines, in order, and report each one.
#'
#' The shared data load happens once, above this call, rather than once per
#' panel; that is what one script per figure buys. A selection narrows the run:
#'
#'   Rscript figures/landscape/FIG2.R            # every panel
#'   Rscript figures/landscape/FIG2.R FIG2B      # one panel
#'
#' Each panel runs inside tryCatch, so a panel that fails costs its own PDF and
#' not the figure's. The run exits 1 if any panel failed, which is how a
#' non-interactive caller learns about a failure it did not read.
#'
#' @param panels Named list of zero-argument functions; names are panel ids, in
#'   draw order.
#' @param selection Panel ids to run. Empty means all of them.
#' @return Invisibly, a data frame of panel, status, output and message.
run_panels <- function(panels, selection = commandArgs(trailingOnly = TRUE)) {
  if (!is.list(panels) || length(panels) == 0) {
    stop("run_panels() needs a non-empty list of panel functions", call. = FALSE)
  }
  ids <- names(panels)
  if (is.null(ids) || any(!nzchar(ids)) || anyDuplicated(ids) > 0) {
    stop("run_panels() needs every element named, once, with its panel id",
         call. = FALSE)
  }
  not_function <- ids[!vapply(panels, is.function, logical(1))]
  if (length(not_function) > 0) {
    stop("run_panels(): not a function: ", paste(not_function, collapse = ", "),
         call. = FALSE)
  }

  selection <- selection[nzchar(selection)]
  unknown <- setdiff(selection, ids)
  if (length(unknown) > 0) {
    stop("unknown panel ", paste(unknown, collapse = ", "),
         ". This script builds: ", paste(ids, collapse = ", "), call. = FALSE)
  }
  run <- if (length(selection) > 0) selection else ids

  rows <- lapply(run, function(panel) {
    outcome <- tryCatch({
      value <- panels[[panel]]()
      if (is_panel_skip(value)) {
        list(status = "SKIPPED", output = "", message = value$reason)
      } else {
        list(status = "OK", output = panel_result_path(panel, value), message = "")
      }
    }, error = function(e) {
      list(status = "FAIL", output = "", message = conditionMessage(e),
           why = explain_failure(e))
    })
    tag <- c(OK = " OK ", FAIL = "FAIL", SKIPPED = "SKIP")[[outcome$status]]
    detail <- if (nzchar(outcome$output)) outcome$output else outcome$message
    # One line, so the id and the raw message stay greppable and parseable; the
    # explanation follows indented and is not part of the outcome line.
    message(sprintf("[%s] %-8s %s", tag, panel, detail))
    for (line in outcome$why) message("        ", line)
    data.frame(panel = panel, status = outcome$status, output = outcome$output,
               message = outcome$message, stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  message(sprintf("        %d panel(s): %d ok, %d failed, %d skipped",
                  nrow(out), sum(out$status == "OK"), sum(out$status == "FAIL"),
                  sum(out$status == "SKIPPED")))
  if (any(out$status == "FAIL") && !interactive()) {
    quit(status = 1)
  }
  invisible(out)
}

# ---- rasterisation ---------------------------------------------------------

#' Rasterise a layer only when its point count would bloat the PDF.
#'
#' A fully vector scatter of 10^6 points is a 100 MB PDF that Illustrator cannot
#' open; a rasterised one is a few MB with the axes, text and legend still
#' vector. Panels call this instead of reaching for ggrastr directly so the dpi
#' comes from panel_spec() (RASTER_DPI, else the manifest) and one threshold
#' governs every panel.
#'
#' THE THRESHOLD IS SOMEWHAT ARBITRARY. 3,000 marks is a judgement call, not a
#' break-even: it sits well below where rasterising starts to pay on size alone,
#' and it is set there deliberately, to rasterise sooner rather than later.
#'
#' What the size arithmetic says, for context. The comparison that matters is
#' UNCOMPRESSED vector against the raster as it ships, because Illustrator
#' expands and parses every vector mark while a rasterised layer stays a
#' compressed bitmap it decodes once. Measured on a full-page scatter, with
#' pdf(compress = FALSE) for the vector side:
#'
#'   vector          204 bytes per mark, linear
#'   raster 450 dpi  ~2.5 MB, near flat in mark count
#'   raster 600 dpi  ~3.8 MB, near flat in mark count
#'
#' which puts break-even near 13,000 marks at 450 dpi and 19,400 at 600. The
#' default here is far below both, so between 3,000 and ~19,000 marks it trades
#' a larger file for a panel that opens and redraws. That trade is the arbitrary
#' part: file size is measurable and the editing experience is not, so there is
#' no number to derive.
#'
#' Two things the per-layer arithmetic cannot see, both of which argue for
#' rasterising earlier than break-even. ggrastr emits one image PER FACET, and
#' those per-facet images are sparse and compress far better than the dense
#' full-page scatter measured above - ED3A's 36 images total 1.3 MB, not
#' 36 x 3.8 MB. And `n` is one layer's mark count, so a panel drawing three dense
#' layers off one data frame carries three times what any single call sees.
#'
#' @param layer A ggplot layer, e.g. geom_point(...).
#' @param n Number of points the layer will draw.
#' @param panel Panel id, for the dpi.
#' @param threshold Mark count above which rasterising wins. 0 forces it, which
#'   is what a panel reproducing a legacy artifact that rasterised
#'   unconditionally passes.
maybe_rasterise <- function(layer, n, panel, threshold = 3000) {
  if (n <= threshold) {
    return(layer)
  }
  spec <- panel_spec(panel)
  message(sprintf(
    "        rasterising a %d-point layer at %d dpi (a vector layer this dense is a PDF Illustrator cannot open)",
    n, spec$raster_dpi
  ))
  # dev = "ragg" rather than ggrastr's cairo default: output is indistinguishable
  # at these resolutions, it is faster, and it keeps rasterisation from silently
  # depending on cairo support that only the cairo_pdf panels check for.
  ggrastr::rasterise(layer, dpi = spec$raster_dpi, dev = "ragg")
}

# ---- pheatmap colour-bar titles ---------------------------------------------

# pheatmap has no legend title. These two put one above the colour bar, and are
# here rather than beside one figure's panels because FIG4H, FIG5C, FIG5H,
# ED6D, ED6F and ED6I all want it.

#' Where pheatmap put its colour bar: the left edge of the legend column, in
#' inches from the left of the current viewport. Read off the gtable rather than
#' assumed, because a heatmap with column annotation lays out different columns
#' from one without.
#'
#' Must be called inside the viewport the gtable is drawn in: the matrix column
#' is a null unit, so its width is only defined there.
heatmap_legend_left <- function(gtable) {
  cell <- which(gtable$layout$name == "legend")
  if (length(cell) != 1L) {
    stop("expected exactly one legend cell in the pheatmap gtable, found ",
         length(cell), call. = FALSE)
  }
  left_of <- seq_len(gtable$layout$l[cell] - 1L)
  sum(grid::convertWidth(gtable$widths[left_of], "in", valueOnly = TRUE))
}

#' Draw a colour-bar title at `left` inches, pulled back if the word is wider
#' than the column and would run off the page.
#'
#' Call it in the viewport `left` was measured in.
draw_heatmap_legend_title <- function(left, title, fontsize = 10, top_mm = 1.5) {
  label <- grid::textGrob(title, gp = grid::gpar(fontsize = fontsize))
  width <- grid::convertWidth(grid::grobWidth(label), "in", valueOnly = TRUE)
  page <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE)
  grid::grid.text(
    title,
    x = grid::unit(min(left, page - width - 0.04), "in"),
    y = grid::unit(1, "npc") - grid::unit(top_mm, "mm"),
    just = c("left", "top"), gp = grid::gpar(fontsize = fontsize)
  )
}

#' Draw a pheatmap gtable with a title over its colour bar.
#'
#' The heatmap is held `title_mm` short of the top to make the room.
draw_heatmap_with_legend_title <- function(gtable, title, title_mm = 5,
                                           fontsize = 10) {
  grid::pushViewport(grid::viewport(
    y = grid::unit(0, "npc"), just = "bottom",
    height = grid::unit(1, "npc") - grid::unit(title_mm, "mm")))
  grid::grid.draw(gtable)
  left <- heatmap_legend_left(gtable)
  grid::popViewport()
  draw_heatmap_legend_title(left, title, fontsize = fontsize)
}
