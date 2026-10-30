#!/usr/bin/env Rscript
# figures/landscape/setup/install_packages.R — install what config/ says a run needs.
#
# Reads config/packages.tsv (every package the figure scripts, fits and lib/
# load, with the version a panel needs where one is pinned), installs what is
# missing, and reports what is still off. It is the
# executable form of docs/INSTALL.md, and the Dockerfile runs it too, so a native
# install and the image are built by the same code.
#
#   Rscript figures/landscape/setup/install_packages.R                      # public set
#   Rscript figures/landscape/setup/install_packages.R --with-data          # plus the gated Data package
#   Rscript figures/landscape/setup/install_packages.R --data-source=../MotrpacHumanPreSuspensionData
#
# Options:
#   --with-data            install MotrpacHumanPreSuspensionData from GitHub. The
#                          repository is private: GITHUB_PAT must hold a token with
#                          read access to it (see docs/INSTALL.md).
#   --data-source=<dir>    install the Data package from a local checkout or an
#                          unzipped copy instead (what reviewers are sent).
#   --bioc=<release>       Bioconductor release to install from. Default: the
#                          release BiocManager already uses, or the one that pairs
#                          with this R (4.4 -> 3.20, 4.5 -> 3.22, 4.6 -> 3.23).
#   --stage=<stage>        cran-bioc | github | all (default). The Dockerfile runs
#                          the two stages as separate layers.
#   --root=<dir>           figures/landscape/, when not derivable from this file.
#   --dry-run              print what would be installed and stop.
#
# Packages install at their current release: CRAN and Bioconductor packages from
# the repositories, the GitHub packages from the default branch of each
# repository. Already-installed packages are left alone, whatever their
# version, except one below a ">=" floor in packages.tsv (ggpubr < 1.0.0),
# which is moved to the current release. Nothing else is upgraded, because the
# versions a panel was drawn from are recorded per PDF in
# outputs/panels/export_manifest.tsv and an unasked-for upgrade is the one way
# to make that record wrong without noticing.
#
# Exit status: non-zero when a package in packages.tsv is still missing, or
# still below its ">=" floor, at the end. A missing Data package is a note
# unless it was asked for.

options(warn = 1)

# ---- Arguments ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
flag <- function(name) paste0("--", name) %in% args
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}

known <- c("--with-data", "--dry-run", "--bioc=", "--stage=", "--root=", "--data-source=")
bad <- args[!vapply(args, function(a) any(startsWith(a, known)), logical(1))]
if (length(bad)) stop("unknown argument(s): ", paste(bad, collapse = " "))

stage <- opt("stage", "all")
if (!stage %in% c("cran-bioc", "github", "all")) stop("--stage must be cran-bioc, github or all")
with_data <- flag("with-data")
data_source <- opt("data-source")
dry_run <- flag("dry-run")

# The GitHub repositories, installed from their default branch.
ANALYSIS_REPO <- "MoTrPAC/MotrpacHumanPreSuspensionAnalysis"
DATA_REPO <- "MoTrPAC/MotrpacHumanPreSuspensionData"
BICQC_REPO <- "MoTrPAC/MotrpacBicQC"
PLIER_REPO <- "wgmao/PLIER"

# ---- Locate config/ ----------------------------------------------------------

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) return(NULL)
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
}
root <- opt("root", if (!is.null(script_dir())) file.path(script_dir(), "..") else getwd())
root <- normalizePath(root, mustWork = TRUE)
packages_tsv <- file.path(root, "config", "packages.tsv")
if (!file.exists(packages_tsv)) stop("not found: ", packages_tsv, "\n  pass --root=<path to figures/landscape>")

say <- function(...) cat("[install] ", paste0(...), "\n", sep = "")

installed <- function(p) requireNamespace(p, quietly = TRUE)
version_of <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)

# ---- R and Bioconductor ------------------------------------------------------

r_minor <- paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")
if (utils::compareVersion(paste(R.version$major, R.version$minor, sep = "."), "4.4.0") < 0) {
  stop("R ", R.version$major, ".", R.version$minor, " is too old: the Analysis package needs R >= 4.4")
}
if (r_minor != "4.4") {
  say("note: R ", r_minor, ". Every panel in assembled/ was drawn under R 4.4.1")
}

recommended_bioc <- switch(r_minor, "4.4" = "3.20", "4.5" = "3.22", "4.6" = "3.23", NA_character_)
if (!installed("BiocManager")) {
  say("installing BiocManager")
  if (!dry_run) utils::install.packages("BiocManager")
}
bioc <- opt("bioc")
if (is.null(bioc)) {
  current <- tryCatch(as.character(BiocManager::version()), error = function(e) NA_character_)
  # A release already in use that pairs with this R is kept: switching releases
  # reinstalls every Bioconductor package on the machine, which nobody asked for.
  valid <- switch(r_minor, "4.4" = c("3.19", "3.20"), "4.5" = c("3.21", "3.22"), "4.6" = "3.23", character(0))
  bioc <- if (!is.na(current) && current %in% valid) current else recommended_bioc
}
if (is.na(bioc)) stop("no Bioconductor release pairs with R ", r_minor, "; pass --bioc=")
say("R ", R.version$major, ".", R.version$minor, ", Bioconductor ", bioc)

options(Ncpus = max(1L, parallel::detectCores() - 1L))
# PLIER's tarball is ~156 MB and the Data package ~0.5 GB; the 60 s default
# aborts both on a slow link.
options(timeout = max(3600, getOption("timeout")))

# ---- What the manifest says --------------------------------------------------

# quote = "": the install column holds R calls with double quotes in them.
manifest <- read.csv(packages_tsv, sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
manifest$required_version <- trimws(as.character(manifest$required_version))
# required_version is empty (any release), ">= x" (a floor: below it the package
# is moved to the current release) or "== x" (the version the panels were drawn
# under: reported, never changed).
manifest$op <- ifelse(startsWith(manifest$required_version, ">="), ">=",
                      ifelse(startsWith(manifest$required_version, "=="), "==", NA_character_))
manifest$want <- trimws(sub("^(>=|==)", "", manifest$required_version))
manifest$github <- grepl("install_github", manifest$install, fixed = TRUE)

github_pkgs <- manifest$package[manifest$github]
cran_bioc_pkgs <- manifest$package[!manifest$github]
pins <- manifest[!is.na(manifest$op), , drop = FALSE]

pin_ok <- function(p, want, op) {
  got <- version_of(p)
  if (is.na(got)) return(FALSE)
  cmp <- utils::compareVersion(got, want)
  if (op == "==") cmp == 0 else if (op == ">=") cmp >= 0 else FALSE
}

# ---- Stage: cran-bioc --------------------------------------------------------
# One BiocManager::install() over the whole missing set. It resolves CRAN and
# Bioconductor (software and annotation) together against the one release, which
# is what keeps limma, edgeR, TMSig and their S4 neighbours on the same ABI.

if (stage %in% c("cran-bioc", "all")) {
  cran_bioc <- cran_bioc_pkgs
  missing <- cran_bioc[!vapply(cran_bioc, installed, logical(1))]
  if (!length(missing)) {
    say("every CRAN/Bioconductor package in packages.tsv is installed (", length(cran_bioc), ")")
  } else {
    say(length(missing), " to install from CRAN/Bioconductor ", bioc, ": ", paste(missing, collapse = " "))
    if (!dry_run) {
      options(BiocManager.check_repositories = FALSE)
      if (!identical(as.character(BiocManager::version()), bioc)) {
        BiocManager::install(version = bioc, ask = FALSE, update = FALSE)
      }
      BiocManager::install(missing, ask = FALSE, update = FALSE, version = bioc)
      still <- missing[!vapply(missing, installed, logical(1))]
      # install.packages() signals a failed compile with a warning, so check the
      # library rather than trusting the call. One retry, one package at a time,
      # because a transient download failure is indistinguishable from a real one.
      if (length(still)) {
        say("retrying ", length(still), " that did not install: ", paste(still, collapse = " "))
        for (p in still) tryCatch(BiocManager::install(p, ask = FALSE, update = FALSE, version = bioc),
                                  error = function(e) warning(p, ": ", conditionMessage(e), call. = FALSE))
      }
    }
  }
}

# ---- Stage: github -----------------------------------------------------------

install_gh <- function(ref, ...) {
  say("install_github(\"", ref, "\")")
  if (dry_run) return(invisible())
  # upgrade = "never": the default would move already-installed CRAN packages,
  # and nothing here asks for that.
  remotes::install_github(ref, upgrade = "never", build_vignettes = FALSE, ...)
}

if (stage %in% c("github", "all")) {
  if (!installed("remotes")) {
    say("installing remotes")
    if (!dry_run) utils::install.packages("remotes")
  }
  if (!nzchar(Sys.getenv("GITHUB_PAT")) && !nzchar(Sys.getenv("GITHUB_TOKEN"))) {
    say("note: GITHUB_PAT is unset. Anonymous GitHub access is rate-limited to 60 requests/hour,",
        " enough for these installs; the private Data package needs a token regardless")
  }

  # MotrpacBicQC: both data packages require >= 2.0.0, so a build below that is
  # replaced by the current one.
  bicqc <- pins[pins$package == "MotrpacBicQC", , drop = FALSE]
  if (nrow(bicqc) && !pin_ok("MotrpacBicQC", bicqc$want, bicqc$op)) {
    install_gh(BICQC_REPO)
  } else if (!installed("MotrpacBicQC")) {
    install_gh(BICQC_REPO)
  } else say("MotrpacBicQC ", version_of("MotrpacBicQC"), " already installed")

  if (!installed("PLIER")) install_gh(PLIER_REPO) else say("PLIER ", version_of("PLIER"), " already installed")

  if (installed("MotrpacHumanPreSuspensionAnalysis")) {
    say("MotrpacHumanPreSuspensionAnalysis ", version_of("MotrpacHumanPreSuspensionAnalysis"), " already installed")
  } else {
    install_gh(ANALYSIS_REPO)
  }

  if (!is.null(data_source)) {
    if (!file.exists(file.path(data_source, "DESCRIPTION"))) stop("--data-source has no DESCRIPTION: ", data_source)
    say("installing MotrpacHumanPreSuspensionData from ", data_source)
    if (!dry_run) remotes::install_local(data_source, upgrade = "never", build_vignettes = FALSE, force = TRUE)
  } else if (with_data) {
    if (installed("MotrpacHumanPreSuspensionData")) {
      say("MotrpacHumanPreSuspensionData ", version_of("MotrpacHumanPreSuspensionData"), " already installed")
    } else {
      if (!nzchar(Sys.getenv("GITHUB_PAT")) && !nzchar(Sys.getenv("GITHUB_TOKEN"))) {
        say("MotrpacHumanPreSuspensionData is a private repository and GITHUB_PAT is unset; this install will fail",
            " unless gitcreds holds a token")
      }
      install_gh(DATA_REPO)
    }
  } else if (!installed("MotrpacHumanPreSuspensionData")) {
    say("MotrpacHumanPreSuspensionData not installed and not requested (--with-data or --data-source=).",
        " 76 panels build without it; see docs/INSTALL.md for access")
  }
}

if (dry_run) { say("dry run; nothing installed"); quit(status = 0) }

# ---- Pins: move a CRAN package below its floor, then report every pin --------

if (stage == "all") {
  for (i in seq_len(nrow(pins))) {
    p <- pins$package[i]
    if (pins$op[i] != ">=" || pins$github[i]) next
    if (!installed(p) || pin_ok(p, pins$want[i], pins$op[i])) next
    say(p, " ", version_of(p), " is below its floor (", pins$op[i], " ", pins$want[i], "); installing the current CRAN release")
    utils::install.packages(p)
  }
}

cat("\n")
say("pinned versions in packages.tsv:")
problems <- character(0)
for (i in seq_len(nrow(pins))) {
  p <- pins$package[i]; got <- version_of(p)
  status <- if (is.na(got)) "MISSING" else if (pin_ok(p, pins$want[i], pins$op[i])) "ok" else "OFF"
  cat(sprintf("  %-36s %-8s installed %-8s required %s %s%s\n", p, status, ifelse(is.na(got), "-", got),
              pins$op[i], pins$want[i], if (pins$op[i] == ">=") "  (blocking)" else ""))
  if (status == "MISSING" && p == "MotrpacHumanPreSuspensionData" && !with_data && is.null(data_source)) next
  if (status == "MISSING" && stage == "cran-bioc" && pins$github[i]) next
  if (status == "MISSING" || (status == "OFF" && pins$op[i] == ">=")) problems <- c(problems, p)
}

check_set <- switch(stage, "cran-bioc" = cran_bioc_pkgs, github = github_pkgs, manifest$package)
if (!with_data && is.null(data_source)) check_set <- setdiff(check_set, "MotrpacHumanPreSuspensionData")
still_missing <- check_set[!vapply(check_set, installed, logical(1))]
if (length(still_missing)) {
  say("still missing from packages.tsv: ", paste(still_missing, collapse = " "))
  problems <- c(problems, still_missing)
}

if (length(problems)) {
  say("FAILED: ", length(unique(problems)), " package(s) are not where config/ needs them: ", paste(unique(problems), collapse = " "))
  quit(status = 1)
}
say("done. Next: bash figures/landscape/setup/preflight.sh")
quit(status = 0)
