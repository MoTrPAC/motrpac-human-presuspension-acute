#!/usr/bin/env bash
# figures/landscape/setup/preflight.sh — check a machine before the first run.
#
# The README's "Before the first run" table says what a run needs; this script
# checks it, and a few things the table leaves to the first failure. It follows
# motrpac-human-presuspension-repro's Stage 0 (scripts/00_preflight/00_preflight.sh):
# one row per check, a TSV report, and a non-zero exit only when a REQUIRED
# check fails. It changes nothing except the report file.
#
#   bash figures/landscape/setup/preflight.sh                # everything
#   bash figures/landscape/setup/preflight.sh --no-network   # skip the bucket and HTTPS probes
#
# Output: figures/landscape/outputs/preflight_report.tsv  (status<TAB>check<TAB>detail)
#         PREFLIGHT_REPORT=<path> writes it elsewhere.
# Exit:   1 if any check FAILs, 0 otherwise. WARN and INFO never change it.
#
# Statuses:
#   FAIL  nothing can run: R, a panel package, the Analysis package, a config
#         file, a figure script or a vendored input is missing.
#   WARN  some panels cannot build, or a number may differ from assembled/:
#         the Data package is absent, a pin is off, gsutil or bucket access is
#         absent, the machine is small.
#   INFO  a state the scripts handle by themselves: a fit not yet run, a gated
#         input absent. Reported so a SKIP later is not a surprise.
#
# Every variable config/landscape.env declares is honoured, so the checks look
# where the panels will look: EPIGEN_QC_DIR, PLIER_DIR, QUANT_BUCKET, GSUTIL and
# the fit overrides.

set -uo pipefail

# ---- Locate the tree --------------------------------------------------------
SETUP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LANDSCAPE_ROOT="${LANDSCAPE_ROOT:-$(cd "${SETUP_DIR}/.." && pwd)}"
export LANDSCAPE_ROOT
REPO_ROOT="$(cd "${LANDSCAPE_ROOT}/../.." && pwd)"

if [[ ! -f "${LANDSCAPE_ROOT}/lib/panel_export.R" || ! -f "${LANDSCAPE_ROOT}/config/panel_map.json" ]]; then
  printf '[FAIL] %s does not look like figures/landscape/ (no lib/panel_export.R and config/panel_map.json)\n' "${LANDSCAPE_ROOT}" >&2
  exit 1
fi

# The run-time variables, with the defaults the panels use. Values already
# exported win, exactly as they do for a panel.
# shellcheck source=../config/landscape.env
source "${LANDSCAPE_ROOT}/config/landscape.env"

NETWORK=1
for arg in "$@"; do
  case "${arg}" in
    --no-network) NETWORK=0 ;;
    -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) printf 'unknown argument: %s\n' "${arg}" >&2; exit 2 ;;
  esac
done

RSCRIPT="${RSCRIPT:-Rscript}"
GCLOUD="${GCLOUD:-gcloud}"
REPORT_TSV="${PREFLIGHT_REPORT:-${LANDSCAPE_ROOT}/outputs/preflight_report.tsv}"
mkdir -p "$(dirname "${REPORT_TSV}")"

# ---- Reporting --------------------------------------------------------------
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  _c_reset=$'\033[0m'; _c_red=$'\033[31m'; _c_grn=$'\033[32m'; _c_yel=$'\033[33m'; _c_blu=$'\033[34m'
else
  _c_reset=''; _c_red=''; _c_grn=''; _c_yel=''; _c_blu=''
fi

CHECK_FAILS=0; CHECK_WARNS=0; CHECK_INFOS=0; CHECK_PASSES=0
printf 'status\tcheck\tdetail\n' > "${REPORT_TSV}"

section() { printf '\n%s── %s%s\n' "${_c_blu}" "$*" "${_c_reset}"; }

# record_check STATUS CHECK DETAIL   (STATUS in PASS|WARN|FAIL|INFO)
record_check() {
  local status="$1" check="$2" detail="${3:-}"
  detail="${detail//$'\t'/ }"; detail="${detail//$'\n'/ }"
  printf '%s\t%s\t%s\n' "${status}" "${check}" "${detail}" >> "${REPORT_TSV}"
  case "${status}" in
    PASS) printf '%s[ OK ]%s %s — %s\n' "${_c_grn}" "${_c_reset}" "${check}" "${detail}"; CHECK_PASSES=$((CHECK_PASSES+1)) ;;
    WARN) printf '%s[WARN]%s %s — %s\n' "${_c_yel}" "${_c_reset}" "${check}" "${detail}"; CHECK_WARNS=$((CHECK_WARNS+1)) ;;
    FAIL) printf '%s[FAIL]%s %s — %s\n' "${_c_red}" "${_c_reset}" "${check}" "${detail}"; CHECK_FAILS=$((CHECK_FAILS+1)) ;;
    INFO) printf '%s[INFO]%s %s — %s\n' "${_c_blu}" "${_c_reset}" "${check}" "${detail}"; CHECK_INFOS=$((CHECK_INFOS+1)) ;;
  esac
}

have() { command -v "$1" >/dev/null 2>&1; }

# Run an R expression and keep only the lines it marks with "@@". A user's
# .Rprofile may print on startup, and a package may chatter on load; the marker
# is what separates the answer from that.
r_lines() {
  "${RSCRIPT}" -e "$1" "${@:2}" 2>/dev/null | sed -n 's/^@@//p'
}

url_reachable() {
  local url="$1"
  curl -fsS -m 20 -I "${url}" >/dev/null 2>&1 && return 0
  curl -fsS -m 20 -r 0-0 -o /dev/null "${url}" >/dev/null 2>&1
}

gcs_can_read() { "${GSUTIL}" ls "$1" >/dev/null 2>&1; }

rel() { printf '%s' "${1#"${REPO_ROOT}"/}"; }

printf 'Landscape preflight — %s\n' "${LANDSCAPE_ROOT}"
printf 'Report: %s\n' "${REPORT_TSV}"

# ---- 1. Tools ----------------------------------------------------------------
section "1. Tools"

HAVE_R=0
if have "${RSCRIPT}"; then
  HAVE_R=1
  record_check PASS "tool:Rscript" "found at $(command -v "${RSCRIPT}")"
  r_version="$(r_lines 'cat("@@", R.version$major, ".", R.version$minor, "\n", sep = "")')"
  r_version="${r_version// /}"
  case "${r_version}" in
    4.4.1) record_check PASS "r_version" "R 4.4.1, the version every panel in assembled/ was drawn under" ;;
    4.[4-9].*|4.[1-9][0-9].*|[5-9].*)
      record_check WARN "r_version" "R ${r_version}; assembled/ was drawn under 4.4.1. Panels build, but a number may differ in the last digit" ;;
    "") record_check FAIL "r_version" "could not run R to read its version" ;;
    *) record_check FAIL "r_version" "R ${r_version}; MotrpacHumanPreSuspensionAnalysis requires R >= 4.4" ;;
  esac
else
  record_check FAIL "tool:Rscript" "not on PATH; every script here runs under Rscript"
fi

if have "${GSUTIL}"; then
  record_check PASS "tool:gsutil" "found at $(command -v "${GSUTIL}")"
else
  record_check WARN "tool:gsutil" "not on PATH (GSUTIL=${GSUTIL}). FIG1C, ED1B, ED1D, SF1A-SF1F, ED1A and the sex and cell-type fits read the consortium bucket through it; every other panel runs without it"
fi
if have "${GCLOUD}"; then
  record_check PASS "tool:gcloud" "found at $(command -v "${GCLOUD}")"
else
  record_check WARN "tool:gcloud" "not on PATH; needed once, to authenticate gsutil (gcloud auth login)"
fi
if have curl; then
  record_check PASS "tool:curl" "found at $(command -v curl)"
else
  record_check WARN "tool:curl" "not on PATH; only this script's HTTPS probe uses it (R downloads the epigenomics DA itself)"
fi

# ---- 2. The packages in config/packages.tsv ---------------------------------
# Every row but the two data packages, which section 3 reports with the
# message their absence deserves.
section "2. R packages (config/packages.tsv)"

PACKAGES_TSV="${LANDSCAPE_ROOT}/config/packages.tsv"
if [[ ! -f "${PACKAGES_TSV}" ]]; then
  record_check FAIL "packages" "config/packages.tsv not found"
elif (( HAVE_R )); then
  # A while loop, not mapfile: macOS ships bash 3.2, which has no mapfile.
  PANEL_PKGS=()
  while IFS= read -r pkg; do
    [[ -n "${pkg}" ]] && PANEL_PKGS+=("${pkg}")
  done < <(awk -F'\t' '/^#/ || $1 == "package" || $1 ~ /^MotrpacHumanPreSuspension/ {next} {print $1}' "${PACKAGES_TSV}")
  missing_pkgs="$(r_lines '
    pkgs <- commandArgs(trailingOnly = TRUE)
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    cat("@@", paste(miss, collapse = " "), "\n", sep = "")
  ' "${PANEL_PKGS[@]}")"
  if [[ -z "${missing_pkgs// /}" ]]; then
    record_check PASS "packages" "all ${#PANEL_PKGS[@]} packages installed"
  else
    record_check FAIL "packages" "missing: ${missing_pkgs# } — run: Rscript figures/landscape/setup/install_packages.R"
  fi
fi

# ---- 3. The pinned rows of config/packages.tsv ------------------------------
# The rows with a required_version: the two data packages, MotrpacBicQC, ggpubr
# and rstatix. A ">=" floor that is not met is a FAIL because that version
# decides numbers, not appearance (ggpubr 0.6.3 and 1.0.0 draw different FIG4G
# and ED6E bubbles). An "==" row that is off is a WARN: the panel builds, and
# export_manifest.tsv records what it was drawn from. The Data package is gated,
# so its absence is a WARN that names what it costs rather than a FAIL.
section "3. Pinned versions (config/packages.tsv)"

DATA_PRESENT=0; ANALYSIS_PRESENT=0
if [[ -f "${PACKAGES_TSV}" ]] && (( HAVE_R )); then
  while IFS=$'\t' read -r pkg status installed want op blocking; do
    [[ -z "${pkg}" ]] && continue
    case "${pkg}:${status}" in
      *:OK)
        record_check PASS "pin:${pkg}" "${installed} (${op} ${want})"
        [[ "${pkg}" == MotrpacHumanPreSuspensionData ]] && DATA_PRESENT=1
        [[ "${pkg}" == MotrpacHumanPreSuspensionAnalysis ]] && ANALYSIS_PRESENT=1 ;;
      MotrpacHumanPreSuspensionData:MISSING)
        record_check WARN "pin:${pkg}" "not installed. It is access-gated (consortium data-access request). 17 panels, the four fits and their tables need it; 76 panels build from MotrpacHumanPreSuspensionAnalysis alone" ;;
      MotrpacHumanPreSuspensionAnalysis:MISSING)
        record_check FAIL "pin:${pkg}" "not installed; every panel reads it. See docs/INSTALL.md" ;;
      *:MISSING)
        record_check FAIL "pin:${pkg}" "not installed (want ${op} ${want})" ;;
      *:MISMATCH)
        [[ "${pkg}" == MotrpacHumanPreSuspensionData ]] && DATA_PRESENT=1
        [[ "${pkg}" == MotrpacHumanPreSuspensionAnalysis ]] && ANALYSIS_PRESENT=1
        if [[ "${blocking}" == "yes" ]]; then
          record_check FAIL "pin:${pkg}" "installed ${installed}, need ${op} ${want} (blocking: this version decides numbers)"
        else
          record_check WARN "pin:${pkg}" "installed ${installed}, pinned ${op} ${want}; panels build, export_manifest.tsv records the version used"
        fi ;;
      *) record_check WARN "pin:${pkg}" "unrecognised status '${status}'" ;;
    esac
  done < <(r_lines '
    tsv <- commandArgs(trailingOnly = TRUE)[1]
    rows <- read.csv(tsv, sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
    req <- trimws(as.character(rows$required_version))
    rows <- rows[nzchar(req), , drop = FALSE]; req <- req[nzchar(req)]
    for (i in seq_len(nrow(rows))) {
      p <- rows$package[i]; op <- substr(req[i], 1, 2); want <- trimws(substring(req[i], 3))
      got <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
      status <- if (is.na(got)) "MISSING" else {
        cmp <- utils::compareVersion(got, want)
        ok <- if (op == "==") cmp == 0 else if (op == ">=") cmp >= 0 else NA
        if (isTRUE(ok)) "OK" else "MISMATCH"
      }
      blocking <- if (op == ">=") "yes" else "no"
      cat("@@", paste(p, status, ifelse(is.na(got), "-", got), want, op, blocking, sep = "\t"), "\n", sep = "")
    }
  ' "${PACKAGES_TSV}")
fi

# ---- 4. The packages load and are the build the panels expect ---------------
# Attaching is the test, not requireNamespace(): the figure scripts attach both
# data packages (lib/panel_export.R, attach_data_packages), and loading
# MotrpacHumanPreSuspensionAnalysis pulls Mfuzz and through it tcltk, which is
# where a macOS install without XQuartz fails. The two formals checks prove the
# installed build is the one the pin names in content: Analysis 2.0.8 added
# plot_enrich_heatmap(return_drawing = ), which FIG4A, FIG4B and ED5A call, and
# Data 2.0.3's load_qc(epigen = ) is what FIG1C, ED1B, ED1D and SF1 call.
section "4. Package load"

if (( HAVE_R )); then
  while IFS=$'\t' read -r key value; do
    [[ -z "${key}" ]] && continue
    case "${key}:${value}" in
      analysis_attach:ok) record_check PASS "load:MotrpacHumanPreSuspensionAnalysis" "attaches" ;;
      analysis_attach:*) record_check FAIL "load:MotrpacHumanPreSuspensionAnalysis" "does not attach: ${value}" ;;
      analysis_objects:*) record_check PASS "load:analysis_objects" "${value} data objects (FCM_ORA, COVARIATES_FILE, the *_DA and *_SUM_STATS tables)" ;;
      analysis_api:ok) record_check PASS "load:analysis_api" "plot_enrich_heatmap() takes return_drawing (Analysis >= 2.0.8)" ;;
      analysis_api:*) record_check FAIL "load:analysis_api" "plot_enrich_heatmap() has no return_drawing argument; FIG4A, FIG4B and ED5A need Analysis >= 2.0.8" ;;
      data_attach:absent) : ;;
      data_attach:ok) record_check PASS "load:MotrpacHumanPreSuspensionData" "attaches" ;;
      data_attach:*) record_check FAIL "load:MotrpacHumanPreSuspensionData" "installed but does not attach: ${value}" ;;
      data_api:ok) record_check PASS "load:data_api" "load_qc() takes epigen (Data >= 2.0.3)" ;;
      data_api:absent) : ;;
      data_api:*) record_check WARN "load:data_api" "load_qc() has no epigen argument; FIG1C, ED1B, ED1D and SF1 need Data >= 2.0.3" ;;
      bicqc:ok) record_check PASS "load:MotrpacBicQC" "loads; assay_codes present (FIG4, ED5)" ;;
      bicqc:*) record_check FAIL "load:MotrpacBicQC" "${value}" ;;
      tcltk:TRUE) record_check PASS "capability:tcltk" "available (Mfuzz, imported by the Analysis package, needs it)" ;;
      tcltk:*) record_check WARN "capability:tcltk" "R was built without Tcl/Tk support; on macOS install XQuartz (see docs/INSTALL.md)" ;;
      cairo:TRUE) record_check PASS "capability:cairo" "available (cairo_pdf for the enrichment heatmaps)" ;;
      cairo:*) record_check WARN "capability:cairo" "no cairo device; plot_enrich_heatmap() writes PDFs with cairo_pdf(): FIG4A, FIG4B, ED5A" ;;
      bioc:*) record_check "${value%%|*}" "bioc_version" "${value#*|}" ;;
    esac
  done < <(r_lines '
    say <- function(k, v) cat("@@", k, "\t", v, "\n", sep = "")
    att <- function(p) tryCatch({
      suppressPackageStartupMessages(library(p, character.only = TRUE, warn.conflicts = FALSE)); "ok"
    }, error = function(e) gsub("[\t\n]+", " ", conditionMessage(e)))
    a <- att("MotrpacHumanPreSuspensionAnalysis"); say("analysis_attach", a)
    if (a == "ok") {
      say("analysis_objects", nrow(utils::data(package = "MotrpacHumanPreSuspensionAnalysis")$results))
      say("analysis_api", if ("return_drawing" %in% names(formals(MotrpacHumanPreSuspensionAnalysis::plot_enrich_heatmap))) "ok" else "missing")
    }
    if (requireNamespace("MotrpacHumanPreSuspensionData", quietly = TRUE)) {
      d <- att("MotrpacHumanPreSuspensionData"); say("data_attach", d)
      if (d == "ok") say("data_api", if ("epigen" %in% names(formals(MotrpacHumanPreSuspensionData::load_qc))) "ok" else "missing")
    } else { say("data_attach", "absent"); say("data_api", "absent") }
    say("bicqc", tryCatch({
      if (!requireNamespace("MotrpacBicQC", quietly = TRUE)) stop("not installed")
      if (!"assay_codes" %in% utils::data(package = "MotrpacBicQC")$results[, "Item"]) stop("loads but has no assay_codes object")
      "ok"
    }, error = function(e) conditionMessage(e)))
    say("tcltk", capabilities("tcltk")); say("cairo", capabilities("cairo"))
    rv <- paste(R.version$major, sub("[.].*$", "", R.version$minor), sep = ".")
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      bv <- as.character(BiocManager::version())
      ok <- switch(rv, "4.4" = c("3.19", "3.20"), "4.5" = c("3.21", "3.22"), "4.6" = "3.23", character(0))
      if (bv %in% ok) say("bioc", paste0("PASS|Bioconductor ", bv, " on R ", rv))
      else say("bioc", paste0("WARN|Bioconductor ", bv, " on R ", rv, "; R 4.4 pairs with 3.19 or 3.20, 4.5 with 3.21 or 3.22, 4.6 with 3.23. Mixed releases are where limma, edgeR and TMSig stop loading together"))
    } else say("bioc", "WARN|BiocManager is not installed, so the Bioconductor release cannot be checked")
  ')
fi

# ---- 5. Configuration and the scripts it names -----------------------------
section "5. Configuration"

for f in config/panel_map.json config/table_map.json config/highlights.json; do
  if [[ ! -f "${LANDSCAPE_ROOT}/${f}" ]]; then
    record_check FAIL "config:${f##*/}" "missing"
  elif (( HAVE_R )); then
    parse="$(r_lines 'x <- tryCatch({ jsonlite::read_json(commandArgs(TRUE)[1]); "ok" }, error = function(e) conditionMessage(e)); cat("@@", x, "\n", sep = "")' "${LANDSCAPE_ROOT}/${f}")"
    [[ "${parse}" == "ok" ]] && record_check PASS "config:${f##*/}" "parses" || record_check FAIL "config:${f##*/}" "does not parse: ${parse}"
  fi
done

for f in config/landscape.env extended_data_3/sex_da.env extended_data_3/celltype_da.env \
         figure_3/clinical_omics.env figure_4/plier.env figure_4/plier_lvs.env figure_4/cmeans.env; do
  [[ -f "${LANDSCAPE_ROOT}/${f}" ]] \
    && record_check PASS "env:${f}" "present" \
    || record_check FAIL "env:${f}" "missing; the script that sources it stops"
done

# Every script the two manifests name must exist: run_all.R builds its list from
# them, and a figure whose script is gone fails in setup rather than per panel.
if (( HAVE_R )) && [[ -f "${LANDSCAPE_ROOT}/config/panel_map.json" ]]; then
  while IFS=$'\t' read -r key value; do
    case "${key}" in
      scripts_ok) record_check PASS "manifest:scripts" "${value}" ;;
      scripts_missing) record_check FAIL "manifest:scripts" "named by a manifest but not on disk: ${value}" ;;
      panels) record_check PASS "manifest:panels" "${value}" ;;
    esac
  done < <(r_lines '
    root <- commandArgs(TRUE)[1]
    pm <- jsonlite::read_json(file.path(root, "config", "panel_map.json"), simplifyVector = FALSE)
    tm <- jsonlite::read_json(file.path(root, "config", "table_map.json"), simplifyVector = FALSE)
    scripts <- unique(c(vapply(pm$panels, function(p) p$script, character(1)),
                        unlist(lapply(tm$tables, function(t) t$script))))
    miss <- scripts[!file.exists(file.path(root, scripts))]
    cat("@@panels\t", length(pm$panels), " panels in panel_map.json, ", length(tm$tables), " tables in table_map.json\n", sep = "")
    if (length(miss)) cat("@@scripts_missing\t", paste(miss, collapse = " "), "\n", sep = "")
    else cat("@@scripts_ok\tall ", length(scripts), " scripts the manifests name are on disk\n", sep = "")
  ' "${LANDSCAPE_ROOT}")
fi

# ---- 6. Inputs vendored under <figure>/sources/ ------------------------------
# docs/external_dependencies.md is the list. Every file here ships with the
# repository, so a missing one is an incomplete checkout, not a configuration
# problem. The two gated inputs are INFO rows: their absence is the normal
# public state and the panels that read them say so and skip.
section "6. Vendored inputs"

VENDORED=(
  figure_1/sources/1kg_samples.tsv
  figure_1/sources/1kg_population_colors.csv
  figure_1/sources/motrpac_named_features.csv
  figure_2/sources/landscape_figureED2_PTMSigDB_ORA_results_v2_Aug2026.csv
  figure_5/sources/tfproanno.RDS
  figure_6/sources/Gambardella-et-al-2020-table-s2.xlsx
  figure_6/sources/TFEB_TARGET_GENES.v2024.1.Hs.tsv
  figure_6/sources/GSM2354032_TFEB_e5_peaks.bed
  figure_6/sources/TFEB_targets_v2.xlsx
  figure_6/sources/SCION_supplemental_table_edges.csv
  figure_6/sources/Merged_Network_combined_NMS_v2.csv
  figure_7/sources/compartments_extracellular.tsv
)
missing_vendored=()
for f in "${VENDORED[@]}"; do
  [[ -s "${LANDSCAPE_ROOT}/${f}" ]] || missing_vendored+=("${f}")
done
if (( ${#missing_vendored[@]} == 0 )); then
  record_check PASS "vendored:files" "all ${#VENDORED[@]} tracked inputs present"
else
  record_check FAIL "vendored:files" "missing from the checkout: ${missing_vendored[*]}"
fi

homer_dir="${LANDSCAPE_ROOT}/figure_5/sources/Precovid_DEG_HOMER_KnownTF_Results"
homer_n="$(find "${homer_dir}" -name '*_knownResults.txt' 2>/dev/null | wc -l | tr -d ' ')"
if [[ "${homer_n}" == "22" ]]; then
  record_check PASS "vendored:homer" "22 knownResults.txt files (FIG5B-FIG5H, ED7, ST5a)"
else
  record_check FAIL "vendored:homer" "expected 22 *_knownResults.txt under figure_5/sources/Precovid_DEG_HOMER_KnownTF_Results, found ${homer_n}"
fi

gated_input() {
  local check="$1" default_rel="$2" panels="$3"; shift 3
  local path="" var
  for var in "$@"; do
    [[ -n "${!var:-}" ]] && { path="${!var}"; break; }
  done
  [[ -z "${path}" ]] && path="${LANDSCAPE_ROOT}/${default_rel}"
  if [[ -s "${path}" ]]; then
    record_check INFO "gated:${check}" "present at $(rel "${path}"); ${panels} can be drawn. Keep it out of git"
  else
    record_check INFO "gated:${check}" "absent (individual-level, not distributed): ${panels} skip with CANNOT BE RUN BY THE PUBLIC. Override: $*"
  fi
}
gated_input ED1A_pca "figure_1/sources/1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS" "ED1A" ED1A_PCA_RDS
gated_input CIBERSORTx "extended_data_3/sources/CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv" "ED3D, ED3E and celltype_da_fit.R" ED3D_CIBERSORTX_CSV CELLTYPE_DA_CIBERSORTX_CSV

# ---- 7. The four fits --------------------------------------------------------
# Not required: a panel whose fit is absent reports SKIPPED and names the fit.
# Listed so the SKIP is expected. The manifest TSV beside each fit is what marks
# it complete, so that is the file looked for, at the path the fit's .env (or
# an exported override) names.
section "7. Fits (cached computations)"

fit_status() {
  local check="$1" manifest="$2" panels="$3" cmd="$4"
  if [[ -s "${manifest}" ]]; then
    record_check INFO "fit:${check}" "present ($(rel "${manifest}")); ${panels} build"
  else
    record_check INFO "fit:${check}" "not run; ${panels} skip until: ${cmd}"
  fi
}
fit_status sex_da "${SEX_DA_MANIFEST_TSV:-${SEX_STRATIFIED_DA_DIR}/sex_da_manifest.tsv}" \
  "ED3A, ST2g, ST2h" "Rscript figures/landscape/extended_data_3/sex_da_fit.R (hours; Data package and bucket access)"
fit_status celltype_da "${CELLTYPE_DA_MANIFEST_TSV:-${CELLTYPE_DA_DIR}/celltype_da_manifest.tsv}" \
  "ED3E" "Rscript figures/landscape/extended_data_3/celltype_da_fit.R (needs the gated CIBERSORTx input)"
fit_status clinical_omics "${CLINICAL_OMICS_MANIFEST_TSV:-$(dirname "${CLINICAL_OMICS_TABLE}")/clinical_omics_manifest.tsv}" \
  "FIG3EFG, ED4G, ST3f, ST3g, ST3h, EDT1" "Rscript figures/landscape/figure_3/clinical_omics_fit.R (tens of minutes; Data package)"
fit_status plier "${PLIER_MANIFEST_TSV:-${PLIER_DIR}/plier_manifest.tsv}" \
  "FIG4G, FIG4H, ED6A-ED6I" "Rscript figures/landscape/figure_4/plier_fit.R (~2 h, ~8.5 GB; Data package)"

# ---- 8. Consortium access ----------------------------------------------------
# Public panels never touch a bucket. These rows say which gated reads will
# work, so a WARN here is a scope statement, not a broken install.
section "8. Consortium access (gcloud, buckets)"

EPIGEN_COLLECTION="gs://motrpac-data-hub/analysis/human-precovid-sed-adu/c2.0"
if have "${GCLOUD}"; then
  active_acct="$("${GCLOUD}" auth list --filter=status:ACTIVE --format='value(account)' 2>/dev/null | head -n1)"
  if [[ -n "${active_acct}" ]]; then
    record_check PASS "gcloud_account" "active: ${active_acct}"
  else
    record_check WARN "gcloud_account" "no active account; run: gcloud auth login"
  fi
fi

epigen_cached="$(find "${EPIGEN_QC_DIR}" -type f 2>/dev/null | wc -l | tr -d ' ')"
if (( NETWORK )) && have "${GSUTIL}"; then
  if gcs_can_read "${EPIGEN_COLLECTION}/"; then
    record_check PASS "bucket:epigen_qc" "readable: ${EPIGEN_COLLECTION} (FIG1C, ED1B, ED1D, SF1 cache it under $(rel "${EPIGEN_QC_DIR}"))"
  elif (( epigen_cached > 0 )); then
    record_check INFO "bucket:epigen_qc" "not readable, but ${epigen_cached} files are already cached under $(rel "${EPIGEN_QC_DIR}"); panels read the cache"
  else
    record_check WARN "bucket:epigen_qc" "cannot list ${EPIGEN_COLLECTION}: FIG1C, ED1B, ED1D and SF1A-SF1F need consortium read access (or an EPIGEN_QC_DIR filled elsewhere)"
  fi
  if gcs_can_read "${QUANT_BUCKET}/"; then
    record_check PASS "bucket:quant" "readable: ${QUANT_BUCKET} (SF1C-SF1F, the sex and cell-type fits)"
  else
    record_check WARN "bucket:quant" "cannot list ${QUANT_BUCKET}: SF1C-SF1F and the two DA fits read raw counts from it"
  fi
elif (( NETWORK )); then
  record_check WARN "bucket:epigen_qc" "no gsutil, so bucket access was not probed"
else
  record_check INFO "bucket:epigen_qc" "skipped (--no-network); ${epigen_cached} files cached under $(rel "${EPIGEN_QC_DIR}")"
fi

# ---- 9. Public downloads -----------------------------------------------------
# The epigenomics DA tables (~7.7 GB over five files) come from the public c2.0
# CloudFront release on every call that asks for them, through
# load_differential_analysis(epigen = TRUE). No credentials; a HEAD request on
# the smallest file is the probe.
section "9. Public downloads"

EPIGEN_DA_URL="https://d1yw74buhe0ts0.cloudfront.net/data/analysis/human_presuspension_sed_adu/c2.0/epigenomics/da/human-precovid-sed-adu_t11-adipose_epigen-methylcap-seq_da_malax-glmm-acute_v1.2.txt"
if (( ! NETWORK )); then
  record_check INFO "https:epigen_da" "skipped (--no-network)"
elif ! have curl; then
  record_check WARN "https:epigen_da" "no curl, so the CloudFront release was not probed"
elif url_reachable "${EPIGEN_DA_URL}"; then
  record_check PASS "https:epigen_da" "public epigenomics DA release reachable (FIG2, ED2, ED4A, SF1G, ST2d download from it)"
else
  record_check WARN "https:epigen_da" "cannot reach ${EPIGEN_DA_URL%/*}/ — FIG2A-FIG2F, ED2, ED4A, SF1G and ST2d download the epigenomics DA from it on every run"
fi

# ---- 10. Machine -------------------------------------------------------------
# What the expensive paths need. The epigen QC cache is ~6 GB, the two sex and
# cell-type tables ~0.6 GB, and the PLIER RNA fit holds ~8.5 GB in memory.
section "10. Machine"

free_kb="$(df -Pk "${LANDSCAPE_ROOT}" 2>/dev/null | awk 'NR==2 {print $4}')"
if [[ -n "${free_kb}" ]]; then
  free_gb=$(( free_kb / 1024 / 1024 ))
  if (( free_gb >= 25 )); then
    record_check PASS "disk" "${free_gb} GB free under $(rel "${LANDSCAPE_ROOT}")"
  else
    record_check WARN "disk" "${free_gb} GB free under $(rel "${LANDSCAPE_ROOT}"); the epigen QC cache alone is ~6 GB and a full run writes several GB under outputs/"
  fi
fi

mem_kb=""
if [[ -r /proc/meminfo ]]; then
  mem_kb="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
elif have sysctl; then
  mem_b="$(sysctl -n hw.memsize 2>/dev/null || true)"
  [[ -n "${mem_b}" ]] && mem_kb=$(( mem_b / 1024 ))
fi
if [[ -n "${mem_kb}" ]]; then
  mem_gb=$(( mem_kb / 1024 / 1024 ))
  if (( mem_gb >= 16 )); then
    record_check PASS "memory" "${mem_gb} GB"
  else
    record_check WARN "memory" "${mem_gb} GB; the PLIER fit needs ~8.5 GB and FIG2 holds several GB of differential analysis. In Docker, raise the VM's memory limit"
  fi
fi

# ---- Summary -----------------------------------------------------------------
echo
printf 'Preflight: %d passed, %d failed, %d warnings, %d notes. Report: %s\n' \
  "${CHECK_PASSES}" "${CHECK_FAILS}" "${CHECK_WARNS}" "${CHECK_INFOS}" "$(rel "${REPORT_TSV}")"
if (( CHECK_FAILS > 0 )); then
  printf '%s[FAIL]%s %d required check(s) failed. Fix those before running a figure script.\n' "${_c_red}" "${_c_reset}" "${CHECK_FAILS}" >&2
  exit 1
fi
if (( CHECK_WARNS > 0 )); then
  printf '%s[ OK ]%s Ready for the panels the warnings above do not name.\n' "${_c_grn}" "${_c_reset}"
else
  printf '%s[ OK ]%s Ready.\n' "${_c_grn}" "${_c_reset}"
fi
exit 0
