# Installing the landscape figure environment

What `figures/landscape/` needs before its first run, how to install it, how to check
the result, and the Docker alternative. Paths here are relative to `figures/landscape/`, as
in `external_dependencies.md` beside this file; the commands are given from the repository
root. The package list lives in two files under `config/` and everything here reads them
rather than restating them:

| File | Holds |
|---|---|
| `config/panel_packages.txt` | every package the figure scripts, the four fits and `lib/` load: 39 from CRAN and Bioconductor, plus `MotrpacBicQC` and `PLIER` from GitHub. Each line names the package's source and the one-line R call that installs it on its own |
| `config/required_packages.tsv` | the five pins whose version changes a number, not just an appearance: the two data packages, `MotrpacBicQC`, `ggpubr >= 1.0.0` and its `rstatix`. Each row names its source (CRAN, Bioconductor or GitHub) and the one-line R call that installs it on its own |

Three commands, run from the repository root:

```bash
Rscript figures/landscape/setup/install_packages.R     # install what the two files name
bash    figures/landscape/setup/preflight.sh           # check the machine, panel by panel
Rscript figures/landscape/figure_2/FIG2.R              # build something
```

## 1. What you need

| | | Who needs it |
|---|---|---|
| **R 4.4.1** | the version every panel in `assembled/` was drawn under. R >= 4.4 is a hard floor (the Analysis package requires it); a later 4.x builds every panel but a number may differ in its last digit | everyone |
| **Bioconductor 3.20** | the release that pairs with R 4.4 (3.19 also does; 4.5 pairs with 3.22, 4.6 with 3.23). Mixing releases is where `limma`, `edgeR` and `TMSig` stop loading together | everyone |
| **System libraries** | macOS: Xcode command-line tools and **XQuartz**. Linux: the development headers listed in §2 | everyone |
| **`MotrpacHumanPreSuspensionAnalysis` 2.0.8** | public, from GitHub. The aggregate results every panel reads, and the palettes and labels | everyone |
| **`MotrpacBicQC` >= 2.0.0** | public, from GitHub. Both data packages require it | everyone |
| **`MotrpacHumanPreSuspensionData` 2.0.3** | access-gated: a formal MoTrPAC consortium data-access request, then read access to a private GitHub repository. Of the 94 panels, 76 build without it; 17 panels, the four fits and their tables need it | the fits, FIG1C, ED1B, ED1D, SF1, and the single-feature and sex-stratified panels |
| **`gsutil`, authenticated** | the Google Cloud SDK, logged in with an account that has consortium read access to the data hub | FIG1C, ED1B, ED1D, SF1A-SF1F, ED1A, and the sex and cell-type fits |
| **Disk and memory** | ~6 GB for the epigenomics QC cache the first time FIG1C or ED1 runs; ~8.5 GB of memory for the PLIER fit | those panels and that fit |

The epigenomics differential-analysis tables (~7.7 GB over five files, read by FIG2, ED2,
ED4A, SF1G and ST2d) are a public HTTPS download and need neither `gsutil` nor credentials,
only a network connection. They are not cached: every call downloads what it asks for.

## 2. System setup

### macOS

```bash
xcode-select --install            # compilers, for the packages that build from source
brew install --cask xquartz       # then log out and back in
```

XQuartz is not optional. The macOS builds of R link `tcltk` and the cairo device against the
X11 libraries it provides; `MotrpacHumanPreSuspensionAnalysis` imports `Mfuzz`, which imports
`tcltk`, so without XQuartz the Analysis package fails to install with

```
unable to load shared object '.../library/tcltk/libs/tcltk.so':
  Library not loaded: /opt/X11/lib/libX11.6.dylib
```

and FIG4A, FIG4B and ED5A, which write their heatmaps through `cairo_pdf()`, fail to draw.

If `gsutil` is on your shell's `PATH` but a script cannot find it, set `GSUTIL` to the full
path (`config/landscape.env` reads it); RStudio launched from the Dock does not inherit the
shell's `PATH`, which is one reason the scripts are run with `Rscript` from a terminal.

### Linux (Ubuntu 22.04 names)

```bash
sudo apt-get install -y \
  libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libicu-dev libsqlite3-dev \
  zlib1g-dev libbz2-dev liblzma-dev libcairo2-dev libx11-dev libxt-dev \
  libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
  libpng-dev libjpeg-dev libtiff-dev libwebp-dev libglpk-dev libgmp-dev libfftw3-dev \
  tcl-dev tk-dev
```

This is the list `docker/Dockerfile` installs (there as the runtime `libtcl8.6`/`libtk8.6`,
since the rocker image's R is already built against them), which is the one place it has been
exercised end to end. Other distributions: the same libraries under their own names.

### Google Cloud SDK

Only for the panels and fits the table in §1 names. Install from
<https://cloud.google.com/sdk/docs/install>, then authenticate once on the machine:

```bash
gcloud auth login
gsutil ls gs://motrpac-data-hub/analysis/human-precovid-sed-adu/c2.0/    # should list, not 403
```

A `403` or `AccessDenied` is access, not installation: the data-hub buckets are consortium
resources. The preflight reports it as a WARN naming the panels it costs.

## 3. Install the R packages

```bash
Rscript figures/landscape/setup/install_packages.R
```

From a clean R this installs, in order: `BiocManager` and the Bioconductor release that pairs
with your R; every missing CRAN and Bioconductor package in `config/panel_packages.txt` in one
resolver pass; and `MotrpacBicQC`, `PLIER` and `MotrpacHumanPreSuspensionAnalysis` from the
current release on GitHub. It then prints every row of `config/required_packages.tsv` with the
installed version beside the version the panels were drawn under, and exits non-zero if
anything a panel loads is still missing.

What it does not do: upgrade anything already installed. The versions each PDF was drawn from
are recorded per panel in `outputs/panels/export_manifest.tsv`, and an unasked-for upgrade is
the one way to make that record wrong without noticing. The single exception is `ggpubr`:
below 1.0.0 it adjusts p-values across a different set of rows, FIG4G and ED6E move, and so
it is the one blocking pin the installer corrects.

`--dry-run` prints the plan and installs nothing. `--bioc=3.19` keeps an existing 3.19
library on R 4.4 rather than moving it to 3.20 (the installer already does this when 3.19 is
what `BiocManager` reports).

### The Data package

Two routes, once your data-access request is approved:

```bash
# (a) from GitHub, with a personal access token that can read the private repository
GITHUB_PAT=ghp_... Rscript figures/landscape/setup/install_packages.R --with-data

# (b) from a checkout or an unzipped copy, which is what reviewers are sent
Rscript figures/landscape/setup/install_packages.R --data-source=../MotrpacHumanPreSuspensionData
```

The package is about 0.5 GB to download and needs roughly 2 GB of free memory to install.
Without it the installer says so and carries on; `preflight.sh` reports it as a WARN; and a
panel that needs it fails with a message naming the package rather than something later.

### By hand, if you would rather

The installer is a convenience over these calls. In R 4.4:

```r
install.packages("BiocManager")
BiocManager::install(version = "3.20")
pkgs <- readLines("figures/landscape/config/panel_packages.txt")
pkgs <- trimws(sub("#.*$", "", pkgs)); pkgs <- setdiff(pkgs[nzchar(pkgs)], c("MotrpacBicQC", "PLIER"))
BiocManager::install(pkgs, update = FALSE, ask = FALSE)
install.packages("remotes")
remotes::install_github("MoTrPAC/MotrpacBicQC", upgrade = "never")
remotes::install_github("wgmao/PLIER", upgrade = "never")
remotes::install_github("MoTrPAC/MotrpacHumanPreSuspensionAnalysis", upgrade = "never")
# with access:
remotes::install_github("MoTrPAC/MotrpacHumanPreSuspensionData", upgrade = "never")
```

`pak::pak("MoTrPAC/MotrpacHumanPreSuspensionAnalysis")` also works, as each package's own
README says. Either way the package lands at its current release; `config/required_packages.tsv`
records the version each panel was drawn under, and `preflight.sh` reports the difference.

## 4. Check the machine: `setup/preflight.sh`

```bash
bash figures/landscape/setup/preflight.sh                 # everything
bash figures/landscape/setup/preflight.sh --no-network    # no bucket or HTTPS probes
```

Modelled on the Stage 0 preflight of `motrpac-human-presuspension-repro`: one line per check,
a TSV report at `outputs/preflight_report.tsv` (`status`, `check`, `detail`), and a non-zero
exit only when something **required** is missing. It changes nothing but that report. Ten
sections:

| | Checks | FAILs when |
|---|---|---|
| 1 | `Rscript`, R version, `gsutil`, `gcloud`, `curl` on `PATH` | no R, or R < 4.4 |
| 2 | every package in `config/panel_packages.txt` is installed | any is missing |
| 3 | every row of `config/required_packages.tsv` against its pin | the Analysis package is missing, or a **blocking** pin is off |
| 4 | the data packages attach; the installed builds carry the API the panels call (`plot_enrich_heatmap(return_drawing=)`, `load_qc(epigen=)`); `MotrpacBicQC` loads; Tcl/Tk and cairo capabilities; Bioconductor release pairs with R | a package that is installed does not attach, or the Analysis build predates 2.0.8 |
| 5 | the three `config/*.json` parse, the seven `.env` files exist, every script the manifests name is on disk | any is missing |
| 6 | the twelve tracked inputs under `*/sources/` and the 22 HOMER files; whether the two gated inputs are present | a tracked input is missing (an incomplete checkout) |
| 7 | which of the four fits have run | never: a missing fit is INFO, and the panels that read it SKIP |
| 8 | an active `gcloud` account; read access to the epigenomics collection and the quantification bucket | never: this is access, reported as WARN with the panels it costs |
| 9 | the public CloudFront release of the epigenomics DA tables answers | never (WARN) |
| 10 | free disk under `figures/landscape/`, total memory | never (WARN under 25 GB / 16 GB) |

Four statuses. **FAIL**: nothing runs until it is fixed. **WARN**: some panels will not
build, or a number may differ from `assembled/`; the detail names which. **INFO**: a state
the scripts handle on their own (a fit not yet run, a gated input absent) reported so the
SKIP later is not a surprise. **PASS**.

On a public machine with every package installed, the expected result is zero FAIL, a WARN
for the Data package and each bucket, and INFO rows for the four fits and two gated inputs.
That machine builds 76 panels.

Every variable `config/landscape.env` declares is honoured, so the checks look where the
panels will: `EPIGEN_QC_DIR`, `PLIER_DIR`, `QUANT_BUCKET`, `GSUTIL` and the fit overrides.
`PREFLIGHT_REPORT=<path>` moves the report.

## 5. Docker

The same environment as an image, built from the same two manifests by the same installer,
so the container and §3 cannot disagree. R 4.4.1 on `rocker/r-ver`, Bioconductor 3.20, a
dated CRAN snapshot, the Google Cloud SDK, every panel package, and the Analysis package at
its current release. `docker-compose.yml` at the repository root is the entry point:

```bash
docker compose build landscape                                          # once
docker compose run --rm landscape bash figures/landscape/setup/preflight.sh
docker compose run --rm landscape Rscript figures/landscape/figure_2/FIG2.R
docker compose run --rm landscape Rscript figures/landscape/run_all.R
docker compose run --rm landscape                                       # a shell
```

Every run starts with a banner of what it is built from: R, Bioconductor, the three
package versions, and the branch and commit of the repository it is mounted over.

**What is mounted, and what is not.** The repository is bind-mounted read-write at
`/github/motrpac-human-presuspension-acute`, so `outputs/`, the epigenomics QC cache under
`staging/` and the fits land on the host and persist between runs. `~/.config/gcloud` is
mounted so `gcloud auth login` on the host serves the container. `/local-lib`, a named volume,
is the one writable R library. Nothing else on your machine is visible.

**The Data package in Docker.** It is not in the image. Either

- uncomment the `../MotrpacHumanPreSuspensionData` mount in `docker-compose.yml`: the
  entrypoint installs it into `/local-lib` on the first run and reinstalls only when the
  checkout's version or commit changes (`FORCE=1` reinstalls regardless), or
- bake it in at build time with a token that can read the private repository:

  ```bash
  GITHUB_PAT=ghp_... docker compose build landscape \
      --build-arg INSTALL_DATA=1 --secret id=github_pat,env=GITHUB_PAT
  ```

  The token is a BuildKit secret: it is available to that one build step and is not in any
  layer.

**Memory.** The PLIER fit holds ~8.5 GB; Docker Desktop's default VM limit is often lower.
Raise it in the Docker Desktop settings before `plier_fit.R`, or run that fit natively.
`VARIANCEPARTITION_PARALLEL_CORES=4 docker compose run ...` hands the two mixed-model fits
more cores.

**Build time and architecture.** On `linux/amd64` most CRAN packages arrive as binaries and
the build is tens of minutes; on Apple silicon (`linux/arm64`) they compile from source and
it is over an hour. `docker compose build landscape` builds for the host; add
`DOCKER_DEFAULT_PLATFORM=linux/amd64` for the binary path under emulation. The build context
is `figures/landscape/` with a whitelist `.dockerignore`, so `staging/` and `assembled/` are
never sent to the daemon.

**What the preflight shows inside the container.** With the default mounts: every tool and
package PASS, the Data package WARN, and the bucket rows PASS or WARN according to the
host's `gcloud` login. That is the scope boundary, not a broken image.

## 6. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `Library not loaded: /opt/X11/lib/libX11.6.dylib` installing the Analysis package | macOS without XQuartz. Install it (§2), log out and in, retry |
| `there is no package called 'TMSig'` or `'TMSig' is not available for Bioconductor version 3.18` | Bioconductor release older than 3.20. `BiocManager::install(version = "3.20")` on R 4.4, then re-run the installer |
| `package 'limma' was installed before R 4.4.0: please re-install it` | a library carried over from an older R. Reinstall the named package from the current release |
| `HTTP error 403` / `API rate limit exceeded` from `install_github()` | anonymous GitHub access is 60 requests/hour. Set `GITHUB_PAT` to any token and retry |
| `HTTP error 404` installing `MotrpacHumanPreSuspensionData` | the repository is private and the token cannot read it. The data-access request has to be approved and the account added |
| preflight: `pin:ggpubr ... blocking` | `ggpubr < 1.0.0`. `install.packages("ggpubr")`; the installer does this itself on a full run |
| preflight: `load:analysis_api` FAIL with `pin:MotrpacHumanPreSuspensionAnalysis` PASS | a 2.0.8 build from before `plot_enrich_heatmap(return_drawing = )` landed. Reinstall from GitHub |
| `[FAIL] FIG1C  ... CommandException` or `AccessDenied` | a bucket read without consortium access. §1 says which panels need it; everything else still builds |
| `cannot open URL 'https://d1yw74buhe0ts0.cloudfront.net/...'` | the public epigenomics DA download failed: network, or a proxy that blocks it. No credentials are involved |
| A panel reports `SKIPPED` naming a `*_fit.R` | the fit has not run. See the README's *The four fits* |
