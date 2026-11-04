# figures/landscape — the R environment the figure scripts run in, as an image.
#
# Built from the same manifest a native install reads (config/packages.tsv)
# by the same script (setup/install_packages.R), so the image and
# docs/INSTALL.md cannot disagree about
# what is installed. R 4.4.1, Bioconductor 3.20, a CRAN snapshot, the Google
# Cloud SDK, and every package the panels load, including
# MotrpacHumanPreSuspensionAnalysis at its current release.
#
# What this image is NOT: the data. The repository is bind-mounted at run time
# (see docker-compose.yml at the repository root), and the access-gated
# MotrpacHumanPreSuspensionData package is installed into a named volume from a
# mounted checkout by the entrypoint, or at build time when a GitHub token with
# access is passed as a secret. Consortium bucket access comes from the host's
# gcloud credentials, mounted in.
#
# Build (from the repository root):
#   docker compose build landscape
#   docker compose build landscape --build-arg INSTALL_DATA=1 \
#       --secret id=github_pat,env=GITHUB_PAT          # bake the Data package in
#
# rocker/r-ver is Ubuntu jammy (22.04); every apt name below is a jammy name.

ARG R_VERSION=4.4.1
FROM rocker/r-ver:${R_VERSION}

# The release that pairs with R 4.4 and the one the Analysis package's README
# names. TMSig, which that package imports, first appeared here.
ARG BIOC_VERSION=3.20
ENV BIOC_VERSION=${BIOC_VERSION}

# The CRAN snapshot every CRAN package comes from. Any date after ggpubr 1.0.0
# reached CRAN satisfies the one blocking pin; the base image's own default
# (2024-10-30) does not.
ARG CRAN_SNAPSHOT=2026-10-01

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8
ENV TZ=America/Chicago

# ---- System libraries --------------------------------------------------------
# What the package set compiles and links against. The non-obvious ones:
#   curl                  the binary, not just libcurl: rocker purges it, and
#                         setup/preflight.sh and gcloud's installer both call it.
#   libtcl8.6 libtk8.6    Mfuzz (imported by the Analysis package) imports tcltk;
#   xvfb                  loading it wants a display. The entrypoint starts Xvfb.
#   libglpk-dev           igraph, which annotatr pulls in through its Bioconductor
#                         neighbours.
#   libfontconfig/free-   ragg, textshaping and systemfonts, which ggrastr
#   type/harfbuzz/fribidi rasterises through.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates gnupg apt-transport-https \
      curl git make cmake perl python3 \
      libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
      libicu-dev libsqlite3-dev \
      zlib1g-dev libbz2-dev liblzma-dev \
      libcairo2-dev libx11-dev libxt-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libjpeg-dev libtiff-dev libwebp-dev \
      libglpk-dev libgmp-dev libfftw3-dev \
      libtcl8.6 libtk8.6 xvfb \
    && rm -rf /var/lib/apt/lists/*

# Headless display for tcltk; Xvfb is started by the entrypoint.
ENV DISPLAY=:99

# ---- Google Cloud SDK --------------------------------------------------------
# gsutil is what load_qc(epigen = TRUE) and the SF1 raw-count resolver shell
# out to. Published for amd64 and arm64; gsutil ships inside google-cloud-cli.
RUN curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
      | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
      > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# ---- CRAN snapshot -----------------------------------------------------------
# rocker bakes its default snapshot into Rprofile.site; append ours so it wins.
# P3M serves jammy binaries for amd64 only, so arm64 takes the source path.
RUN . /etc/os-release \
    && if [ "$(uname -m)" = "aarch64" ]; then \
         REPO="https://p3m.dev/cran/${CRAN_SNAPSHOT}"; \
       else \
         REPO="https://p3m.dev/cran/__linux__/${UBUNTU_CODENAME}/${CRAN_SNAPSHOT}"; \
       fi \
    && echo "options(repos = c(CRAN = '${REPO}'), download.file.method = 'libcurl')" \
       >> "${R_HOME}/etc/Rprofile.site" \
    && echo "CRAN snapshot: ${REPO}"

# ---- R packages, from the manifest -------------------------------------------
# Two layers so a change to the GitHub pins does not rebuild the ~40-package
# CRAN/Bioconductor layer beneath it. The cache mount keeps downloaded tarballs
# across retries.
COPY config/packages.tsv /opt/landscape/config/packages.tsv
COPY setup/install_packages.R /opt/landscape/setup/install_packages.R

RUN --mount=type=cache,target=/tmp/Rdownloads,sharing=locked \
    TMPDIR=/tmp/Rdownloads Rscript /opt/landscape/setup/install_packages.R \
      --root=/opt/landscape --stage=cran-bioc --bioc=${BIOC_VERSION}

# The secret is optional: without it the three public packages install
# anonymously, and the Data package is left to the entrypoint.
ARG INSTALL_DATA=0
RUN --mount=type=cache,target=/tmp/Rdownloads,sharing=locked \
    --mount=type=secret,id=github_pat \
    GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null || true)" \
    TMPDIR=/tmp/Rdownloads Rscript /opt/landscape/setup/install_packages.R \
      --root=/opt/landscape --stage=github --bioc=${BIOC_VERSION} \
      $([ "${INSTALL_DATA}" = "1" ] && echo --with-data)

# ---- Runtime -----------------------------------------------------------------
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Everything above is the image's read-only site library. /local-lib is the one
# writable library, backed by a named volume: the Data package installed from a
# mounted checkout lands there and survives between runs. Set last so nothing
# above installs into it.
RUN mkdir -p /local-lib
ENV R_LIBS_USER=/local-lib

ENV GITHUB_ROOT=/github
WORKDIR /github/motrpac-human-presuspension-acute

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
