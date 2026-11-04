#!/usr/bin/env bash
# Container entrypoint: say what this run is built from, install the gated Data
# package from a mounted checkout if one is there, then hand off.
#
# The banner exists because the image pins everything except the repository it
# runs over and the Data checkout it may be handed, and those two are what
# decide whether a panel matches assembled/. So it says out loud what it got.
set -euo pipefail

GITHUB_ROOT="${GITHUB_ROOT:-/github}"
LOCAL_LIB="${R_LIBS_USER:-/local-lib}"

git_state() {
  local dir="$1"
  # -e, not -d: in a linked worktree .git is a file pointing at the real git dir.
  [[ -e "${dir}/.git" ]] || { printf 'not a git checkout'; return 0; }
  local branch sha dirty
  branch="$(git -C "${dir}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
  sha="$(git -C "${dir}" rev-parse --short HEAD 2>/dev/null || echo '?')"
  dirty="$(git -C "${dir}" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  if [[ "${dirty}" == "0" ]]; then printf '%s @ %s (clean)' "${branch}" "${sha}"
  else printf '%s @ %s (%s uncommitted)' "${branch}" "${sha}" "${dirty}"; fi
}

# git refuses to read a repository owned by another uid, which a bind mount
# from the host is.
git config --global --add safe.directory '*' 2>/dev/null || true

# Mfuzz imports tcltk, and attaching the Analysis package loads Mfuzz. tcltk
# wants a display; give it a headless one if DISPLAY names a local screen that
# is not already up.
if [[ -n "${DISPLAY:-}" && "${DISPLAY}" == :* ]] && command -v Xvfb >/dev/null 2>&1; then
  screen="${DISPLAY#:}"; screen="${screen%%.*}"
  if [[ ! -e "/tmp/.X11-unix/X${screen}" ]]; then
    Xvfb "${DISPLAY}" -screen 0 1280x1024x24 >/dev/null 2>&1 &
    for _ in 1 2 3 4 5; do
      [[ -e "/tmp/.X11-unix/X${screen}" ]] && break
      sleep 0.3
    done
  fi
fi

pkg_version() {
  Rscript -e "cat(tryCatch(as.character(packageVersion('$1')), error = function(e) 'not installed'))" 2>/dev/null || echo '?'
}

# ---- The Data package, from a mounted checkout -------------------------------
# docker-compose.yml has the mount commented out; when it is enabled, install
# from it into the writable library, and reinstall only when the checkout's
# version, commit or dirty state changes. FORCE=1 reinstalls regardless.
DATA_REPO="${GITHUB_ROOT}/MotrpacHumanPreSuspensionData"
if [[ "${SKIP_LOCAL_PKG_INSTALL:-0}" != "1" && -f "${DATA_REPO}/DESCRIPTION" ]]; then
  version="$(awk '/^Version:/ {print $2; exit}' "${DATA_REPO}/DESCRIPTION")"
  if [[ -e "${DATA_REPO}/.git" ]]; then
    want="${version} $(git -C "${DATA_REPO}" rev-parse --short HEAD 2>/dev/null || echo '?') dirty=$(git -C "${DATA_REPO}" status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
  else
    want="${version} no-git"
  fi
  stamp="${LOCAL_LIB}/.installed-MotrpacHumanPreSuspensionData"
  have="$(cat "${stamp}" 2>/dev/null || true)"
  if [[ "${FORCE:-0}" != "1" && "${want}" == "${have}" && -d "${LOCAL_LIB}/MotrpacHumanPreSuspensionData" ]]; then
    echo "[data-pkg] MotrpacHumanPreSuspensionData up to date (${want})"
  else
    echo "[data-pkg] installing MotrpacHumanPreSuspensionData (${want}) -> ${LOCAL_LIB}"
    mkdir -p "${LOCAL_LIB}"
    if R CMD INSTALL --no-docs --library="${LOCAL_LIB}" "${DATA_REPO}"; then
      printf '%s' "${want}" > "${stamp}"
    else
      rm -f "${stamp}"
      echo "[data-pkg] FAILED to install from ${DATA_REPO}" >&2
      exit 1
    fi
  fi
fi

echo "─── motrpac-human-presuspension-acute · figures/landscape ─────────────────────"
printf '  %-36s %s\n' "R" "$(R --version 2>/dev/null | head -n1 | sed 's/ --.*//')"
printf '  %-36s %s\n' "Bioconductor" "$(Rscript -e 'cat(as.character(BiocManager::version()))' 2>/dev/null || echo '?')"
printf '  %-36s %s\n' "MotrpacHumanPreSuspensionAnalysis" "$(pkg_version MotrpacHumanPreSuspensionAnalysis)"
printf '  %-36s %s\n' "MotrpacHumanPreSuspensionData" "$(pkg_version MotrpacHumanPreSuspensionData)"
printf '  %-36s %s\n' "MotrpacBicQC" "$(pkg_version MotrpacBicQC)"
printf '  %-36s %s\n' "repository" "$(git_state "${GITHUB_ROOT}/motrpac-human-presuspension-acute")"
[[ -d "${DATA_REPO}" ]] && printf '  %-36s %s\n' "Data checkout" "$(git_state "${DATA_REPO}")"
echo "──────────────────────────────────────────────────────────────────────────────"

exec "$@"
