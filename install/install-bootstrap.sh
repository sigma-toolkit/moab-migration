#!/usr/bin/env bash
#
# install-bootstrap.sh -- self-fetching wrapper for MOAB's install-moab.sh.
#
# Designed for the curl|bash handoff pattern, so a user can install MOAB
# (and its three TPLs eigen3/zoltan/tempestremap) with a single invocation
# and no manual file shuffling. The bootstrap fetches the rest of the
# install/ directory from a configurable URL into a local cache dir, then
# exec's install-moab.sh with whatever args were passed through.
#
# Usage (one-liner):
#
#     curl -fsSL https://bitbucket.org/fathomteam/moab/raw/master/install/install-bootstrap.sh \
#         | bash -s -- --machine=auto --dry-run
#
# Or save and run:
#
#     curl -fsSLo install-bootstrap.sh \
#         https://bitbucket.org/fathomteam/moab/raw/master/install/install-bootstrap.sh
#     bash install-bootstrap.sh --machine=auto --hdf5-root=$HDF5_ROOT \
#         --netcdf-root=$NETCDF_C_PATH --pnetcdf-root=$PNETCDF_PATH
#
# Env-var configuration (all optional):
#
#     INSTALL_MOAB_BRANCH         Branch / tag to fetch     (default: master)
#     INSTALL_MOAB_REPO_RAW_URL   Raw-URL base for the repo (default: bitbucket fathomteam/moab)
#     INSTALL_MOAB_DIR            Where to download install/ (default: $PWD/install)
#     INSTALL_MOAB_REFRESH        Force re-download (yes|no) (default: no)
#
# Examples:
#
#     # Fetch from a fork
#     INSTALL_MOAB_REPO_RAW_URL=https://bitbucket.org/myuser/moab/raw \
#         curl -fsSL .../install-bootstrap.sh | bash -s -- --machine=auto
#
#     # Pin to a release branch
#     INSTALL_MOAB_BRANCH=release-v5.6.0 \
#         curl -fsSL .../install-bootstrap.sh | bash -s -- --machine=auto
#
#     # Cache install/ in $HOME so re-runs skip the download
#     INSTALL_MOAB_DIR=$HOME/.local/share/moab-installer \
#         curl -fsSL .../install-bootstrap.sh | bash -s -- --machine=auto

set -euo pipefail

BRANCH="${INSTALL_MOAB_BRANCH:-master}"
REPO_RAW_URL="${INSTALL_MOAB_REPO_RAW_URL:-https://bitbucket.org/fathomteam/moab/raw}"
INSTALL_DIR="${INSTALL_MOAB_DIR:-$PWD/install}"
REFRESH="${INSTALL_MOAB_REFRESH:-no}"

URL_BASE="$REPO_RAW_URL/$BRANCH/install"

# Files that constitute the install/ directory. Keep in sync with the
# install_moab_files variable in the top-level Makefile.am.
FILES=(
    "install-moab.sh"
    "install-moab-e3sm.sh"
    "INSTALL-MOAB.md"
    "CONTRIBUTING-MACHINES.md"
    "workflow.sh"
    "scripts/e3sm_env.py"
)

# Subset that should end up executable.
EXEC_FILES=(
    "install-moab.sh"
    "install-moab-e3sm.sh"
    "workflow.sh"
    "scripts/e3sm_env.py"
)

log() { printf '[install-bootstrap] %s\n' "$*" >&2; }
die() { printf '[install-bootstrap] ERROR: %s\n' "$*" >&2; exit 1; }

# Pick a downloader. Prefer curl (more common on HPC); fall back to wget.
if command -v curl >/dev/null 2>&1; then
    fetch() { curl -fsSL --retry 3 --retry-delay 2 -o "$1" "$2"; }
elif command -v wget >/dev/null 2>&1; then
    fetch() { wget -q --tries=3 --waitretry=2 -O "$1" "$2"; }
else
    die "neither 'curl' nor 'wget' is on PATH; install one and retry"
fi

# Bash version sanity (install-moab.sh requires >= 3.2)
if (( BASH_VERSINFO[0] < 3 )) || { (( BASH_VERSINFO[0] == 3 )) && (( BASH_VERSINFO[1] < 2 )); }; then
    die "bash >= 3.2 required (you have $BASH_VERSION)"
fi

# If the cache dir already has install-moab.sh and we weren't asked to refresh,
# skip the download and exec the cached copy directly. Lets re-runs avoid the
# network round-trip.
if [[ -x "$INSTALL_DIR/install-moab.sh" && "$REFRESH" != "yes" ]]; then
    log "Reusing cached install scripts at $INSTALL_DIR"
    log "  (set INSTALL_MOAB_REFRESH=yes to re-download from $URL_BASE)"
else
    log "Fetching install scripts from $URL_BASE"
    log "  -> $INSTALL_DIR"
    mkdir -p "$INSTALL_DIR/scripts"
    for f in "${FILES[@]}"; do
        log "  $f"
        fetch "$INSTALL_DIR/$f" "$URL_BASE/$f" \
            || die "failed to fetch $f (check URL: $URL_BASE/$f)"
    done
    for f in "${EXEC_FILES[@]}"; do
        chmod +x "$INSTALL_DIR/$f"
    done
    log "Bootstrap complete."
fi

log "Exec: $INSTALL_DIR/install-moab.sh $*"
exec "$INSTALL_DIR/install-moab.sh" "$@"
