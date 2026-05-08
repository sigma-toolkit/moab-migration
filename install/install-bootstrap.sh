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
# Configuration: CLI flags (preferred for the curl|bash pattern -- args
# survive the pipe naturally) OR env vars. Flags take precedence over env.
#
#     --bootstrap-branch=NAME     Branch / tag to fetch       (default: master)
#     --bootstrap-repo=URL        Raw-URL base for the repo   (default: bitbucket fathomteam/moab)
#     --bootstrap-dir=PATH        Where to download install/  (default: $PWD/install)
#     --bootstrap-refresh         Force re-download even if cache present
#
#     INSTALL_MOAB_BRANCH         (same; env-var fallback)
#     INSTALL_MOAB_REPO_RAW_URL   (same)
#     INSTALL_MOAB_DIR            (same)
#     INSTALL_MOAB_REFRESH=yes    (same)
#
# Important: env vars set BEFORE 'curl' do NOT propagate to the piped bash.
# Use --bootstrap-* flags for the curl|bash pattern, e.g.:
#
#     curl -fsSL https://bitbucket.org/fathomteam/moab/raw/myfork-branch/install/install-bootstrap.sh \
#         | bash -s -- --bootstrap-branch=myfork-branch --machine=auto
#
# (The bootstrap intercepts its own --bootstrap-* flags and forwards the
# rest to install-moab.sh.)
#
# Examples:
#
#     # Fetch from a fork's branch (one-liner, args survive the pipe)
#     curl -fsSL https://bitbucket.org/myuser/moab/raw/master/install/install-bootstrap.sh \
#         | bash -s -- --bootstrap-repo=https://bitbucket.org/myuser/moab/raw \
#                      --machine=auto
#
#     # Pin to a release branch
#     curl -fsSL .../install-bootstrap.sh \
#         | bash -s -- --bootstrap-branch=release-v5.6.0 --machine=auto
#
#     # Cache install/ in $HOME so re-runs across build dirs share one copy
#     curl -fsSL .../install-bootstrap.sh \
#         | bash -s -- --bootstrap-dir=$HOME/.local/share/moab-installer \
#                      --machine=auto

set -euo pipefail

log() { printf '[install-bootstrap] %s\n' "$*" >&2; }
die() { printf '[install-bootstrap] ERROR: %s\n' "$*" >&2; exit 1; }

# Defaults (env vars seed the value; CLI flags override below)
BRANCH="${INSTALL_MOAB_BRANCH:-master}"
REPO_RAW_URL="${INSTALL_MOAB_REPO_RAW_URL:-https://bitbucket.org/fathomteam/moab/raw}"
INSTALL_DIR="${INSTALL_MOAB_DIR:-$PWD/install}"
REFRESH="${INSTALL_MOAB_REFRESH:-no}"

# Intercept --bootstrap-* flags; pass the rest to install-moab.sh via exec.
FORWARDED_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bootstrap-branch=*)   BRANCH="${1#*=}" ;;
        --bootstrap-repo=*)     REPO_RAW_URL="${1#*=}" ;;
        --bootstrap-dir=*)      INSTALL_DIR="${1#*=}" ;;
        --bootstrap-refresh)    REFRESH=yes ;;
        --bootstrap-help)
            sed -n '1,/^set -euo pipefail$/p' "${BASH_SOURCE[0]:-$0}" \
                | sed 's/^# \?//; /^set -euo/d; /^#!/d'
            exit 0
            ;;
        *) FORWARDED_ARGS+=("$1") ;;
    esac
    shift
done

URL_BASE="$REPO_RAW_URL/$BRANCH/install"

# Print resolved config up-front so the user can spot "wrong branch" /
# "wrong repo" mistakes before the fetch fails. Particularly useful when
# the curl|bash pattern silently drops env vars (set before 'curl', not
# 'bash') -- you'll see BRANCH=master here when you expected something
# else, making it obvious to add --bootstrap-branch=NAME.
log "Bootstrap config:"
log "  branch     : $BRANCH"
log "  repo (raw) : $REPO_RAW_URL"
log "  cache dir  : $INSTALL_DIR"
log "  refresh    : $REFRESH"
log ""

# Files that constitute the install/ directory. Keep in sync with the
# install_moab_files variable in the top-level Makefile.am.
FILES=(
    "install-moab.sh"
    "INSTALL-MOAB.md"
    "CONTRIBUTING-MACHINES.md"
    "scripts/e3sm_env.py"
)

# Subset that should end up executable.
EXEC_FILES=(
    "install-moab.sh"
    "scripts/e3sm_env.py"
)

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

# Decide whether to reuse the cache. Reuse skips the download AND avoids
# silently running stale scripts that were fetched from a different branch
# / repo than the user is now requesting. The .bootstrap-source stamp is
# written after each successful fetch and read here for source-of-truth
# comparison.
STAMP_FILE="$INSTALL_DIR/.bootstrap-source"
SHOULD_REUSE="no"
if [[ -x "$INSTALL_DIR/install-moab.sh" && "$REFRESH" != "yes" ]]; then
    if [[ -f "$STAMP_FILE" ]]; then
        cached_branch="$(awk -F= '/^BRANCH=/{print $2; exit}' "$STAMP_FILE" 2>/dev/null || true)"
        cached_repo="$(awk -F= '/^REPO=/{print $2; exit}' "$STAMP_FILE" 2>/dev/null || true)"
        if [[ "$cached_branch" == "$BRANCH" && "$cached_repo" == "$REPO_RAW_URL" ]]; then
            SHOULD_REUSE="yes"
        else
            log "Cache mismatch detected:"
            log "  cached:    branch=$cached_branch repo=$cached_repo"
            log "  requested: branch=$BRANCH repo=$REPO_RAW_URL"
            log "  -> re-downloading"
        fi
    else
        # No stamp file = fetched by an older bootstrap version, OR the
        # directory was populated by hand. Treat as untrusted and re-download
        # once; the new fetch writes a stamp for next time.
        log "Cache present at $INSTALL_DIR but has no source stamp"
        log "  (probably from an older bootstrap version)"
        log "  -> re-downloading once to write the stamp"
    fi
fi

if [[ "$SHOULD_REUSE" == "yes" ]]; then
    log "Reusing cached install scripts at $INSTALL_DIR"
    log "  (pass --bootstrap-refresh to force re-download from $URL_BASE)"
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
    # Stamp the cache with the source so future invocations can detect
    # mismatch (different branch, different fork, etc.).
    {
        echo "BRANCH=$BRANCH"
        echo "REPO=$REPO_RAW_URL"
        echo "FETCHED_AT=$(date -u +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || date -u)"
    } > "$STAMP_FILE"
    log "Bootstrap complete."
fi

log "Exec: $INSTALL_DIR/install-moab.sh ${FORWARDED_ARGS[*]:-}"
exec "$INSTALL_DIR/install-moab.sh" ${FORWARDED_ARGS[@]+"${FORWARDED_ARGS[@]}"}
