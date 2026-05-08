#!/usr/bin/env bash
#
# install-moab-e3sm.sh
#
# Robust orchestrator that builds MOAB for use with E3SM. The orchestrator
# owns the full configure/build/install lifecycle for the three TPLs that
# E3SM-MOAB needs (eigen3, zoltan, tempestremap) and hands the resulting
# install paths to MOAB. MOAB itself is therefore built generically against
# pre-built TPLs, with no dependency on its --download-* macros.
#
# This decoupling means:
#   * MOAB can be driven via either autotools (--build-system=autotools,
#     default) or CMake (--build-system=cmake) using exactly the same TPL
#     stack
#   * each TPL has its own download/configure/build/install retry loop with
#     targeted remediations (network, GCC>=10 -fcommon, missing autoreconf,
#     NetCDF/HDF5 hint vars, C++14, LAPACK)
#   * a transient failure in one TPL never poisons MOAB's configure
#
# Required environment from the E3SM env (or pass via flags):
#   MPI_ROOT, HDF5_ROOT, NETCDF_C_PATH, PNETCDF_PATH
#
# Run with --help for full usage.

set -euo pipefail

# Distinguish env-set hints from explicit-flag bindings.
#   _ENV_HAS_*  -- variable came in from the user's environment (just a hint)
#   _FLAG_HAS_* -- user explicitly passed --cc/--cxx/--fc/--f77 on the CLI
# Env hints are overridden by MPI wrappers from MPI_ROOT or PATH; only the
# flag form is treated as a hard, MPI-wrapper-overriding override.
_ENV_HAS_CC=${CC+yes}
_ENV_HAS_CXX=${CXX+yes}
_ENV_HAS_FC=${FC+yes}
_ENV_HAS_F77=${F77+yes}
_FLAG_HAS_CC=no
_FLAG_HAS_CXX=no
_FLAG_HAS_FC=no
_FLAG_HAS_F77=no
_USER_SET_MOAB_SRC=${MOAB_SRC_DIR+yes}
_USER_SET_BUILD_DIR=${BUILD_DIR+yes}

#-----------------------------------------------------------------------------
# Defaults (override via env or CLI)
#-----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PREFIX_PATH="${PREFIX_PATH:-$HOME/install/MOAB}"
TPL_PREFIX="${TPL_PREFIX:-$PREFIX_PATH/tpls}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

# MOAB source management. The script clones MOAB if absent and fast-forwards
# to the latest $MOAB_BRANCH HEAD on every run, so a re-invocation also picks
# up upstream changes. To opt out (use a private dev checkout), pass
# --src-dir=PATH (no git operations are performed against that path).
MOAB_REPO_URL="${MOAB_REPO_URL:-https://bitbucket.org/fathomteam/moab.git}"
MOAB_BRANCH="${MOAB_BRANCH:-master}"
NO_SOURCE_UPDATE="${NO_SOURCE_UPDATE:-no}"

BUILD_SYSTEM="${BUILD_SYSTEM:-autotools}"   # autotools | cmake
ENABLE_DEBUG="${ENABLE_DEBUG:-yes}"
ENABLE_OPTIMIZE="${ENABLE_OPTIMIZE:-yes}"
ENABLE_SHARED="${ENABLE_SHARED:-no}"
ENABLE_STATIC="${ENABLE_STATIC:-yes}"
RUN_CHECK="${RUN_CHECK:-no}"
CLEAN_BUILD="${CLEAN_BUILD:-no}"
CLEAN_TPLS="${CLEAN_TPLS:-no}"
DRY_RUN="${DRY_RUN:-no}"
PRINT_MODE="${PRINT_MODE:-no}"
RECONFIGURE="${RECONFIGURE:-no}"
SKIP_MPI_VALIDATION="${SKIP_MPI_VALIDATION:-no}"

# Skip an individual TPL build if it is already installed at $TPL_PREFIX/<tpl>
# (set to "no" to always rebuild)
REUSE_TPLS="${REUSE_TPLS:-yes}"

# Robustness knobs
DOWNLOAD_RETRIES="${DOWNLOAD_RETRIES:-5}"
DOWNLOAD_RETRY_DELAY="${DOWNLOAD_RETRY_DELAY:-10}"
TPL_MAX_ATTEMPTS="${TPL_MAX_ATTEMPTS:-3}"
TAIL_LOGS="${TAIL_LOGS:-yes}"
TPL_MIRROR="${TPL_MIRROR:-https://web.cels.anl.gov/projects/sigma/downloads/TPL}"

# Compiler wrappers
CC_BIN="${CC:-mpicc}"
CXX_BIN="${CXX:-mpicxx}"
FC_BIN="${FC:-mpif90}"
F77_BIN="${F77:-mpif77}"

EXTRA_MOAB_ARGS="${EXTRA_MOAB_ARGS:-}"
EXTRA_ZOLTAN_ARGS="${EXTRA_ZOLTAN_ARGS:-}"
EXTRA_TEMPESTREMAP_ARGS="${EXTRA_TEMPESTREMAP_ARGS:-}"

# Machine selection (resolved later, in apply_machine_defaults)
MACHINE_NAME="${MACHINE_NAME:-auto}"
COMPILER_FAMILY="${COMPILER_FAMILY:-}"

# Profile + E3SM checkout
PROFILE="${PROFILE:-e3sm}"           # e3sm | standalone
E3SM_ROOT="${E3SM_ROOT:-}"           # Required for --profile=e3sm
ASSUME_YES="${ASSUME_YES:-no}"       # Bypass the e3sm-profile env confirmation prompt

# TPL versions (must match MOAB's expected sources for consistency with E3SM)
EIGEN3_VERSION="3.4.0"
ZOLTAN_VERSION="3.9.1"
TEMPESTREMAP_VERSION="2.2.x"

#-----------------------------------------------------------------------------
# Machine database
#-----------------------------------------------------------------------------
# Per-machine entries are pure metadata (no hardcoded TPL paths) so the script
# stays maintainable across path migrations and avoids leaking site-specific
# layouts. Push 2 will use the metadata to drive E3SM CIME-based environment
# resolution (--profile=e3sm); for now (Push 1), the metadata is informational
# only -- the script reports the detected machine in the orchestration banner
# and prints a module-load hint, but does not auto-load anything.
#
# Adding a new machine: define machine_<name>_match() and machine_<name>_meta(),
# then append <name> to MACHINE_REGISTRY. Mark new entries with
# last_validated=TBD until you have personally confirmed a successful build.
#
# meta() sets these globals:
#   MACHINE_META_E3SM_NAME         -- name in E3SM's config_machines.xml (may differ)
#   MACHINE_META_DEFAULT_COMPILER  -- compiler family used if --compiler is unset
#   MACHINE_META_SUPPORTED_COMPILERS -- comma-separated list (informational)
#   MACHINE_META_STANDALONE_HINT   -- module-load string printed for standalone users
#   MACHINE_META_LAST_VALIDATED    -- YYYY-MM-DD or TBD
#   MACHINE_META_NOTES             -- optional free-form (Cray PrgEnv quirks, etc.)

MACHINE_REGISTRY="bebop improv crux gce perlmutter"

# Default-detection hostname helpers
_hn() { printf '%s' "${HOSTNAME:-$(hostname 2>/dev/null || true)}"; }

#---------- Bebop (LCRC, ANL) ----------
machine_bebop_match() {
    [[ "${LMOD_SYSTEM_NAME:-}" == "bebop" ]] && return 0
    [[ "$(_hn)" == bebop* ]]
}
machine_bebop_meta() {
    MACHINE_META_E3SM_NAME="bebop"          # matches MACH= in cime_config/machines/config_machines.xml
    MACHINE_META_DEFAULT_COMPILER="gnu"
    MACHINE_META_SUPPORTED_COMPILERS="gnu,intel"
    MACHINE_META_STANDALONE_HINT="module load gcc/13.2.0 openmpi/4.1.8 hdf5/1.12.3 netcdf-c parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="2026-05-08"
    MACHINE_META_NOTES="config_machines.xml does NOT export HDF5_ROOT; pass --hdf5-root=PATH or load a hdf5 module that sets it. BLAS/LAPACK come from \$LAPACK_ROOT/\$BLAS_ROOT in the e3sm profile."
}

#---------- Improv (LCRC, ANL) ----------
machine_improv_match() {
    [[ "${LMOD_SYSTEM_NAME:-}" == "improv" ]] && return 0
    [[ "$(_hn)" == improv* ]]
}
machine_improv_meta() {
    MACHINE_META_E3SM_NAME="improv"         # TODO verify
    MACHINE_META_DEFAULT_COMPILER="gnu"
    MACHINE_META_SUPPORTED_COMPILERS="gnu,intel,aocc"
    MACHINE_META_STANDALONE_HINT="module load gcc/12.3.0 openmpi hdf5 netcdf-c parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="TBD"
    MACHINE_META_NOTES="Newer LCRC cluster (replaces Bebop). Validate before promoting last_validated."
}

#---------- Crux (ALCF, ANL) ----------
machine_crux_match() {
    [[ "${LMOD_SYSTEM_NAME:-}" == "crux" ]] && return 0
    [[ "$(_hn)" == crux* ]]
}
machine_crux_meta() {
    MACHINE_META_E3SM_NAME="crux"           # TODO verify
    MACHINE_META_DEFAULT_COMPILER="gnu"
    MACHINE_META_SUPPORTED_COMPILERS="gnu,cray,nvhpc"
    MACHINE_META_STANDALONE_HINT="module load PrgEnv-gnu cray-hdf5-parallel cray-netcdf-hdf5parallel cray-parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="TBD"
    MACHINE_META_NOTES="Cray PrgEnv. Compiler wrappers are cc/CC/ftn (auto-detected); --mpi-root not needed."
}

#---------- ANL/GCE (Linux Ubuntu) ----------
machine_gce_match() {
    [[ "${LMOD_SYSTEM_NAME:-}" == "gce" ]] && return 0
    case "$(_hn)" in
        gce*|*.gce.anl.gov|gce-*) return 0 ;;
    esac
    return 1
}
machine_gce_meta() {
    MACHINE_META_E3SM_NAME="anlgce-ub22"      # ANL/GCE Ubuntu 22 entry in config_machines.xml
    MACHINE_META_DEFAULT_COMPILER="gnu"
    MACHINE_META_SUPPORTED_COMPILERS="gnu,intel"
    MACHINE_META_STANDALONE_HINT="module load gcc mpich hdf5 netcdf-c parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="2026-05-08"
    MACHINE_META_NOTES="ANL Climate dev cluster. Spack-managed modules; ensure HDF5/NetCDF/PNetCDF are loaded before invoking."
}

#---------- Perlmutter (NERSC) ----------
machine_perlmutter_match() {
    [[ "${NERSC_HOST:-}" == "perlmutter" ]] && return 0
    [[ "${LMOD_SYSTEM_NAME:-}" == "perlmutter" ]] && return 0
    [[ "$(_hn)" == nid* || "$(_hn)" == login* ]] && [[ -d /opt/cray/pe ]] && return 0
    return 1
}
machine_perlmutter_meta() {
    MACHINE_META_E3SM_NAME="pm-cpu"          # E3SM uses pm-cpu / pm-gpu (override with --machine and --e3sm-name=pm-gpu manually for now)
    MACHINE_META_DEFAULT_COMPILER="gnu"
    MACHINE_META_SUPPORTED_COMPILERS="gnu,intel,nvidia,aocc"
    MACHINE_META_STANDALONE_HINT="module load PrgEnv-gnu cray-hdf5-parallel cray-netcdf-hdf5parallel cray-parallel-netcdf"
    MACHINE_META_LAST_VALIDATED="2026-05-08"
    MACHINE_META_NOTES="Cray PrgEnv. Compiler wrappers cc/CC/ftn auto-detected. Pass --machine=perlmutter --compiler=nvidia for GPU builds."
}

# Reset the meta globals before invoking a machine_<name>_meta() function
reset_machine_meta() {
    MACHINE_META_E3SM_NAME=""
    MACHINE_META_DEFAULT_COMPILER=""
    MACHINE_META_SUPPORTED_COMPILERS=""
    MACHINE_META_STANDALONE_HINT=""
    MACHINE_META_LAST_VALIDATED=""
    MACHINE_META_NOTES=""
}

# Returns the name of the first registered machine whose match() returns 0,
# or empty if none match.
detect_machine() {
    local m
    for m in $MACHINE_REGISTRY; do
        if "machine_${m}_match" 2>/dev/null; then
            printf '%s' "$m"
            return 0
        fi
    done
    printf ''
}

# Validates that $1 is a registered machine name. Returns 0 if so.
is_known_machine() {
    local m="$1" r
    for r in $MACHINE_REGISTRY; do
        [[ "$r" == "$m" ]] && return 0
    done
    return 1
}

#-----------------------------------------------------------------------------
# Pretty output
#-----------------------------------------------------------------------------
COLOR_RED=$'\033[0;31m'
COLOR_GREEN=$'\033[0;32m'
COLOR_YELLOW=$'\033[0;33m'
COLOR_BLUE=$'\033[0;34m'
COLOR_RESET=$'\033[0m'
[[ -t 1 ]] || { COLOR_RED=""; COLOR_GREEN=""; COLOR_YELLOW=""; COLOR_BLUE=""; COLOR_RESET=""; }

# All progress output goes to stderr so functions can return values via stdout
# (e.g. fetch_tpl_archive prints the archive path; capturing with $(...) must
# not pick up log lines).
log()    { printf '%s[install-moab]%s %s\n' "$COLOR_BLUE" "$COLOR_RESET" "$*" >&2; }
ok()     { printf '%s[install-moab]%s %s%s%s\n' "$COLOR_BLUE" "$COLOR_RESET" "$COLOR_GREEN" "$*" "$COLOR_RESET" >&2; }
warn()   { printf '%s[install-moab]%s %sWARNING:%s %s\n' "$COLOR_BLUE" "$COLOR_RESET" "$COLOR_YELLOW" "$COLOR_RESET" "$*" >&2; }
die()    { printf '%s[install-moab]%s %sERROR:%s %s\n' "$COLOR_BLUE" "$COLOR_RESET" "$COLOR_RED" "$COLOR_RESET" "$*" >&2; exit "${2:-1}"; }
section(){ printf '\n%s========== %s ==========%s\n\n' "$COLOR_BLUE" "$*" "$COLOR_RESET" >&2; }

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Builds MOAB and its TPL stack (eigen3, zoltan, tempestremap) for E3SM.
Each TPL is built standalone and then handed to MOAB as a pre-built
dependency, so MOAB can be driven via either autotools or CMake.

Required environment variables (or flags):
  HDF5_ROOT         HDF5 installation root       (--with-hdf5)
  NETCDF_C_PATH     NetCDF-C installation root   (--with-netcdf)
  PNETCDF_PATH      Parallel-NetCDF root         (--with-pnetcdf)

Optional but recommended:
  MPI_ROOT          MPI installation root. If set, compiler wrappers are taken
                    from \$MPI_ROOT/bin/mpicc etc., and --with-mpi=\$MPI_ROOT is
                    passed to MOAB. If unset (e.g. Cray PrgEnv where wrappers
                    are cc/CC/ftn), the script falls back to PATH lookup of
                    mpicc/mpicxx/mpif90/mpif77 with sensible alternates
                    (mpic++, mpifort, mpiifort, ftn).
                    User-set CC/CXX/FC/F77 (env or --cc etc.) always win.

Common options:
  --prefix=PATH         MOAB install prefix          [PREFIX_PATH=$PREFIX_PATH]
  --tpl-prefix=PATH     TPL install root             [TPL_PREFIX=\$PREFIX_PATH/tpls]
  --build-dir=PATH      Build/work dir               [BUILD_DIR=\$PWD/moab-build]
  --src-dir=PATH        Use this MOAB checkout (NO git operations performed).
                        If omitted, the script clones into \$BUILD_DIR/moab-src
                        and updates it on every run.
  --moab-repo=URL       MOAB git remote              [MOAB_REPO_URL=$MOAB_REPO_URL]
  --moab-branch=NAME    MOAB branch                  [MOAB_BRANCH=$MOAB_BRANCH]
  --no-source-update    Skip git fetch+reset on the script-managed clone
  --jobs=N              Parallel make jobs           [JOBS=$JOBS]
  --build-system=BS     autotools | cmake            [BUILD_SYSTEM=$BUILD_SYSTEM]
  --mpi-root=PATH       MPI root
  --hdf5-root=PATH      HDF5 root
  --netcdf-root=PATH    NetCDF-C root
  --pnetcdf-root=PATH   Parallel-NetCDF root
  --cc=BIN              C compiler                   [CC=$CC_BIN]
  --cxx=BIN             C++ compiler                 [CXX=$CXX_BIN]
  --fc=BIN              Fortran compiler             [FC=$FC_BIN]
  --f77=BIN             F77 compiler                 [F77=$F77_BIN]
  --shared              Build shared libraries
  --no-static           Disable static (implies --shared)
  --no-debug            Disable debug symbols
  --no-optimize         Disable compiler optimization
  --check               Run "make check" after build
  --clean               Wipe BUILD_DIR before configuring (full restart)
  --clean-tpls          Wipe TPL_PREFIX before TPL builds
  --no-reuse-tpls       Always rebuild TPLs even if already installed
  --reconfigure         Force MOAB re-configure even if already configured
  --skip-mpi-validation Skip the MPI compile/link sanity test (only if you
                        cannot compile on this host)
  --no-tail             Suppress live tail of TPL build logs
  --dry-run             Print all commands, do nothing
  --print               Print only the resolved MOAB configure/cmake command
                        (copy-pasteable, no banner). Implies --dry-run.
  --profile=NAME        e3sm (default) | standalone. With e3sm, the script
                        sources modules + env vars from
                        \$E3SM_ROOT/cime_config/machines/config_machines.xml
                        for the resolved (machine, compiler) pair. With
                        standalone, the script trusts whatever environment
                        you have already loaded (suitable for MOAB downstream
                        users who don't have an E3SM checkout).
  --e3sm-root=PATH      Path to E3SM checkout (or set \$E3SM_ROOT). Required
                        for --profile=e3sm. Must contain
                        cime_config/machines/config_machines.xml.
  --yes, -y             Skip the e3sm-profile env confirmation prompt. Required
                        when stdin is not a tty (CI, piped input). Equivalent
                        to ASSUME_YES=yes.
  --machine=NAME        Use the named machine entry from the database. NAME=auto
                        (default) auto-detects via hostname/NERSC_HOST/
                        LMOD_SYSTEM_NAME. Pass --list-machines to see all entries.
                        Machine entries currently provide informational hints only;
                        Push 2 will use them to drive E3SM CIME env resolution.
  --compiler=NAME       Compiler family on the chosen machine (gnu, intel, cray,
                        nvidia, aocc, ...). Defaults to the machine entry's
                        default_compiler. Informational in Push 1.
  --list-machines       Print the registered machine entries with last-validated
                        dates, mark the auto-detected one, and exit.
  --extra=ARGS          Extra args appended to MOAB invocation
  --extra-zoltan=ARGS   Extra args appended to Zoltan's configure (verbatim)
  --extra-tempestremap=ARGS
                        Extra args appended to TempestRemap's configure (verbatim).
                        Use this for site-specific BLAS/LAPACK selection, e.g.
                          --extra-tempestremap="--with-blas=-L/opt/mkl/lib -lmkl_rt \
                                                --with-lapack=-L/opt/mkl/lib"
                        or pass LIBS=/LDFLAGS=/CPPFLAGS= overrides directly.
                        Changing this between runs auto-triggers a rebuild of that
                        TPL even if the install marker exists.

Resume behavior (default; opt out with --clean / --no-reuse-tpls / --reconfigure):
  * TPLs whose install marker exists are skipped
  * MOAB is not re-configured if its build dir is already configured AND the
    fingerprint of the configure args matches the last successful run
  * make / cmake --build / make install run incrementally on every invocation,
    so a transient compiler error is fixed by editing the source and re-running

Robustness knobs (env vars):
  DOWNLOAD_RETRIES=$DOWNLOAD_RETRIES         Tarball download attempts
  DOWNLOAD_RETRY_DELAY=$DOWNLOAD_RETRY_DELAY  Backoff seconds per retry
  TPL_MAX_ATTEMPTS=$TPL_MAX_ATTEMPTS              Per-TPL configure retry budget
  TPL_MIRROR=URL                              Override default TPL mirror
  EXTRA_MOAB_ARGS="..."                       Extra args appended to MOAB invocation
  EXTRA_ZOLTAN_ARGS="..."                     Extra args appended to Zoltan configure
  EXTRA_TEMPESTREMAP_ARGS="..."               Extra args appended to TempestRemap configure

  -h, --help            Show this help
EOF
}

#-----------------------------------------------------------------------------
# Argument parsing
#-----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix=*)        PREFIX_PATH="${1#*=}" ;;
        --tpl-prefix=*)    TPL_PREFIX="${1#*=}" ;;
        --build-dir=*)     BUILD_DIR="${1#*=}";  _USER_SET_BUILD_DIR=yes ;;
        --src-dir=*)       MOAB_SRC_DIR="${1#*=}"; _USER_SET_MOAB_SRC=yes ;;
        --moab-repo=*)     MOAB_REPO_URL="${1#*=}" ;;
        --moab-branch=*)   MOAB_BRANCH="${1#*=}" ;;
        --no-source-update) NO_SOURCE_UPDATE=yes ;;
        --jobs=*)          JOBS="${1#*=}" ;;
        --build-system=*)  BUILD_SYSTEM="${1#*=}" ;;
        --mpi-root=*)      MPI_ROOT="${1#*=}" ;;
        --hdf5-root=*)     HDF5_ROOT="${1#*=}" ;;
        --netcdf-root=*)   NETCDF_C_PATH="${1#*=}" ;;
        --pnetcdf-root=*)  PNETCDF_PATH="${1#*=}" ;;
        --cc=*)            CC_BIN="${1#*=}";  _FLAG_HAS_CC=yes ;;
        --cxx=*)           CXX_BIN="${1#*=}"; _FLAG_HAS_CXX=yes ;;
        --fc=*)            FC_BIN="${1#*=}";  _FLAG_HAS_FC=yes ;;
        --f77=*)           F77_BIN="${1#*=}"; _FLAG_HAS_F77=yes ;;
        --shared)          ENABLE_SHARED=yes ;;
        --no-static)       ENABLE_STATIC=no; ENABLE_SHARED=yes ;;
        --no-debug)        ENABLE_DEBUG=no ;;
        --no-optimize)     ENABLE_OPTIMIZE=no ;;
        --check)           RUN_CHECK=yes ;;
        --clean)           CLEAN_BUILD=yes ;;
        --clean-tpls)      CLEAN_TPLS=yes ;;
        --no-reuse-tpls)   REUSE_TPLS=no ;;
        --reconfigure)     RECONFIGURE=yes ;;
        --skip-mpi-validation) SKIP_MPI_VALIDATION=yes ;;
        --no-tail)         TAIL_LOGS=no ;;
        --dry-run)         DRY_RUN=yes ;;
        --print)           PRINT_MODE=yes; DRY_RUN=yes ;;
        --extra=*)               EXTRA_MOAB_ARGS="${1#*=}" ;;
        --extra-zoltan=*)        EXTRA_ZOLTAN_ARGS="${1#*=}" ;;
        --extra-tempestremap=*)  EXTRA_TEMPESTREMAP_ARGS="${1#*=}" ;;
        --machine=*)             MACHINE_NAME="${1#*=}" ;;
        --compiler=*)            COMPILER_FAMILY="${1#*=}" ;;
        --list-machines)         _LIST_MACHINES=yes ;;
        --profile=*)             PROFILE="${1#*=}" ;;
        --e3sm-root=*)           E3SM_ROOT="${1#*=}" ;;
        --yes|-y|--no-confirm)   ASSUME_YES=yes ;;
        -h|--help)         usage; exit 0 ;;
        *) die "unknown option: $1 (use --help)" 1 ;;
    esac
    shift
done

#-----------------------------------------------------------------------------
# --list-machines: print the registry and exit (does not require any env)
#-----------------------------------------------------------------------------
list_machines() {
    printf 'Registered machines:\n\n'
    printf '  %-12s %-7s %-30s %s\n' "NAME" "DEFAULT" "SUPPORTED" "LAST-VALIDATED"
    printf '  %-12s %-7s %-30s %s\n' "----" "-------" "---------" "--------------"
    local m detected
    detected="$(detect_machine)"
    for m in $MACHINE_REGISTRY; do
        reset_machine_meta
        "machine_${m}_meta"
        local mark=" "
        [[ -n "$detected" && "$detected" == "$m" ]] && mark="*"
        printf ' %s%-12s %-7s %-30s %s\n' "$mark" "$m" \
            "$MACHINE_META_DEFAULT_COMPILER" \
            "$MACHINE_META_SUPPORTED_COMPILERS" \
            "$MACHINE_META_LAST_VALIDATED"
    done
    printf '\n'
    [[ -n "$detected" ]] && printf '  * = auto-detected on this host (--machine=auto resolves to this)\n\n'
    printf 'For details on a single machine: --machine=NAME --dry-run\n'
}
if [[ "${_LIST_MACHINES:-no}" == "yes" ]]; then
    list_machines
    exit 0
fi

#-----------------------------------------------------------------------------
# Resolve --machine and --compiler
#-----------------------------------------------------------------------------
# Resolution order:
#   --machine=NAME       -> exact lookup, error if unknown
#   --machine=auto       -> detect_machine (matches first registered match() that returns 0)
#   <unset>              -> same as auto
# After resolution, the chosen machine's meta() is invoked, populating
# MACHINE_META_* globals. --compiler=NAME, if not passed, defaults to
# MACHINE_META_DEFAULT_COMPILER. The compiler family is informational in
# Push 1; Push 2 will use it to drive E3SM CIME env resolution.
apply_machine_defaults() {
    local requested="${MACHINE_NAME:-auto}"
    local resolved=""

    if [[ "$requested" == "auto" || -z "$requested" ]]; then
        resolved="$(detect_machine)"
        if [[ -z "$resolved" ]]; then
            # No registered match -- the orchestration banner will say so.
            MACHINE_NAME=""
            return 0
        fi
        MACHINE_NAME="$resolved"
    else
        if ! is_known_machine "$requested"; then
            die "unknown machine: $requested (try --list-machines)" 1
        fi
        MACHINE_NAME="$requested"
    fi

    reset_machine_meta
    "machine_${MACHINE_NAME}_meta"

    if [[ -z "$COMPILER_FAMILY" ]]; then
        COMPILER_FAMILY="$MACHINE_META_DEFAULT_COMPILER"
    fi

    # Validate compiler against supported list (warning only; user may know better)
    local supported="$MACHINE_META_SUPPORTED_COMPILERS"
    if [[ -n "$supported" && ",$supported," != *",$COMPILER_FAMILY,"* ]]; then
        warn "compiler '$COMPILER_FAMILY' not in $MACHINE_NAME's supported list ($supported); proceeding anyway"
    fi
}
apply_machine_defaults

#-----------------------------------------------------------------------------
# Profile resolution + E3SM environment loading
#-----------------------------------------------------------------------------
# --profile=e3sm        (default) Source modules + env vars from
#                       $E3SM_ROOT/cime_config/machines/config_machines.xml
#                       for the resolved machine + compiler. Requires that
#                       --e3sm-root or $E3SM_ROOT be set, and that
#                       MACHINE_META_E3SM_NAME is non-empty for the chosen
#                       machine entry.
# --profile=standalone  Trust whatever environment the user already has
#                       loaded. No CIME interaction. Same behavior as Push 1
#                       and prior, suitable for MOAB downstream users who
#                       don't have an E3SM checkout.
case "$PROFILE" in
    e3sm|standalone) ;;
    *) die "--profile must be one of: e3sm, standalone (got: $PROFILE)" 1 ;;
esac

# Adapter: some E3SM machines export aliases that install-moab.sh doesn't
# look for directly. Normalize so the rest of the script just sees the names
# it expects.
adapt_e3sm_env_to_install_moab_vars() {
    # NETCDF_C_PATH: bebop/improv/pm-cpu use NETCDF_C_PATH; anlgce-ub22 + crux
    # use NETCDF_PATH. If only the latter is set, alias it.
    if [[ -z "${NETCDF_C_PATH:-}" && -n "${NETCDF_PATH:-}" ]]; then
        export NETCDF_C_PATH="$NETCDF_PATH"
    fi
    # HDF5_ROOT: Cray PrgEnv sets CRAY_HDF5_PARALLEL_PREFIX via cray-hdf5-parallel.
    # MOAB / install-moab.sh wants HDF5_ROOT.
    if [[ -z "${HDF5_ROOT:-}" && -n "${CRAY_HDF5_PARALLEL_PREFIX:-}" ]]; then
        export HDF5_ROOT="$CRAY_HDF5_PARALLEL_PREFIX"
    fi
    # HDF5_ROOT sibling-derivation: bebop and improv configs do not export
    # HDF5_ROOT but DO export NETCDF_C_PATH like
    #     /lcrc/group/e3sm/soft/<machine>/netcdf-c/<ver>/<compiler>/<mpi>
    # The corresponding hdf5 install lives at the sibling
    #     /lcrc/group/e3sm/soft/<machine>/hdf5/<ver>/<compiler>/<mpi>
    # Sniff for it -- conservative: only adopt if the sibling exists AND has
    # include/hdf5.h (i.e. it's a real HDF5 install, not a coincidental dir).
    if [[ -z "${HDF5_ROOT:-}" && -n "${NETCDF_C_PATH:-}" ]]; then
        local derived
        derived="$(_derive_hdf5_root_from_netcdf "$NETCDF_C_PATH")"
        if [[ -n "$derived" ]]; then
            export HDF5_ROOT="$derived"
            log "Derived HDF5_ROOT from NETCDF_C_PATH sibling: $HDF5_ROOT"
        fi
    fi
}

# Given a NETCDF_C_PATH that follows the .../netcdf-c/<ver>/<rest> convention,
# look for a sibling .../hdf5/<*>/<rest> with include/hdf5.h. Prints the path
# on success, empty on failure. Conservative: tries exact-version-suffix match
# first (rare to align), then any-version match for the same compiler+mpi tail.
_derive_hdf5_root_from_netcdf() {
    local nc="$1"
    [[ "$nc" == */netcdf-c/* ]] || return 0
    # Split into parent (.../soft/<machine>) and tail (compiler/mpi or whatever
    # follows the netcdf version segment).
    local parent="${nc%/netcdf-c/*}"            # .../soft/<machine>
    local nc_tail="${nc#*/netcdf-c/}"           # <ver>/<compiler>/<mpi>
    local nc_after_ver="${nc_tail#*/}"          # <compiler>/<mpi>  (drop <ver>)
    [[ "$nc_after_ver" != "$nc_tail" ]] || return 0   # malformed: no <ver> segment

    local hdf5_parent="$parent/hdf5"
    [[ -d "$hdf5_parent" ]] || return 0

    # Iterate hdf5 versions present and pick the first one whose tail matches
    # the netcdf compiler/mpi suffix.
    local ver d
    for ver in "$hdf5_parent"/*; do
        [[ -d "$ver" ]] || continue
        d="$ver/$nc_after_ver"
        if [[ -f "$d/include/hdf5.h" ]]; then
            printf '%s' "$d"
            return 0
        fi
    done
    return 0
}

apply_e3sm_profile() {
    [[ "$PROFILE" == "e3sm" ]] || return 0

    # Resolve E3SM checkout
    if [[ -z "$E3SM_ROOT" ]]; then
        die "--profile=e3sm requires --e3sm-root=PATH or \$E3SM_ROOT (must point at an E3SM checkout containing cime_config/machines/config_machines.xml)" 2
    fi
    if [[ ! -f "$E3SM_ROOT/cime_config/machines/config_machines.xml" ]]; then
        die "--e3sm-root=$E3SM_ROOT is not an E3SM checkout (no cime_config/machines/config_machines.xml found)" 2
    fi

    # Need a machine entry with a non-empty E3SM_NAME mapping
    if [[ -z "${MACHINE_NAME:-}" ]]; then
        die "--profile=e3sm requires a known machine; pass --machine=NAME (see --list-machines) or run on a registered host" 2
    fi
    if [[ -z "${MACHINE_META_E3SM_NAME:-}" ]]; then
        die "machine '$MACHINE_NAME' has no E3SM config_machines.xml mapping (MACHINE_META_E3SM_NAME is empty); use --profile=standalone or extend the entry" 2
    fi

    # Locate the helper next to this script
    local helper
    helper="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts/e3sm_env.py"
    [[ -f "$helper" ]] || die "e3sm_env.py not found at $helper" 2

    log "Loading E3SM env: $MACHINE_META_E3SM_NAME / $COMPILER_FAMILY (via $helper)"

    # Generate to a temp file. Always run the helper (validates it works for
    # this machine + compiler), but only source the result for non-dry-run
    # invocations. Dry-run shouldn't pollute the calling shell's env, and
    # sourcing module loads on a host without the target's Lmod would fail.
    local tmp_env
    tmp_env="$(mktemp -t e3sm_env.XXXXXX.sh)"
    if ! python3 "$helper" \
            --e3sm-root="$E3SM_ROOT" \
            --machine="$MACHINE_META_E3SM_NAME" \
            --compiler="$COMPILER_FAMILY" \
            > "$tmp_env" 2> >(while IFS= read -r line; do warn "e3sm_env: $line"; done); then
        warn "e3sm_env.py exited non-zero; see $tmp_env for partial output"
        die "failed to derive E3SM environment for $MACHINE_META_E3SM_NAME / $COMPILER_FAMILY" 2
    fi

    log "Sourcing E3SM env snippet: $tmp_env"
    # The script runs in its own subshell, so sourcing here doesn't pollute the
    # user's parent shell. Sourcing under --dry-run is safe and necessary --
    # otherwise HDF5_ROOT / NETCDF_C_PATH / PNETCDF_PATH stay unset and the
    # downstream require_var validation fails before we can preview anything.
    # shellcheck disable=SC1090
    source "$tmp_env" || die "sourcing E3SM env snippet failed (see $tmp_env for the snippet content)" 2

    adapt_e3sm_env_to_install_moab_vars

    log "E3SM env applied. Resolved roots:"
    log "  HDF5_ROOT     = ${HDF5_ROOT:-(unset; pass --hdf5-root=PATH)}"
    log "  NETCDF_C_PATH = ${NETCDF_C_PATH:-(unset; pass --netcdf-root=PATH)}"
    log "  PNETCDF_PATH  = ${PNETCDF_PATH:-(unset; pass --pnetcdf-root=PATH)}"
    [[ -n "${BLAS_ROOT:-}"   ]] && log "  BLAS_ROOT     = $BLAS_ROOT"
    [[ -n "${LAPACK_ROOT:-}" ]] && log "  LAPACK_ROOT   = $LAPACK_ROOT"

    # Auto-populate --extra* with BLAS/LAPACK if user didn't set them and
    # the e3sm env exposed them. Lets the e3sm profile DTRT out of the box.
    if [[ -n "${BLAS_ROOT:-}" && -n "${LAPACK_ROOT:-}" ]]; then
        local blas_spec lapack_spec
        blas_spec="$(_resolve_blas_lapack_spec "$BLAS_ROOT" blas '-lgfortran')"
        lapack_spec="$(_resolve_blas_lapack_spec "$LAPACK_ROOT" lapack '-lm')"
        if [[ -z "$blas_spec" || -z "$lapack_spec" ]]; then
            warn "BLAS_ROOT=$BLAS_ROOT / LAPACK_ROOT=$LAPACK_ROOT did not resolve to a known library layout"
            warn "(probed: \$ROOT/, \$ROOT/lib/, \$ROOT/lib64/, \$ROOT/lib/intel64/ for libblas.a/liblapack.a/libmkl_*.a)"
            warn "Pass --extra and --extra-tempestremap manually, or extend _resolve_blas_lapack_spec"
        else
            local autodetect_extras="--with-blas=\"$blas_spec\" --with-lapack=\"$lapack_spec\""
            if [[ -z "$EXTRA_TEMPESTREMAP_ARGS" ]]; then
                EXTRA_TEMPESTREMAP_ARGS="$autodetect_extras"
                log "Auto-set --extra-tempestremap from BLAS_ROOT/LAPACK_ROOT"
            fi
            if [[ -z "$EXTRA_MOAB_ARGS" ]]; then
                EXTRA_MOAB_ARGS="$autodetect_extras"
                log "Auto-set --extra (MOAB) from BLAS_ROOT/LAPACK_ROOT"
            fi
        fi
    fi

    confirm_e3sm_env
}

# Resolve a BLAS/LAPACK library spec from a root path. Probes common layouts:
#   $ROOT/                   (Bebop pattern: libs directly in root, e.g.
#                             /lcrc/group/e3sm/soft/bebop/netlib-lapack/.../libblas.a)
#   $ROOT/lib/               (standard Unix prefix)
#   $ROOT/lib64/             (RHEL-style 64-bit)
#   $ROOT/lib/intel64/       (Intel MKL: $MKLROOT/lib/intel64/libmkl_intel_lp64.a)
#   $ROOT/lib/intel64_lin/   (older MKL variant)
#
# Args:
#   $1 = root path (e.g. $BLAS_ROOT, $LAPACK_ROOT, $MKLROOT)
#   $2 = library base name family: "blas" -> tries libblas.a then libmkl_intel_lp64.a;
#                                   "lapack" -> tries liblapack.a then libmkl_lapack95_lp64.a
#   $3 = trailing runtime libs to append (e.g. "-lgfortran" for blas, "-lm" for lapack)
#
# Prints (stdout): a single LIBS-style argument string suitable for AX_BLAS /
# AX_LAPACK's --with-blas=/--with-lapack=. Empty stdout if nothing found.
_resolve_blas_lapack_spec() {
    local root="$1" family="$2" runtime="$3"
    local sub d candidate

    # Build the (libname-list, libname-family) probe order. For MKL we'd need
    # multiple libs in the link line; treat that as a hint and emit the full
    # MKL combo if the marker file is found.
    local primary_names
    case "$family" in
        blas)   primary_names=("blas" "openblas" "mkl_intel_lp64" "sci_gnu_82_mp" "sci_gnu_82" "essl") ;;
        lapack) primary_names=("lapack" "openblas" "mkl_lapack95_lp64" "sci_gnu_82_mp" "sci_gnu_82") ;;
        *) printf ''; return 1 ;;
    esac

    # Static archive preferred (better for portability + avoids LD_LIBRARY_PATH gotchas)
    for sub in "" "/lib" "/lib64" "/lib/intel64" "/lib/intel64_lin"; do
        d="${root}${sub}"
        [[ -d "$d" ]] || continue
        for libname in "${primary_names[@]}"; do
            candidate="$d/lib${libname}.a"
            if [[ -f "$candidate" ]]; then
                # MKL family needs a multi-lib spec
                if [[ "$libname" == "mkl_intel_lp64" ]]; then
                    printf '%s/libmkl_intel_lp64.a %s/libmkl_sequential.a %s/libmkl_core.a -lpthread -lm -ldl' "$d" "$d" "$d"
                elif [[ "$libname" == "mkl_lapack95_lp64" ]]; then
                    # MKL has LAPACK in the same intel_lp64 trio
                    printf '%s/libmkl_intel_lp64.a %s/libmkl_sequential.a %s/libmkl_core.a -lpthread -lm -ldl' "$d" "$d" "$d"
                else
                    printf '%s %s' "$candidate" "$runtime"
                fi
                return 0
            fi
        done
    done

    # Fall back to shared (.so / .dylib) with -L/-l form
    for sub in "" "/lib" "/lib64" "/lib/intel64" "/lib/intel64_lin"; do
        d="${root}${sub}"
        [[ -d "$d" ]] || continue
        for libname in "${primary_names[@]}"; do
            if [[ -f "$d/lib${libname}.so" || -f "$d/lib${libname}.dylib" ]]; then
                printf -- '-L%s -l%s %s' "$d" "$libname" "$runtime"
                return 0
            fi
        done
    done

    printf ''
    return 1
}

# Show the resolved E3SM env and ask the user to confirm it's correct before
# proceeding. Bypass with --yes / -y / ASSUME_YES=yes. Refuse to proceed
# silently when stdin isn't a tty.
confirm_e3sm_env() {
    if [[ "$ASSUME_YES" == "yes" ]]; then
        log "Skipping confirmation prompt (--yes / ASSUME_YES=yes)"
        return 0
    fi
    if [[ ! -t 0 ]]; then
        die "stdin is not a tty -- pass --yes (or set ASSUME_YES=yes) to skip the e3sm env confirmation prompt in non-interactive mode" 1
    fi

    # Print to stderr (matches the rest of the banner) so prompt isn't lost in pipes.
    printf '\n' >&2
    printf '%s[install-moab]%s Does this E3SM environment look correct? [y/N] ' \
        "$COLOR_BLUE" "$COLOR_RESET" >&2

    local reply=""
    IFS= read -r reply || true
    case "$reply" in
        y|Y|yes|YES) log "Confirmed -- continuing" ;;
        *)
            log "Aborted at user request. To re-run with a different machine/compiler:"
            log "  install-moab.sh --machine=NAME --compiler=NAME --profile=e3sm --e3sm-root=$E3SM_ROOT [...]"
            log "Or pass --profile=standalone to skip env loading entirely."
            exit 0
            ;;
    esac
}
apply_e3sm_profile

# Recompute TPL_PREFIX default if --prefix changed and TPL_PREFIX wasn't set
if [[ "$TPL_PREFIX" == "$HOME/install/MOAB/tpls" && "$PREFIX_PATH" != "$HOME/install/MOAB" ]]; then
    TPL_PREFIX="$PREFIX_PATH/tpls"
fi

# Resolve BUILD_DIR / MOAB_SRC_DIR defaults now that arg parsing is complete.
# Two scenarios:
#   1. User passed --src-dir (private dev checkout): BUILD_DIR defaults to
#      $MOAB_SRC_DIR/build-e3sm so artifacts live next to the source.
#   2. Otherwise (script-managed source): BUILD_DIR defaults to $PWD/moab-build
#      and MOAB_SRC_DIR defaults to $BUILD_DIR/moab-src.
if [[ "${_USER_SET_BUILD_DIR:-no}" != "yes" ]]; then
    if [[ "${_USER_SET_MOAB_SRC:-no}" == "yes" ]]; then
        BUILD_DIR="$MOAB_SRC_DIR/build-e3sm"
    else
        BUILD_DIR="$PWD/moab-build"
    fi
fi
if [[ "${_USER_SET_MOAB_SRC:-no}" != "yes" ]]; then
    MOAB_SRC_DIR="$BUILD_DIR/moab-src"
fi

case "$BUILD_SYSTEM" in
    autotools|cmake) ;;
    *) die "--build-system must be autotools or cmake (got: $BUILD_SYSTEM)" 1 ;;
esac

#-----------------------------------------------------------------------------
# Validation
#-----------------------------------------------------------------------------
# Up-front tooling checks. The MOAB source tree itself is validated later,
# after prepare_moab_source has had a chance to clone it.
if [[ "$BUILD_SYSTEM" == "cmake" ]]; then
    command -v cmake >/dev/null 2>&1 || die "cmake not found on PATH (required by --build-system=cmake)" 2
fi

require_var() {
    local name="$1"
    local val="${!name:-}"
    [[ -n "$val" ]] || die "required env var $name is not set (load the E3SM env first or pass via flag)" 2
    [[ -d "$val" ]] || warn "$name=$val does not point to an existing directory"
}
require_var HDF5_ROOT
require_var NETCDF_C_PATH
require_var PNETCDF_PATH

# MPI_ROOT is OPTIONAL. On Cray PrgEnv (and similar) the wrappers are cc/CC/ftn
# with no meaningful root. If MPI_ROOT is set, we use it to resolve the wrapper
# binaries. If not, we rely on the user's PATH and explicit --cc/--cxx/--fc/--f77.
#
# Empty value (e.g. user invoked `--mpi-root=$MPI_ROOT` where $MPI_ROOT was
# unset in their shell) is treated as "not set" rather than "MPI is at /".
if [[ -n "${MPI_ROOT+x}" ]] && [[ -z "$MPI_ROOT" ]]; then
    warn "MPI_ROOT is set but empty (did you forget to load the MPI module before invoking the script?). Treating as unset."
    unset MPI_ROOT
fi
if [[ -n "${MPI_ROOT:-}" ]]; then
    [[ -d "$MPI_ROOT" ]] || warn "MPI_ROOT=$MPI_ROOT does not point to an existing directory"
fi

# Resolve a single compiler. Priority (first match wins):
#   1. --cc/--cxx/--fc/--f77 flag     (binding)
#   2. $MPI_ROOT/bin/<wrapper>          (binding when MPI_ROOT is set)
#   3. $MPI_ROOT/bin/<alternate>
#   4. PATH lookup of <wrapper>         (preferred over env-set non-MPI value)
#   5. PATH lookup of <alternate>
#   6. env-set CC/CXX/FC/F77            (last-resort hint; may fail validation)
#   7. die
#
# Why env CC is demoted: on Spack-managed systems, `module load gcc-12` sets
# CC=/path/to/gcc, which is NOT an MPI wrapper. Users still expect the script
# to pick up mpicc when MPI is available. Treating env CC as a hard override
# silently bypassed MPI auto-detection. See incident: GCE login node 2026-05.
resolve_compiler() {
    local var="$1" flag_set="$2" env_set="$3" default="$4"; shift 4
    local alternates=( "$@" )
    local val found a

    # 1. Explicit --flag wins absolutely
    if [[ "$flag_set" == "yes" ]]; then
        val="${!var}"
        if [[ "$val" == /* ]]; then
            [[ -x "$val" ]] || die "user-specified $var=$val not executable" 2
        elif command -v "$val" >/dev/null 2>&1; then
            printf -v "$var" '%s' "$(command -v "$val")"
        else
            die "user-specified $var=$val not found on PATH" 2
        fi
        return 0
    fi

    # 2-3. MPI_ROOT-based lookup
    if [[ -n "${MPI_ROOT:-}" ]]; then
        if [[ -x "$MPI_ROOT/bin/$default" ]]; then
            printf -v "$var" '%s' "$MPI_ROOT/bin/$default"
            return 0
        fi
        for a in "${alternates[@]}"; do
            if [[ -x "$MPI_ROOT/bin/$a" ]]; then
                printf -v "$var" '%s' "$MPI_ROOT/bin/$a"
                return 0
            fi
        done
    fi

    # 4-5. PATH lookup of MPI wrapper names (overrides env-set non-MPI value)
    if command -v "$default" >/dev/null 2>&1; then
        found="$(command -v "$default")"
        if [[ "$env_set" == "yes" ]] && [[ "${!var}" != "$found" ]]; then
            log "$var: using $found from PATH (overriding env $var=${!var}; MPI wrapper takes precedence)"
        fi
        printf -v "$var" '%s' "$found"
        return 0
    fi
    for a in "${alternates[@]}"; do
        if command -v "$a" >/dev/null 2>&1; then
            found="$(command -v "$a")"
            if [[ "$env_set" == "yes" ]] && [[ "${!var}" != "$found" ]]; then
                log "$var: using $found from PATH (overriding env $var=${!var})"
            fi
            printf -v "$var" '%s' "$found"
            return 0
        fi
    done

    # 6. Last resort: env-set value (warn -- it almost certainly isn't MPI)
    if [[ "$env_set" == "yes" ]]; then
        val="${!var}"
        warn "$var: no MPI wrapper found via MPI_ROOT or PATH; falling back to env $var=$val (likely will fail MPI validation)"
        if [[ "$val" == /* ]]; then
            [[ -x "$val" ]] || die "env $var=$val not executable" 2
        elif command -v "$val" >/dev/null 2>&1; then
            printf -v "$var" '%s' "$(command -v "$val")"
        else
            die "env $var=$val not found on PATH" 2
        fi
        return 0
    fi

    die "could not resolve $var (tried: $default ${alternates[*]:-}; MPI_ROOT=${MPI_ROOT:-unset}; load MPI module or pass --cc/--cxx/--fc/--f77)" 2
}

# Detect Cray PrgEnv. Any of these signals is sufficient:
#   * PE_ENV         -- set by `module load PrgEnv-*` (GNU/CRAY/INTEL/AOCC/...)
#   * CRAYPE_VERSION -- set whenever craype is loaded
#   * /opt/cray      -- present on every Cray-based system
detect_cray_env() {
    [[ -n "${PE_ENV:-}" ]] || [[ -n "${CRAYPE_VERSION:-}" ]] || [[ -d /opt/cray ]]
}

# Pick the search order. On Cray (without a user-supplied --mpi-root) the
# canonical MPI wrappers are cc/CC/ftn, not mpicc/mpicxx/mpif90 -- prefer them.
if detect_cray_env && [[ -z "${MPI_ROOT:-}" ]]; then
    log "Cray PrgEnv detected (PE_ENV=${PE_ENV:-unset} CRAYPE_VERSION=${CRAYPE_VERSION:-unset}); preferring cc/CC/ftn wrappers"
    resolve_compiler CC_BIN  "${_FLAG_HAS_CC:-no}"  "${_ENV_HAS_CC:-no}"  cc   mpicc
    resolve_compiler CXX_BIN "${_FLAG_HAS_CXX:-no}" "${_ENV_HAS_CXX:-no}" CC   mpicxx mpic++ mpiCC
    resolve_compiler FC_BIN  "${_FLAG_HAS_FC:-no}"  "${_ENV_HAS_FC:-no}"  ftn  mpif90 mpifort mpiifort
    if ! resolve_compiler F77_BIN "${_FLAG_HAS_F77:-no}" "${_ENV_HAS_F77:-no}" ftn mpif77 mpifort mpiifort 2>/dev/null; then
        warn "F77 wrapper not found; reusing FC=$FC_BIN as F77"
        F77_BIN="$FC_BIN"
    fi
else
    resolve_compiler CC_BIN  "${_FLAG_HAS_CC:-no}"  "${_ENV_HAS_CC:-no}"  mpicc   cc
    resolve_compiler CXX_BIN "${_FLAG_HAS_CXX:-no}" "${_ENV_HAS_CXX:-no}" mpicxx  mpic++ mpiCC CC
    resolve_compiler FC_BIN  "${_FLAG_HAS_FC:-no}"  "${_ENV_HAS_FC:-no}"  mpif90  mpifort mpiifort ftn
    if ! resolve_compiler F77_BIN "${_FLAG_HAS_F77:-no}" "${_ENV_HAS_F77:-no}" mpif77 mpifort mpiifort ftn 2>/dev/null; then
        warn "F77 wrapper not found; reusing FC=$FC_BIN as F77"
        F77_BIN="$FC_BIN"
    fi
fi

DOWNLOADER=""
if   command -v curl >/dev/null 2>&1; then DOWNLOADER=curl
elif command -v wget >/dev/null 2>&1; then DOWNLOADER=wget
fi

#-----------------------------------------------------------------------------
# MPI wrapper validation (compile + link a tiny MPI program)
#
# This catches misnamed wrappers, broken modules, missing mpi.h or mpif.h,
# and PrgEnv mismatches BEFORE we spend hours building TPLs that would fail
# the same way deep inside their own configure.
#-----------------------------------------------------------------------------
validate_mpi_wrapper() {
    local lang="$1" compiler="$2"
    local tmpdir; tmpdir="$(mktemp -d)"
    local rc=0 src bin

    case "$lang" in
        c)
            src="$tmpdir/test.c"; bin="$tmpdir/test_c"
            cat > "$src" <<'__SRC__'
#include <mpi.h>
#include <stdio.h>
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Finalize();
    return 0;
}
__SRC__
            ;;
        cxx)
            src="$tmpdir/test.cpp"; bin="$tmpdir/test_cxx"
            cat > "$src" <<'__SRC__'
#include <mpi.h>
#include <iostream>
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Finalize();
    return 0;
}
__SRC__
            ;;
        fortran)
            # Fixed-form to satisfy both F77 and modern Fortran compilers.
            src="$tmpdir/test.f"; bin="$tmpdir/test_f"
            cat > "$src" <<'__SRC__'
      program t
      implicit none
      include 'mpif.h'
      integer ierr, rank
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Finalize(ierr)
      end
__SRC__
            ;;
        *) rm -rf "$tmpdir"; return 1 ;;
    esac

    "$compiler" -o "$bin" "$src" >"$tmpdir/log" 2>&1 || rc=$?
    if (( rc != 0 )) || [[ ! -x "$bin" ]]; then
        warn "MPI $lang validation FAILED: $compiler"
        warn "  command: $compiler -o $bin $src"
        warn "  output:"
        sed 's/^/    /' "$tmpdir/log" >&2 || true
        rm -rf "$tmpdir"
        return 1
    fi
    rm -rf "$tmpdir"
    return 0
}

validate_all_mpi_wrappers() {
    section "MPI wrapper validation"
    if [[ "$SKIP_MPI_VALIDATION" == "yes" ]]; then
        warn "skipping MPI validation (SKIP_MPI_VALIDATION=yes)"
        return 0
    fi
    local fail=0
    log "  C   : $CC_BIN"
    validate_mpi_wrapper c       "$CC_BIN"  && ok "  C   wrapper OK" || fail=1
    log "  C++ : $CXX_BIN"
    validate_mpi_wrapper cxx     "$CXX_BIN" && ok "  C++ wrapper OK" || fail=1
    log "  FC  : $FC_BIN"
    if validate_mpi_wrapper fortran "$FC_BIN"; then
        ok "  FC  wrapper OK"
    else
        warn "  FC  wrapper failed -- MOAB/TempestRemap will likely fail to link MPI Fortran"
        fail=1
    fi
    if [[ "$F77_BIN" != "$FC_BIN" ]]; then
        log "  F77 : $F77_BIN"
        if validate_mpi_wrapper fortran "$F77_BIN"; then
            ok "  F77 wrapper OK"
        else
            warn "  F77 wrapper failed -- you can pass --f77=$FC_BIN as a workaround"
            fail=1
        fi
    fi
    (( fail == 0 )) || die "MPI wrapper validation failed (use --skip-mpi-validation to bypass at your own risk)" 2
}

#-----------------------------------------------------------------------------
# MOAB source management
#
# Default behavior (script-managed source):
#   * If $MOAB_SRC_DIR/.git is absent -> shallow clone $MOAB_REPO_URL@$MOAB_BRANCH
#   * If present -> git fetch + git reset --hard origin/$MOAB_BRANCH (forced FF)
#   * For autotools, run autoreconf -fi if $MOAB_SRC_DIR/configure is missing
#
# User-managed source (--src-dir explicitly set):
#   * No git operations. The directory must already exist with a usable tree.
#
# Opt-out for script-managed (--no-source-update):
#   * Skip git fetch/reset; require directory to already exist.
#-----------------------------------------------------------------------------
prepare_moab_source() {
    section "MOAB source ($MOAB_SRC_DIR)"

    if [[ "${_USER_SET_MOAB_SRC:-no}" == "yes" ]]; then
        log "User-managed source dir (--src-dir) -- no git operations"
        [[ -d "$MOAB_SRC_DIR" ]] || die "user-specified --src-dir not found: $MOAB_SRC_DIR" 1
    elif [[ "$NO_SOURCE_UPDATE" == "yes" ]]; then
        log "--no-source-update: leaving $MOAB_SRC_DIR untouched"
        [[ -d "$MOAB_SRC_DIR" ]] || die "MOAB source dir does not exist and --no-source-update was set: $MOAB_SRC_DIR" 1
    else
        command -v git >/dev/null 2>&1 || die "git is required for source management (or pass --src-dir to use an existing checkout)" 2
        if [[ -d "$MOAB_SRC_DIR/.git" ]]; then
            log "Updating MOAB clone to latest origin/$MOAB_BRANCH"
            ( cd "$MOAB_SRC_DIR" \
                && git fetch --quiet --tags origin "$MOAB_BRANCH" \
                && git checkout --quiet "$MOAB_BRANCH" \
                && git reset --hard --quiet "origin/$MOAB_BRANCH" ) \
                || die "git update failed in $MOAB_SRC_DIR" 3
        elif [[ -d "$MOAB_SRC_DIR" ]] && [[ -n "$(ls -A "$MOAB_SRC_DIR" 2>/dev/null)" ]]; then
            die "$MOAB_SRC_DIR exists but is not a git clone; refusing to nuke it (use --src-dir to point at an existing tree, or remove this directory)" 1
        else
            log "Cloning MOAB from $MOAB_REPO_URL ($MOAB_BRANCH) into $MOAB_SRC_DIR"
            mkdir -p "$(dirname "$MOAB_SRC_DIR")"
            rm -rf "$MOAB_SRC_DIR"   # may exist as empty dir
            git clone --quiet --branch "$MOAB_BRANCH" --single-branch \
                "$MOAB_REPO_URL" "$MOAB_SRC_DIR" \
                || die "git clone failed: $MOAB_REPO_URL" 3
        fi
    fi

    # Post-clone tree validation
    case "$BUILD_SYSTEM" in
        autotools)
            if [[ ! -x "$MOAB_SRC_DIR/configure" ]]; then
                command -v autoreconf >/dev/null 2>&1 \
                    || die "configure missing and autoreconf not on PATH (install autotools)" 2
                log "Bootstrapping with autoreconf -fi"
                ( cd "$MOAB_SRC_DIR" && autoreconf -fi >/dev/null 2>&1 ) \
                    || die "autoreconf -fi failed" 3
            fi
            [[ -x "$MOAB_SRC_DIR/configure" ]] \
                || die "configure script still missing after autoreconf" 3
            ;;
        cmake)
            [[ -f "$MOAB_SRC_DIR/CMakeLists.txt" ]] \
                || die "CMakeLists.txt missing at $MOAB_SRC_DIR" 2
            ;;
    esac

    if [[ -d "$MOAB_SRC_DIR/.git" ]]; then
        local short="$(cd "$MOAB_SRC_DIR" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
        local subj="$(cd "$MOAB_SRC_DIR" && git log -1 --pretty=%s 2>/dev/null || echo '(no log)')"
        ok "MOAB source ready: HEAD=$short  ($subj)"
    else
        ok "MOAB source ready: $MOAB_SRC_DIR (non-git checkout)"
    fi
}

# Helper: short HEAD hash of MOAB src (or "no-git" if not a git checkout).
moab_src_hash() {
    if [[ -d "$MOAB_SRC_DIR/.git" ]]; then
        ( cd "$MOAB_SRC_DIR" && git rev-parse HEAD 2>/dev/null ) || printf 'unknown'
    else
        printf 'no-git'
    fi
}

#-----------------------------------------------------------------------------
# Layout
#-----------------------------------------------------------------------------
ARCHIVES_DIR="$BUILD_DIR/tpl-archives"
TPL_WORK_DIR="$BUILD_DIR/tpl-work"
TPL_LOG_DIR="$BUILD_DIR/tpl-logs"
MOAB_BUILD_DIR="$BUILD_DIR/moab"

#-----------------------------------------------------------------------------
# TPL metadata (Bash 3.2-compatible)
#-----------------------------------------------------------------------------
# Source-fetch metadata. Tempestremap is git-clone-first (shallow), with
# the GitHub branch tarball as a non-git fallback. Zoltan and Eigen3 are
# tarball-only.
tpl_url() {
    case "$1" in
        eigen3)       printf '%s/eigen-3_4_0.tar.gz' "$TPL_MIRROR" ;;
        zoltan)       printf '%s/zoltan-3_9_1.tar.gz' "$TPL_MIRROR" ;;
        tempestremap) printf 'https://github.com/E3SM-Project/tempestremap/archive/refs/heads/master.tar.gz' ;;
    esac
}
tpl_fallback_url() {
    case "$1" in
        eigen3)       printf 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz' ;;
        zoltan)       printf 'https://github.com/sandialabs/Zoltan/archive/refs/tags/v3.901.tar.gz' ;;
        tempestremap) printf '' ;;   # only the GitHub source is allowed
    esac
}
tpl_git_repo() {
    case "$1" in
        tempestremap) printf 'https://github.com/E3SM-Project/tempestremap.git' ;;
        *)            printf '' ;;
    esac
}
tpl_git_branch() {
    case "$1" in
        tempestremap) printf 'master' ;;
        *)            printf '' ;;
    esac
}
tpl_install_dir()  { printf '%s/%s' "$TPL_PREFIX" "$1"; }
tpl_src_dir()      { printf '%s/%s/src' "$TPL_WORK_DIR" "$1"; }
tpl_build_dir()    { printf '%s/%s/build' "$TPL_WORK_DIR" "$1"; }
tpl_log()          { printf '%s/%s_%s.log' "$TPL_LOG_DIR" "$2" "$1"; }   # tpl, phase

# A TPL is "already installed" if the marker file below exists.
tpl_install_marker() {
    case "$1" in
        eigen3)       printf '%s/include/eigen3/Eigen/Dense' "$(tpl_install_dir eigen3)" ;;
        zoltan)       printf '%s/include/zoltan.h' "$(tpl_install_dir zoltan)" ;;
        tempestremap) printf '%s/include/TempestRemapAPI.h' "$(tpl_install_dir tempestremap)" ;;
    esac
}

#-----------------------------------------------------------------------------
# Download with retry + mirror fallback
#-----------------------------------------------------------------------------
download_one() {
    local url="$1" out="$2"
    if [[ "$DOWNLOADER" == "curl" ]]; then
        curl -fSL --connect-timeout 30 -o "$out.part" "$url" && mv "$out.part" "$out"
    else
        wget -q -O "$out.part" "$url" && mv "$out.part" "$out"
    fi
}

fetch_tpl_archive() {
    local name="$1"
    local url="$(tpl_url "$name")"
    local fallback="$(tpl_fallback_url "$name")"
    local archive="$ARCHIVES_DIR/$(basename "$url")"

    if [[ -s "$archive" ]] && tar -tzf "$archive" >/dev/null 2>&1; then
        log "$name: archive cached and valid ($archive)"
        printf '%s' "$archive"; return 0
    fi
    rm -f "$archive" "$archive.part"

    [[ -n "$DOWNLOADER" ]] || { warn "no downloader (curl/wget) on PATH"; return 1; }

    local attempt=1 src_url="$url"
    while (( attempt <= DOWNLOAD_RETRIES )); do
        log "$name: download attempt $attempt/$DOWNLOAD_RETRIES from $src_url"
        if download_one "$src_url" "$archive" \
           && tar -tzf "$archive" >/dev/null 2>&1; then
            ok "$name: archive downloaded ($archive)"
            printf '%s' "$archive"; return 0
        fi
        rm -f "$archive" "$archive.part"
        # Switch to fallback halfway through, only if a fallback exists
        if (( attempt == DOWNLOAD_RETRIES / 2 )) \
           && [[ -n "$fallback" && "$src_url" != "$fallback" ]]; then
            warn "$name: switching to fallback URL ($fallback)"
            src_url="$fallback"
            archive="$ARCHIVES_DIR/$(basename "$src_url")"
        fi
        sleep "$DOWNLOAD_RETRY_DELAY"
        (( attempt++ ))
    done
    return 1
}

extract_tpl_archive() {
    local name="$1" archive="$2"
    local srcdir="$(tpl_src_dir "$name")"
    if [[ -d "$srcdir" ]]; then
        log "$name: removing stale/partial src dir before re-extraction"
        rm -rf "$srcdir"
    fi
    mkdir -p "$srcdir"
    log "$name: extracting $archive into $srcdir"
    tar -xzf "$archive" -C "$srcdir" --strip-components=1
}

# Shallow git clone with retries. Returns 0 on success.
clone_tpl_repo() {
    local name="$1"
    local repo="$(tpl_git_repo "$name")"
    local branch="$(tpl_git_branch "$name")"
    local srcdir="$(tpl_src_dir "$name")"
    [[ -n "$repo" && -n "$branch" ]] || return 1
    command -v git >/dev/null 2>&1 || return 1

    [[ -d "$srcdir" ]] && rm -rf "$srcdir"
    mkdir -p "$(dirname "$srcdir")"

    local attempt=1
    while (( attempt <= DOWNLOAD_RETRIES )); do
        log "$name: shallow clone attempt $attempt/$DOWNLOAD_RETRIES from $repo (branch: $branch)"
        if git clone --depth 1 --branch "$branch" --single-branch "$repo" "$srcdir" >&2; then
            ok "$name: cloned ($srcdir)"
            return 0
        fi
        rm -rf "$srcdir"
        sleep "$DOWNLOAD_RETRY_DELAY"
        (( attempt++ ))
    done
    return 1
}

# Single entry point. Guarantees that on return-0:
#   * $(tpl_src_dir <name>) contains an unpacked source tree
#   * $(tpl_src_dir <name>)/.install-moab.source_ok exists
prepare_tpl_source() {
    local name="$1"
    local srcdir="$(tpl_src_dir "$name")"
    local marker="$srcdir/.install-moab.source_ok"

    if [[ -f "$marker" ]]; then
        log "$name: source already prepared ($srcdir)"
        return 0
    fi

    # Git-clone path (currently only tempestremap). Falls back to tarball
    # if git is unavailable or the clone fails outright.
    if [[ -n "$(tpl_git_repo "$name")" ]]; then
        if clone_tpl_repo "$name"; then
            : > "$marker"
            return 0
        fi
        warn "$name: git clone path failed; trying tarball download"
    fi

    # Tarball path
    local archive
    archive="$(fetch_tpl_archive "$name")" || return 1
    extract_tpl_archive "$name" "$archive" || return 1
    : > "$marker"
    return 0
}

#-----------------------------------------------------------------------------
# Live tail of per-TPL logs
#-----------------------------------------------------------------------------
TAIL_PIDS=()
stop_tails() {
    [[ ${#TAIL_PIDS[@]} -eq 0 ]] && return 0
    local pid
    for pid in "${TAIL_PIDS[@]}"; do kill "$pid" >/dev/null 2>&1 || true; done
    TAIL_PIDS=()
}
trap 'stop_tails' EXIT INT TERM

start_tail() {
    [[ "$TAIL_LOGS" == "yes" ]] || return 0
    local f="$1" tag="$2"
    touch "$f"
    ( tail -n 0 -F "$f" 2>/dev/null \
        | awk -v tag="$tag" '{ printf "  [%s] %s\n", tag, $0; fflush() }' ) &
    TAIL_PIDS+=( $! )
}

#-----------------------------------------------------------------------------
# Failure classification + remediation (per-TPL scoped)
#-----------------------------------------------------------------------------
# Each remediation is recorded as additive flags. We rebuild the TPL's
# configure command from these on each retry.
TPL_EXTRA_CFLAGS=""
TPL_EXTRA_CXXFLAGS=""
TPL_EXTRA_LDFLAGS=""

reset_tpl_remediations() {
    TPL_EXTRA_CFLAGS=""
    TPL_EXTRA_CXXFLAGS=""
    TPL_EXTRA_LDFLAGS=""
}

classify_log() {
    local log="$1"
    [[ -f "$log" ]] || return 0
    if grep -qE 'fatal: unable to access|Could not resolve host|Connection (timed out|refused)|Network is unreachable|wget: unable to resolve' "$log"; then
        printf 'network'; return 0
    fi
    if grep -qE 'multiple definition of|undefined reference to.*Zoltan_|relocation R_X86_64.*against `[A-Za-z_].*'\'' can not be used' "$log"; then
        printf 'fcommon'; return 0
    fi
    if grep -qE 'configure: No such file or directory|aclocal.*not found|cannot find install-sh' "$log"; then
        printf 'autoreconf'; return 0
    fi
    if grep -qE 'netcdf\.h.*No such file|cannot find -lnetcdf|Cannot find NetCDF' "$log"; then
        printf 'netcdfroot'; return 0
    fi
    if grep -qE 'requires (C\+\+1[14])|error: .*[Cc]\+\+1[14]|range-based .for. loops are not allowed' "$log"; then
        printf 'cxx14'; return 0
    fi
    if grep -qE 'cannot find -llapack|undefined reference to .dsyev_|undefined reference to .dgemm_' "$log"; then
        printf 'lapack'; return 0
    fi
    if grep -qE 'Permission denied|Read-only file system|No space left on device' "$log"; then
        printf 'fatal_filesystem'; return 0
    fi
    printf ''
}

# Apply a remediation in the current TPL's scope. Returns 1 if not actionable.
apply_remediation() {
    local tpl="$1" tag="$2"
    case "$tag" in
        network)
            warn "$tpl: network failure detected; will refresh archive"
            rm -f "$ARCHIVES_DIR/$(basename "$(tpl_url "$tpl")")"
            ;;
        fcommon)
            warn "$tpl: applying CFLAGS=-fcommon (GCC>=10 multiple-definition fix)"
            TPL_EXTRA_CFLAGS="$TPL_EXTRA_CFLAGS -fcommon"
            ;;
        autoreconf)
            warn "$tpl: re-bootstrapping with autoreconf -fi"
            local sd="$(tpl_src_dir "$tpl")"
            ( cd "$sd" && autoreconf -fi >/dev/null 2>&1 ) \
                || { warn "$tpl: autoreconf -fi failed"; return 1; }
            ;;
        netcdfroot)
            warn "$tpl: exporting NETCDFROOT/NETCDF_DIR/NETCDF_PATH"
            export NETCDFROOT="$NETCDF_C_PATH"
            export NETCDF_DIR="$NETCDF_C_PATH"
            export NETCDF_PATH="$NETCDF_C_PATH"
            ;;
        cxx14)
            warn "$tpl: forcing -std=c++14"
            TPL_EXTRA_CXXFLAGS="$TPL_EXTRA_CXXFLAGS -std=c++14"
            ;;
        lapack)
            warn "$tpl: appending -llapack -lblas to LDFLAGS"
            TPL_EXTRA_LDFLAGS="$TPL_EXTRA_LDFLAGS -llapack -lblas"
            ;;
        fatal_filesystem)
            warn "$tpl: filesystem error -- not auto-recoverable"
            return 1
            ;;
        *) return 1 ;;
    esac
    return 0
}

# Run a phase command, capturing to the per-TPL log, and tail it live.
run_phase() {
    local tpl="$1" phase="$2"; shift 2
    local logf="$(tpl_log "$tpl" "$phase")"
    : > "$logf"
    start_tail "$logf" "${tpl}-${phase}"
    log "$tpl: $phase ..."
    set +e
    ( "$@" ) >>"$logf" 2>&1
    local rc=$?
    set -e
    stop_tails
    return $rc
}

#-----------------------------------------------------------------------------
# Per-TPL recipes
#-----------------------------------------------------------------------------
# Each recipe is invoked inside build_tpl_with_retries(), which:
#   1. fetches the archive (with retry)
#   2. extracts it
#   3. calls <recipe>_configure / <recipe>_build / <recipe>_install
#   4. on failure, classifies the log and retries up to TPL_MAX_ATTEMPTS

#---------- Eigen3 (header-only) ----------
eigen3_configure() { :; }   # no-op
eigen3_build()     { :; }   # no-op
eigen3_install() {
    local prefix="$(tpl_install_dir eigen3)"
    local srcdir="$(tpl_src_dir eigen3)"
    mkdir -p "$prefix/include/eigen3"
    cp -R "$srcdir/Eigen"  "$prefix/include/eigen3/" \
        || { warn "eigen3: failed copying Eigen/ headers"; return 1; }
    [[ -d "$srcdir/unsupported" ]] && cp -R "$srcdir/unsupported" "$prefix/include/eigen3/" || true
    [[ -f "$srcdir/signature_of_eigen3_matrix_library" ]] \
        && cp "$srcdir/signature_of_eigen3_matrix_library" "$prefix/include/eigen3/" || true
    return 0
}

#---------- Zoltan (autotools) ----------
zoltan_configure() {
    local prefix="$(tpl_install_dir zoltan)"
    local srcdir="$(tpl_src_dir zoltan)"
    local build="$(tpl_build_dir zoltan)"
    mkdir -p "$build"

    # Some Zoltan tarballs ship without configure (git snapshot) -- bootstrap
    if [[ ! -f "$srcdir/configure" ]]; then
        log "zoltan: configure missing; running autoreconf -fi"
        ( cd "$srcdir" && autoreconf -fi >>"$(tpl_log zoltan config)" 2>&1 ) \
            || return 1
    fi

    local cflags="$TPL_EXTRA_CFLAGS"
    local cxxflags="$TPL_EXTRA_CXXFLAGS"
    local ldflags="$TPL_EXTRA_LDFLAGS"
    # User-supplied extras (--extra-zoltan / EXTRA_ZOLTAN_ARGS).
    # `eval` re-parses the string so embedded shell quoting is honored, e.g.
    #     --extra-zoltan='LIBS="-L/path -lopenblas -lgfortran"'
    # produces a single argv entry LIBS=-L/path -lopenblas -lgfortran.
    local extra_user=()
    [[ -n "$EXTRA_ZOLTAN_ARGS" ]] && eval "extra_user=( $EXTRA_ZOLTAN_ARGS )"
    ( cd "$build" && \
        "$srcdir/configure" \
            --prefix="$prefix" \
            --libdir="$prefix/lib" \
            --with-pic=1 \
            --enable-mpi \
            --enable-shared="$ENABLE_SHARED" \
            --enable-static="$ENABLE_STATIC" \
            CC="$CC_BIN" CXX="$CXX_BIN" FC="$FC_BIN" F77="$F77_BIN" \
            ${cflags:+CFLAGS="$cflags"} \
            ${cxxflags:+CXXFLAGS="$cxxflags"} \
            ${ldflags:+LDFLAGS="$ldflags"} \
            ${extra_user[@]+"${extra_user[@]}"} )
}

zoltan_build() {
    local build="$(tpl_build_dir zoltan)"
    make -C "$build" -j"$JOBS"
}

zoltan_install() {
    local build="$(tpl_build_dir zoltan)"
    make -C "$build" install
}

#---------- TempestRemap (autotools) ----------
tempestremap_configure() {
    local prefix="$(tpl_install_dir tempestremap)"
    local srcdir="$(tpl_src_dir tempestremap)"
    local build="$(tpl_build_dir tempestremap)"
    mkdir -p "$build"

    if [[ ! -f "$srcdir/configure" ]]; then
        log "tempestremap: configure missing; running autoreconf -fi"
        ( cd "$srcdir" && autoreconf -fi >>"$(tpl_log tempestremap config)" 2>&1 ) \
            || return 1
    fi

    local cflags="$TPL_EXTRA_CFLAGS"
    local cxxflags="$TPL_EXTRA_CXXFLAGS"
    local ldflags="$TPL_EXTRA_LDFLAGS"
    # User-supplied extras (--extra-tempestremap / EXTRA_TEMPESTREMAP_ARGS).
    # `eval` re-parses the string so embedded shell quoting is honored. Use this
    # for AX_BLAS / AX_LAPACK selection (the values must each be a single argv
    # entry containing the LIBS-style spec):
    #     --extra-tempestremap='--with-blas="/path/libblas.a -lgfortran" \
    #                           --with-lapack="/path/liblapack.a -lm"'
    local extra_user=()
    [[ -n "$EXTRA_TEMPESTREMAP_ARGS" ]] && eval "extra_user=( $EXTRA_TEMPESTREMAP_ARGS )"
    ( cd "$build" && \
        "$srcdir/configure" \
            --prefix="$prefix" \
            --libdir="$prefix/lib" \
            --with-pic=1 \
            --enable-shared="$ENABLE_SHARED" \
            --enable-static="$ENABLE_STATIC" \
            --with-netcdf="$NETCDF_C_PATH" \
            --with-hdf5="$HDF5_ROOT" \
            CC="$CC_BIN" CXX="$CXX_BIN" FC="$FC_BIN" F77="$F77_BIN" \
            ${cflags:+CFLAGS="$cflags"} \
            ${cxxflags:+CXXFLAGS="$cxxflags"} \
            ${ldflags:+LDFLAGS="$ldflags"} \
            ${extra_user[@]+"${extra_user[@]}"} )
}

tempestremap_build() {
    local build="$(tpl_build_dir tempestremap)"
    make -C "$build" -j"$JOBS"
}

tempestremap_install() {
    local build="$(tpl_build_dir tempestremap)"
    make -C "$build" install
}

#-----------------------------------------------------------------------------
# Per-TPL user extras (configure args supplied via --extra-<tpl>)
#-----------------------------------------------------------------------------
# Returns the literal extras string for $1; empty for TPLs that don't accept any.
tpl_extra_args() {
    case "$1" in
        zoltan)       printf '%s' "$EXTRA_ZOLTAN_ARGS" ;;
        tempestremap) printf '%s' "$EXTRA_TEMPESTREMAP_ARGS" ;;
        *)            printf '' ;;
    esac
}

# Path to the per-TPL extras stamp. Lives inside the install prefix so it
# disappears with --clean-tpls and survives a bare `rm <install-marker>`.
tpl_extras_stamp() {
    printf '%s/.install-moab.extra_args' "$(tpl_install_dir "$1")"
}

#-----------------------------------------------------------------------------
# Build-with-retry driver
#-----------------------------------------------------------------------------
build_tpl_with_retries() {
    local tpl="$1"
    section "TPL: $tpl"

    if [[ "$REUSE_TPLS" == "yes" && -e "$(tpl_install_marker "$tpl")" ]]; then
        # If the user passed --extra-<tpl> and it differs from the args used to
        # produce the existing install, force a rebuild (otherwise we'd silently
        # reuse a stale install configured with the wrong BLAS/LAPACK args).
        local want_extras have_extras=""
        want_extras="$(tpl_extra_args "$tpl")"
        [[ -f "$(tpl_extras_stamp "$tpl")" ]] && have_extras="$(cat "$(tpl_extras_stamp "$tpl")")"
        if [[ "$want_extras" != "$have_extras" ]]; then
            warn "$tpl: --extra-$tpl changed since last install (have='$have_extras' want='$want_extras'); rebuilding"
        else
            ok "$tpl: already installed at $(tpl_install_dir "$tpl") (skipping; --no-reuse-tpls to force)"
            return 0
        fi
    fi

    reset_tpl_remediations

    local attempt=1
    while (( attempt <= TPL_MAX_ATTEMPTS )); do
        log "$tpl: attempt $attempt/$TPL_MAX_ATTEMPTS"

        prepare_tpl_source "$tpl" \
            || { warn "$tpl: source preparation failed"; (( attempt++ )); continue; }

        # Wipe build dir on retry to start clean
        if (( attempt > 1 )); then
            rm -rf "$(tpl_build_dir "$tpl")"
        fi

        local rc=0
        run_phase "$tpl" config  "${tpl}_configure" || rc=$?
        if (( rc == 0 )); then run_phase "$tpl" build   "${tpl}_build"   || rc=$?; fi
        if (( rc == 0 )); then run_phase "$tpl" install "${tpl}_install" || rc=$?; fi

        if (( rc == 0 )) && [[ -e "$(tpl_install_marker "$tpl")" ]]; then
            # Record the user-supplied extras so a future run can detect a change.
            mkdir -p "$(tpl_install_dir "$tpl")"
            printf '%s' "$(tpl_extra_args "$tpl")" > "$(tpl_extras_stamp "$tpl")"
            ok "$tpl: installed at $(tpl_install_dir "$tpl")"
            return 0
        fi

        warn "$tpl: attempt $attempt failed (rc=$rc)"
        # Inspect each phase log for known patterns
        local diag="" phase
        for phase in config build install; do
            diag="$(classify_log "$(tpl_log "$tpl" "$phase")")"
            [[ -n "$diag" ]] && { log "$tpl ($phase): pattern=$diag"; break; }
        done
        if [[ -z "$diag" ]]; then
            warn "$tpl: no known remediation pattern; aborting retries"
            return 1
        fi
        apply_remediation "$tpl" "$diag" || {
            warn "$tpl: remediation '$diag' not actionable; aborting"
            return 1
        }
        (( attempt++ ))
    done

    warn "$tpl: exhausted $TPL_MAX_ATTEMPTS attempts"
    return 1
}

#-----------------------------------------------------------------------------
# MOAB invocation (autotools or cmake)
#-----------------------------------------------------------------------------
toggle() {
    local flag="$1" val="$2"
    [[ "$val" == "yes" ]] && printf -- "--enable-%s" "$flag" || printf -- "--disable-%s" "$flag"
}

moab_autotools_args() {
    local with_mpi="--with-mpi"
    [[ -n "${MPI_ROOT:-}" ]] && with_mpi="--with-mpi=$MPI_ROOT"
    local args=(
        "$with_mpi"
        "CC=$CC_BIN" "CXX=$CXX_BIN" "FC=$FC_BIN" "F77=$F77_BIN"
        "--with-hdf5=$HDF5_ROOT"
        "--with-netcdf=$NETCDF_C_PATH"
        "--with-pnetcdf=$PNETCDF_PATH"
        "--with-eigen3=$(tpl_install_dir eigen3)/include/eigen3"
        "--with-zoltan=$(tpl_install_dir zoltan)"
        "--with-tempestremap=$(tpl_install_dir tempestremap)"
        "$(toggle optimize "$ENABLE_OPTIMIZE")"
        "$(toggle debug    "$ENABLE_DEBUG")"
        "$(toggle shared   "$ENABLE_SHARED")"
        "$(toggle static   "$ENABLE_STATIC")"
        "--prefix=$PREFIX_PATH"
    )
    if [[ -n "$EXTRA_MOAB_ARGS" ]]; then
        # eval honors embedded shell quoting, e.g. --extra='--with-blas="path -lgfortran"'
        local extra=()
        eval "extra=( $EXTRA_MOAB_ARGS )"
        args+=( ${extra[@]+"${extra[@]}"} )
    fi
    printf '%s\n' "${args[@]}"
}

moab_cmake_args() {
    local btype="Release"
    [[ "$ENABLE_DEBUG" == "yes" && "$ENABLE_OPTIMIZE" == "yes" ]] && btype="RelWithDebInfo"
    [[ "$ENABLE_DEBUG" == "yes" && "$ENABLE_OPTIMIZE" == "no"  ]] && btype="Debug"
    [[ "$ENABLE_DEBUG" == "no"  && "$ENABLE_OPTIMIZE" == "yes" ]] && btype="Release"
    local shared="OFF"; [[ "$ENABLE_SHARED" == "yes" ]] && shared="ON"
    local args=(
        "-DCMAKE_BUILD_TYPE=$btype"
        "-DCMAKE_INSTALL_PREFIX=$PREFIX_PATH"
        "-DBUILD_SHARED_LIBS=$shared"
        "-DCMAKE_C_COMPILER=$CC_BIN"
        "-DCMAKE_CXX_COMPILER=$CXX_BIN"
        "-DCMAKE_Fortran_COMPILER=$FC_BIN"
        "-DENABLE_MPI=ON"
    )
    [[ -n "${MPI_ROOT:-}" ]] && args+=( "-DMPI_HOME=$MPI_ROOT" )
    args+=(
        "-DENABLE_HDF5=ON"  "-DHDF5_ROOT=$HDF5_ROOT"
        "-DENABLE_NETCDF=ON" "-DNETCDF_DIR=$NETCDF_C_PATH"
        "-DENABLE_PNETCDF=ON" "-DPNETCDF_DIR=$PNETCDF_PATH"
        "-DENABLE_EIGEN3=ON"
        "-DEIGEN3_DIR=$(tpl_install_dir eigen3)/include/eigen3"
        "-DEIGEN3_INCLUDE_DIR=$(tpl_install_dir eigen3)/include/eigen3"
        "-DENABLE_ZOLTAN=ON"
        "-DZOLTAN_DIR=$(tpl_install_dir zoltan)"
        "-DENABLE_TEMPESTREMAP=ON"
        "-DTEMPESTREMAP_DIR=$(tpl_install_dir tempestremap)"
    )
    if [[ -n "$EXTRA_MOAB_ARGS" ]]; then
        local extra=()
        eval "extra=( $EXTRA_MOAB_ARGS )"
        args+=( ${extra[@]+"${extra[@]}"} )
    fi
    printf '%s\n' "${args[@]}"
}

print_moab_command() {
    case "$BUILD_SYSTEM" in
        autotools)
            printf '  %s/configure' "$MOAB_SRC_DIR"
            while IFS= read -r a; do printf ' \\\n      %q' "$a"; done < <(moab_autotools_args)
            ;;
        cmake)
            printf '  cmake -S %q -B %q' "$MOAB_SRC_DIR" "$MOAB_BUILD_DIR"
            while IFS= read -r a; do printf ' \\\n      %q' "$a"; done < <(moab_cmake_args)
            ;;
    esac
    printf '\n'
}

# Hash the planned configure args so we can detect when the user changed flags
# between runs and force a re-configure even with --reconfigure not set.
moab_args_fingerprint() {
    local payload
    case "$BUILD_SYSTEM" in
        autotools) payload="autotools|src=$(moab_src_hash)|$(moab_autotools_args)" ;;
        cmake)     payload="cmake|src=$(moab_src_hash)|$(moab_cmake_args)" ;;
    esac
    if   command -v sha256sum >/dev/null 2>&1; then
        printf '%s' "$payload" | sha256sum | awk '{print $1}'
    elif command -v shasum     >/dev/null 2>&1; then
        printf '%s' "$payload" | shasum -a 256 | awk '{print $1}'
    else
        printf '%s' "$payload" | cksum | awk '{print $1"-"$2}'
    fi
}

moab_already_configured() {
    case "$BUILD_SYSTEM" in
        autotools)
            [[ -f "$MOAB_BUILD_DIR/Makefile" && -f "$MOAB_BUILD_DIR/config.status" ]] || return 1
            ;;
        cmake)
            [[ -f "$MOAB_BUILD_DIR/CMakeCache.txt" ]] || return 1
            ;;
    esac
    local fp_file="$MOAB_BUILD_DIR/.install-moab.config_fingerprint"
    [[ -f "$fp_file" ]] || return 1
    [[ "$(cat "$fp_file" 2>/dev/null)" == "$(moab_args_fingerprint)" ]]
}

run_moab_configure() {
    section "MOAB configure ($BUILD_SYSTEM)"
    mkdir -p "$MOAB_BUILD_DIR"
    if [[ "$RECONFIGURE" != "yes" ]] && moab_already_configured; then
        ok "MOAB: already configured at $MOAB_BUILD_DIR (fingerprint match) -- skipping"
        log "    (use --reconfigure to force, or --clean to wipe BUILD_DIR)"
        return 0
    fi
    print_moab_command
    case "$BUILD_SYSTEM" in
        autotools)
            local args=()
            while IFS= read -r a; do args+=( "$a" ); done < <(moab_autotools_args)
            ( cd "$MOAB_BUILD_DIR" && "$MOAB_SRC_DIR/configure" "${args[@]}" ) \
                || { tail -n 60 "$MOAB_BUILD_DIR/config.log" >&2 || true; return 1; }
            ;;
        cmake)
            local args=()
            while IFS= read -r a; do args+=( "$a" ); done < <(moab_cmake_args)
            cmake -S "$MOAB_SRC_DIR" -B "$MOAB_BUILD_DIR" "${args[@]}" || return 1
            ;;
    esac
    moab_args_fingerprint > "$MOAB_BUILD_DIR/.install-moab.config_fingerprint"
}

run_moab_build_install() {
    section "MOAB build (-j$JOBS) and install"
    case "$BUILD_SYSTEM" in
        autotools)
            ( cd "$MOAB_BUILD_DIR" && make -j"$JOBS" ) || return 1
            [[ "$RUN_CHECK" == "yes" ]] && ( cd "$MOAB_BUILD_DIR" && make -j"$JOBS" check ) || true
            ( cd "$MOAB_BUILD_DIR" && make install )
            ;;
        cmake)
            cmake --build "$MOAB_BUILD_DIR" --parallel "$JOBS" || return 1
            if [[ "$RUN_CHECK" == "yes" ]]; then
                ( cd "$MOAB_BUILD_DIR" && ctest --output-on-failure -j"$JOBS" ) || true
            fi
            cmake --install "$MOAB_BUILD_DIR"
            ;;
    esac
}

#-----------------------------------------------------------------------------
# Verification
#-----------------------------------------------------------------------------
verify_install() {
    section "Verification"
    local ok_all=yes
    must_exist() { if [[ -e "$1" ]]; then ok "  found: $1"; else warn "  missing: $1"; ok_all=no; fi; }
    any_exist()  {
        local f
        for f in "$@"; do [[ -e "$f" ]] && { ok "  found: $f"; return 0; }; done
        warn "  missing all of: $*"; ok_all=no
    }
    log "TPLs:"
    must_exist "$(tpl_install_dir eigen3)/include/eigen3/Eigen/Dense"
    must_exist "$(tpl_install_dir zoltan)/include/zoltan.h"
    any_exist  "$(tpl_install_dir zoltan)/lib/libzoltan.a" \
               "$(tpl_install_dir zoltan)/lib/libzoltan.so" \
               "$(tpl_install_dir zoltan)/lib/libzoltan.dylib"
    must_exist "$(tpl_install_dir tempestremap)/include/TempestRemapAPI.h"
    any_exist  "$(tpl_install_dir tempestremap)/lib/libTempestRemap.a" \
               "$(tpl_install_dir tempestremap)/lib/libTempestRemap.so" \
               "$(tpl_install_dir tempestremap)/lib/libTempestRemap.dylib"
    log "MOAB:"
    must_exist "$PREFIX_PATH/include/moab/Core.hpp"
    any_exist  "$PREFIX_PATH/lib/libMOAB.a" \
               "$PREFIX_PATH/lib/libMOAB.so" \
               "$PREFIX_PATH/lib/libMOAB.dylib"
    [[ "$ok_all" == "yes" ]] && ok "verification PASSED" \
                             || warn "verification flagged missing artifacts"
}

#-----------------------------------------------------------------------------
# Report and dispatch
#-----------------------------------------------------------------------------
section "MOAB orchestration"
log "Profile         : $PROFILE${E3SM_ROOT:+ (E3SM_ROOT=$E3SM_ROOT)}"
if [[ -n "${MACHINE_NAME:-}" ]]; then
    log "Machine         : $MACHINE_NAME (compiler family: $COMPILER_FAMILY; last-validated: $MACHINE_META_LAST_VALIDATED${MACHINE_META_E3SM_NAME:+; e3sm_name=$MACHINE_META_E3SM_NAME})"
    [[ -n "$MACHINE_META_STANDALONE_HINT" ]] && log "Module hint     : $MACHINE_META_STANDALONE_HINT"
    [[ -n "$MACHINE_META_NOTES" ]] && log "Machine notes   : $MACHINE_META_NOTES"
else
    log "Machine         : (none detected; pass --machine=NAME or --list-machines to see options)"
fi
if [[ "${_USER_SET_MOAB_SRC:-no}" == "yes" ]]; then
    log "MOAB source dir : $MOAB_SRC_DIR (user-managed; no git ops)"
else
    log "MOAB source dir : $MOAB_SRC_DIR (script-managed: $MOAB_REPO_URL @ $MOAB_BRANCH)"
fi
log "Build dir       : $BUILD_DIR"
log "MOAB install    : $PREFIX_PATH"
log "TPL install     : $TPL_PREFIX"
log "Build system    : $BUILD_SYSTEM"
log "Parallel jobs   : $JOBS"
log "Compilers       : CC=$CC_BIN CXX=$CXX_BIN FC=$FC_BIN F77=$F77_BIN"
log "MPI_ROOT        : ${MPI_ROOT:-(unset; using PATH-based MPI wrappers)}"
log "HDF5_ROOT       : $HDF5_ROOT"
log "NETCDF_C_PATH   : $NETCDF_C_PATH"
log "PNETCDF_PATH    : $PNETCDF_PATH"
log "TPL mirror      : $TPL_MIRROR"
log "Downloader      : ${DOWNLOADER:-none-found}"
log "Reuse TPLs      : $REUSE_TPLS"
log "Tail TPL logs   : $TAIL_LOGS"
log "Per-TPL retries : $TPL_MAX_ATTEMPTS  (download retries: $DOWNLOAD_RETRIES)"

report_state() {
    section "Resume state"
    local tpl
    for tpl in eigen3 zoltan tempestremap; do
        if [[ -e "$(tpl_install_marker "$tpl")" ]]; then
            local want_extras have_extras=""
            want_extras="$(tpl_extra_args "$tpl")"
            [[ -f "$(tpl_extras_stamp "$tpl")" ]] && have_extras="$(cat "$(tpl_extras_stamp "$tpl")")"
            if [[ "$want_extras" != "$have_extras" ]]; then
                log "  $tpl: installed but --extra-$tpl changed (have='$have_extras' want='$want_extras') -- WILL REBUILD"
            else
                ok "  $tpl: installed at $(tpl_install_dir "$tpl") -- WILL SKIP"
            fi
        else
            log "  $tpl: not installed -- WILL BUILD"
        fi
    done
    if moab_already_configured; then
        if [[ "$RECONFIGURE" == "yes" ]]; then
            log "  moab: configured (fingerprint match) -- WILL RECONFIGURE (--reconfigure)"
        else
            ok "  moab: configured (fingerprint match) -- WILL SKIP CONFIGURE, resume make/install"
        fi
    elif [[ -f "$MOAB_BUILD_DIR/CMakeCache.txt" || -f "$MOAB_BUILD_DIR/Makefile" ]]; then
        log "  moab: configured but args differ (or fingerprint missing) -- WILL RECONFIGURE"
    else
        log "  moab: not configured -- WILL CONFIGURE"
    fi
}

if [[ "$PRINT_MODE" == "yes" ]]; then
    # Quiet copy-pasteable mode (absorbs the legacy suggest_configuration.sh
    # use case). Skip resume state, per-TPL recipes, and the orchestration
    # banner -- emit only the resolved configure / cmake command on stdout.
    print_moab_command
    exit 0
fi

if [[ "$DRY_RUN" == "yes" ]]; then
    mkdir -p "$BUILD_DIR" "$TPL_PREFIX" "$MOAB_BUILD_DIR" 2>/dev/null || true
    report_state
    log ""
    log "Per-TPL recipes (would run, in order):"
    for tpl in eigen3 zoltan tempestremap; do
        if [[ -n "$(tpl_git_repo "$tpl")" ]]; then
            src_desc="git --depth 1 $(tpl_git_repo "$tpl")@$(tpl_git_branch "$tpl")"
        else
            src_desc="$(tpl_url "$tpl")"
        fi
        log "  $tpl: $src_desc -> $(tpl_install_dir "$tpl")"
    done
    unset src_desc
    log ""
    log "MOAB invocation ($BUILD_SYSTEM):"
    print_moab_command
    exit 0
fi

if [[ "$CLEAN_BUILD" == "yes" && -d "$BUILD_DIR" ]]; then
    log "Removing existing build dir: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi
if [[ "$CLEAN_TPLS" == "yes" && -d "$TPL_PREFIX" ]]; then
    log "Removing existing TPL prefix: $TPL_PREFIX"
    rm -rf "$TPL_PREFIX"
fi
mkdir -p "$BUILD_DIR" "$ARCHIVES_DIR" "$TPL_WORK_DIR" "$TPL_LOG_DIR" "$TPL_PREFIX" "$MOAB_BUILD_DIR"

validate_all_mpi_wrappers
prepare_moab_source

report_state

# Build TPLs in order
for tpl in eigen3 zoltan tempestremap; do
    build_tpl_with_retries "$tpl" \
        || die "$tpl build failed; logs at $TPL_LOG_DIR/{config,build,install}_${tpl}.log" 3
done

# Then MOAB
run_moab_configure       || die "MOAB configure failed; see $MOAB_BUILD_DIR/config.log (autotools) or CMake output above" 3
run_moab_build_install   || die "MOAB build/install failed" 3

verify_install

section "Done"
ok "MOAB installed at: $PREFIX_PATH"
ok "TPLs installed under: $TPL_PREFIX"
log "Add the following to your environment for downstream use:"
log "  export MOAB_ROOT=$PREFIX_PATH"
log "  export PATH=\$MOAB_ROOT/bin:\$PATH"
log "  export PKG_CONFIG_PATH=\$MOAB_ROOT/lib/pkgconfig:\${PKG_CONFIG_PATH:-}"
log "  export ZOLTAN_DIR=$(tpl_install_dir zoltan)"
log "  export TEMPESTREMAP_DIR=$(tpl_install_dir tempestremap)"
log "  export EIGEN3_DIR=$(tpl_install_dir eigen3)/include/eigen3"
