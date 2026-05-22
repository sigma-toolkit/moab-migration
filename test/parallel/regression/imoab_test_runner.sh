#!/bin/bash
# imoab_test_runner.sh - top-level dispatcher for the iMOAB parallel-test
# regression harness.
#
# Usage:
#   imoab_test_runner.sh <config-name> [options]
#
# <config-name> is the basename (without .cfg) of a file under configs/.
# Run with --list to see available configs.
#
# Options:
#   --src      <path>       Source mesh file (overrides config default).
#   --tgt      <path>       Target mesh file (overrides config default).
#   --ranks    "N1 N2 ..."  Rank counts to sweep (overrides DEFAULT_RANKS).
#                           Must include 1 (used as the serial baseline).
#   --workdir  <dir>        Use <dir> for intermediate files (default: mktemp,
#                           auto-cleaned). Implies --keep.
#   --keep                  Preserve workdir even if it was mktemp'd.
#   --mbcmp    <path>       Path to mbcmpfiles binary (overrides $MBCMPFILES).
#   --tol      <value>      Tolerance for tol-mode configs (overrides config).
#   --extra    <flag=value> Pass-through to EXTRA_FLAG_1 (e.g. --extra l=/path/to/lnd.h5m).
#   -h|--help               Show this message.

set -u

HERE=$(cd "$(dirname "$0")" && pwd)
CONFIG_DIR="$HERE/configs"

usage() { sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//' >&2; }
list_configs() {
    log_dir=$(ls "$CONFIG_DIR"/*.cfg 2>/dev/null)
    if [ -z "$log_dir" ]; then
        echo "(no configs found in $CONFIG_DIR)" >&2
        return
    fi
    echo "Available configs:"
    for f in "$CONFIG_DIR"/*.cfg; do
        local name="${f##*/}"; name="${name%.cfg}"
        local desc
        desc=$(awk '/^# CONFIG_DESC:/{sub(/^# CONFIG_DESC: */,""); print; exit}' "$f")
        printf '  %-24s %s\n' "$name" "$desc"
    done
}

[ $# -ge 1 ] || { usage; exit 2; }
case "$1" in
    -h|--help) usage; exit 0 ;;
    --list)    list_configs; exit 0 ;;
esac

CONFIG_NAME="$1"; shift
CFG="$CONFIG_DIR/$CONFIG_NAME.cfg"
[ -r "$CFG" ] || { echo "ERROR: config not found: $CFG" >&2; list_configs; exit 2; }

SRC_OVERRIDE=""
TGT_OVERRIDE=""
RANKS_OVERRIDE=""
EXTRA_OVERRIDE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --src)     SRC_OVERRIDE="$2"; shift 2 ;;
        --tgt)     TGT_OVERRIDE="$2"; shift 2 ;;
        --ranks)   RANKS_OVERRIDE="$2"; shift 2 ;;
        --workdir) export WORKDIR="$2"; export KEEP_WORKDIR=1; shift 2 ;;
        --keep)    export KEEP_WORKDIR=1; shift ;;
        --mbcmp)   export MBCMPFILES="$2"; shift 2 ;;
        --tol)     export TOLERANCE="$2"; shift 2 ;;
        --extra)   EXTRA_OVERRIDE="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)         echo "ERROR: unknown option: $1" >&2; usage; exit 2 ;;
    esac
done

source "$CFG"
: "${EXE:?config must set EXE}"
: "${TEST_NAME:?config must set TEST_NAME}"
: "${TEST_MODE:?config must set TEST_MODE (bfb|tol|exitcode)}"
: "${DEFAULT_RANKS:?config must set DEFAULT_RANKS}"

# When IMOAB_BINDIR is set (e.g. by the CMake build, which puts binaries in
# ${CMAKE_BINARY_DIR}/bin rather than the autotools-style
# ${MOAB_BUILD}/test/parallel), rewrite EXE to live under it. Configs only
# need to encode the basename correctly; the binary's location varies by
# build system.
if [ -n "${IMOAB_BINDIR:-}" ]; then
    EXE="$IMOAB_BINDIR/$(basename "$EXE")"
fi

[ -x "$EXE" ] || { echo "ERROR: EXE not executable: $EXE" >&2; exit 2; }

if [ "${TEST_MODE}" != "exitcode" ]; then
    : "${SRC_FLAG:?config must set SRC_FLAG}"
    : "${TGT_FLAG:?config must set TGT_FLAG}"
    : "${DEFAULT_SRC:?config must set DEFAULT_SRC}"
    : "${DEFAULT_TGT:?config must set DEFAULT_TGT}"
    SRC="${SRC_OVERRIDE:-$DEFAULT_SRC}"
    TGT="${TGT_OVERRIDE:-$DEFAULT_TGT}"
    [ -r "$SRC" ] || { echo "ERROR: source mesh not readable: $SRC" >&2; exit 2; }
    [ -r "$TGT" ] || { echo "ERROR: target mesh not readable: $TGT" >&2; exit 2; }
fi

if [ -n "$EXTRA_OVERRIDE" ]; then
    EXTRA_VAL_1="${EXTRA_OVERRIDE#*=}"
    export EXTRA_VAL_1
fi

if [ -n "$RANKS_OVERRIDE" ]; then
    RANKS_ARGS=( $RANKS_OVERRIDE )
else
    RANKS_ARGS=()
fi

# ${RANKS_ARGS[@]+"${RANKS_ARGS[@]}"} is the portable empty-array expansion:
# bash 3.2 (the system bash on macOS) errors on "${ARR[@]}" when ARR=() under
# `set -u`; the +alt form expands to nothing when unset and the array otherwise.
case "$TEST_MODE" in
    bfb)      source "$HERE/lib_bfb_harness.sh";      run_bfb_harness      "$SRC" "$TGT" ${RANKS_ARGS[@]+"${RANKS_ARGS[@]}"} ;;
    tol)      source "$HERE/lib_tol_harness.sh";      run_tol_harness      "$SRC" "$TGT" ${RANKS_ARGS[@]+"${RANKS_ARGS[@]}"} ;;
    exitcode) source "$HERE/lib_exitcode_harness.sh"; run_exitcode_harness              ${RANKS_ARGS[@]+"${RANKS_ARGS[@]}"} ;;
    *)   echo "ERROR: unknown TEST_MODE '$TEST_MODE' (expected bfb|tol|exitcode)" >&2; exit 2 ;;
esac
