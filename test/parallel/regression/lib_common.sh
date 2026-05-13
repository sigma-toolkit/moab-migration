#!/bin/bash
# lib_common.sh - shared helpers for the iMOAB parallel-test regression harness.
# Sourced by lib_bfb_harness.sh and lib_tol_harness.sh; not meant to be run
# standalone.

set -u

log() { printf '%s\n' "$*"; }
die() { printf 'FAIL: %s\n' "$*" >&2; exit 1; }
section() { printf '\n--- %s ---\n' "$*"; }

# Per-config defaults the harness can override but most configs leave alone.
: "${MPIRUN:=mpirun}"
: "${MBCMPFILES:=}"      # path to mbcmpfiles binary; empty disables h5m compare
: "${WORKDIR:=}"         # if empty, mktemp; else use as-is and don't auto-rm
: "${KEEP_WORKDIR:=0}"   # set to 1 (or pass --keep) to preserve intermediate files

# Run a command, tee stdout+stderr to a log file, fail the harness on non-zero.
# Args: <log-file> <cmd> [args...]
run_logged() {
    local logfile="$1"; shift
    log "  + $*"
    if ! "$@" >"$logfile" 2>&1; then
        printf 'FAIL: command exited non-zero (see %s)\n' "$logfile" >&2
        sed 's/^/    /' "$logfile" | tail -40 >&2
        return 1
    fi
    return 0
}

# Compare two text files byte-for-byte.
# Args: <reference> <candidate> <label>
strict_cmp() {
    local ref="$1" cand="$2" label="$3"
    if [ ! -s "$ref" ];  then printf '  FAIL  %s: reference %s missing/empty\n' "$label" "$ref"; return 1; fi
    if [ ! -s "$cand" ]; then printf '  FAIL  %s: candidate %s missing/empty\n' "$label" "$cand"; return 1; fi
    if cmp -s "$ref" "$cand"; then
        printf '  PASS  %s: %s == %s\n' "$label" "$cand" "$ref"
        return 0
    fi
    printf '  FAIL  %s: %s != %s\n' "$label" "$cand" "$ref"
    local ndiff
    ndiff=$(diff "$ref" "$cand" 2>/dev/null | grep -c '^<')
    printf '          (%s/%s lines differ; sample:)\n' "$ndiff" "$(wc -l <"$ref")"
    diff "$ref" "$cand" | head -5 | sed 's/^/            /'
    return 1
}

# Compare two h5m files via mbcmpfiles, asserting L2-norm of every common
# double tag is <= tolerance. tolerance==0 enforces strict BfB.
# Args: <ref.h5m> <cand.h5m> <tolerance> <label>
mbcmp_h5m() {
    local ref="$1" cand="$2" tol="$3" label="$4"
    if [ -z "$MBCMPFILES" ] || [ ! -x "$MBCMPFILES" ]; then
        printf '  SKIP  %s: mbcmpfiles not configured (set MBCMPFILES env var)\n' "$label"
        return 0
    fi
    [ -s "$ref" ]  || { printf '  FAIL  %s: reference %s missing\n' "$label" "$ref"; return 1; }
    [ -s "$cand" ] || { printf '  FAIL  %s: candidate %s missing\n' "$label" "$cand"; return 1; }
    local out
    out=$("$MBCMPFILES" -i "$ref" -j "$cand" 2>&1) || {
        printf '  FAIL  %s: mbcmpfiles errored\n' "$label"; printf '%s\n' "$out" | sed 's/^/    /' | tail -20
        return 1
    }
    # Parse every "l2norm of the diff: <value>" line; fail on the first one
    # that exceeds tolerance.
    local bad=0
    while read -r norm; do
        # awk floating compare; treats tol=="0" as strict (any nonzero norm fails)
        if awk -v n="$norm" -v t="$tol" 'BEGIN{exit !(n+0 > t+0)}'; then
            printf '  FAIL  %s: tag L2-diff = %s (tolerance %s)\n' "$label" "$norm" "$tol"
            bad=1
        fi
    done < <(printf '%s\n' "$out" | awk -F': ' '/l2norm of the diff/{print $NF}')
    if [ $bad -eq 0 ]; then
        printf '  PASS  %s: all tag L2-diffs <= %s\n' "$label" "$tol"
        return 0
    fi
    printf '%s\n' "$out" | sed 's/^/    /' | tail -20
    return 1
}

# Compare two digest text files with a relative tolerance.
# Digest format: "<global_id> <value>" per line.
# Args: <reference> <candidate> <tolerance> <label>
tol_digest_cmp() {
    local ref="$1" cand="$2" tol="$3" label="$4"
    if [ ! -s "$ref" ];  then printf '  FAIL  %s: reference %s missing/empty\n' "$label" "$ref"; return 1; fi
    if [ ! -s "$cand" ]; then printf '  FAIL  %s: candidate %s missing/empty\n' "$label" "$cand"; return 1; fi
    local result
    result=$(awk -v tol="$tol" '
        NR==FNR { ref[$1]=$2; next }
        {
            id=$1; val=$2
            if (!(id in ref)) { bad++; next }
            rv = ref[id]+0; cv = val+0
            denom = (rv == 0 ? 1.0 : (rv < 0 ? -rv : rv))
            rd = ((cv-rv) < 0 ? -(cv-rv) : (cv-rv)) / denom
            if (rd > tol+0) { bad++; if (bad <= 3) printf "    id=%s ref=%s cand=%s reldiff=%.3e\n", id, ref[id], val, rd }
        }
        END { print bad+0 }
    ' "$ref" "$cand")
    local nbad
    nbad=$(printf '%s\n' "$result" | tail -1)
    if [ "$nbad" -eq 0 ]; then
        printf '  PASS  %s: all values within rel-tol %s\n' "$label" "$tol"
        return 0
    fi
    printf '  FAIL  %s: %s values exceed rel-tol %s\n' "$label" "$nbad" "$tol"
    printf '%s\n' "$result" | head -3
    return 1
}

# Resolve the workdir; honours $WORKDIR (preserve) vs mktemp (auto-clean).
# Sets: $RUN_WORKDIR (absolute path); on auto-clean, registers EXIT trap.
setup_workdir() {
    if [ -n "$WORKDIR" ]; then
        mkdir -p "$WORKDIR"
        RUN_WORKDIR=$(cd "$WORKDIR" && pwd)
    else
        RUN_WORKDIR=$(mktemp -d -t imoab_regr_XXXXXX)
        if [ "$KEEP_WORKDIR" -eq 0 ]; then
            trap 'rm -rf "$RUN_WORKDIR"' EXIT
        fi
    fi
    export RUN_WORKDIR
}
