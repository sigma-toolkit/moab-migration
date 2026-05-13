#!/bin/bash
# lib_tol_harness.sh - tolerance-based regression driver.
#
# For tests whose online-computed weights are NOT BfB across rank counts
# (every iMOAB coupler test). Algorithm:
#
#   Phase A: serial run; rename each artifact in $H5M_ARTIFACTS to <name>.n1
#   Phase B: for each n in RANKS (skip 1): run, rename artifacts to <name>.n${n}
#   Phase C: mbcmpfiles each <name>.n1 vs <name>.n${n}, asserting all common
#            double-tag L2-diffs <= TOLERANCE (default 1e-9).
#
# Required from the sourced config:
#   EXE, SRC_FLAG, TGT_FLAG, DEFAULT_SRC, DEFAULT_TGT,
#   H5M_ARTIFACTS  (space-separated artifact basenames the test writes to CWD),
#   DEFAULT_RANKS, TOLERANCE
# Optional:
#   EXTRA_ARGS, EXTRA_FLAG_1/EXTRA_DEFAULT_1 (e.g. -l for land mesh)

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh"

run_tol_harness() {
    local src="$1" tgt="$2"; shift 2
    local ranks=("$@")
    [ ${#ranks[@]} -gt 0 ] || ranks=( $DEFAULT_RANKS )
    local tol="${TOLERANCE:-1e-9}"

    setup_workdir
    cd "$RUN_WORKDIR" || die "cannot cd $RUN_WORKDIR"

    log "=== tolerance regression: $TEST_NAME ==="
    log " EXE       : $EXE"
    log " SRC       : $src"
    log " TGT       : $tgt"
    log " RANKS     : ${ranks[*]}"
    log " TOLERANCE : $tol"
    log " ARTIFACTS : $H5M_ARTIFACTS"
    log " WORKDIR   : $RUN_WORKDIR"

    local extra=""
    [ -n "${EXTRA_FLAG_1:-}" ] && extra="$EXTRA_FLAG_1 ${EXTRA_VAL_1:-$EXTRA_DEFAULT_1}"

    _run_and_capture() {
        local n="$1"
        local logfile="step_n${n}.log"
        run_logged "$logfile" $MPIRUN -n "$n" "$EXE" \
            $SRC_FLAG "$src" $TGT_FLAG "$tgt" $extra ${EXTRA_ARGS:-} \
            || die "run at n=$n exited non-zero"
        # Namespace the renamed copy by $TEST_NAME so concurrent configs
        # writing the same hardcoded artifact name (e.g. recvAtm.h5m) into
        # a shared $WORKDIR do not clobber each other's references.
        for art in $H5M_ARTIFACTS; do
            [ -s "$art" ] || die "expected artifact '$art' missing/empty after n=$n"
            mv "$art" "${TEST_NAME}_${art}.n${n}"
        done
    }

    section "Phase A: serial reference run (n=1)"
    _run_and_capture 1

    section "Phase B: parallel runs"
    for n in "${ranks[@]}"; do
        [ "$n" = "1" ] && continue
        _run_and_capture "$n"
    done

    section "Phase C: mbcmpfiles tolerance compare vs serial reference"
    local fail=0
    for art in $H5M_ARTIFACTS; do
        for n in "${ranks[@]}"; do
            [ "$n" = "1" ] && continue
            mbcmp_h5m "${TEST_NAME}_${art}.n1" "${TEST_NAME}_${art}.n${n}" "$tol" "${art} n=${n}" \
                || fail=$((fail+1))
        done
    done

    log ""
    if [ $fail -eq 0 ]; then
        log "OK: $TEST_NAME within tolerance $tol across rank counts ${ranks[*]}"
        return 0
    fi
    log "FAIL: $TEST_NAME has $fail tolerance regression(s)"
    return 1
}
