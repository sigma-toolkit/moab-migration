#!/bin/bash
# lib_exitcode_harness.sh - exit-code regression driver.
#
# For tests that self-verify internally (check_baseline_file, CHECKIERR)
# and signal pass/fail via exit code. The harness simply runs the test at
# each requested rank count and asserts exit code 0.
#
# Required from the sourced config:
#   EXE, TEST_NAME, DEFAULT_RANKS
# Optional:
#   SRC_FLAG, TGT_FLAG, DEFAULT_SRC, DEFAULT_TGT, EXTRA_ARGS

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh"

run_exitcode_harness() {
    local ranks=("$@")
    [ ${#ranks[@]} -gt 0 ] || ranks=( $DEFAULT_RANKS )

    local src="${DEFAULT_SRC:-}"
    local tgt="${DEFAULT_TGT:-}"

    setup_workdir
    cd "$RUN_WORKDIR" || die "cannot cd $RUN_WORKDIR"

    log "=== exit-code regression: $TEST_NAME ==="
    log " EXE    : $EXE"
    log " RANKS  : ${ranks[*]}"
    log " WORKDIR: $RUN_WORKDIR"

    local src_args=""
    [ -n "${SRC_FLAG:-}" ] && [ -n "$src" ] && src_args="$SRC_FLAG $src"
    local tgt_args=""
    [ -n "${TGT_FLAG:-}" ] && [ -n "$tgt" ] && tgt_args="$TGT_FLAG $tgt"

    local fail=0
    for n in "${ranks[@]}"; do
        local logfile="run_n${n}.log"
        section "Run at n=$n"
        if run_logged "$logfile" $MPIRUN -n "$n" "$EXE" \
                $src_args $tgt_args ${EXTRA_ARGS:-}; then
            printf '  PASS  n=%s: exit 0\n' "$n"
        else
            printf '  FAIL  n=%s: non-zero exit (see %s)\n' "$n" "$logfile"
            fail=$((fail+1))
        fi
    done

    log ""
    if [ $fail -eq 0 ]; then
        log "OK: $TEST_NAME passed at all rank counts ${ranks[*]}"
        return 0
    fi
    log "FAIL: $TEST_NAME failed at $fail rank count(s)"
    return 1
}
