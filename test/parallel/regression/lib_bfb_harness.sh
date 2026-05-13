#!/bin/bash
# lib_bfb_harness.sh - strict byte-for-byte regression driver.
#
# For tests that can (a) write computed maps to disk via $ONLINE_FLAGS, and
# (b) emit per-cell digests of projected fields via $DIGEST_FLAG. Algorithm:
#
#   Phase A: serial run, --compute_online --write_maps  -> baseline_*.nc
#   Phase B: serial run, reload + emit digests           -> ${DIGEST_PREFIX}_${k}_1.txt
#   Phase C: for each n in RANKS (skip 1): same reload   -> ${DIGEST_PREFIX}_${k}_${n}.txt
#   Phase D: cmp -s every (kernel, n) digest pair vs the n=1 baseline
#   Phase E (optional): if config sets $H5M_ARTIFACTS, also strict-compare
#            those h5m outputs across ranks via mbcmpfiles (tolerance 0).
#
# Required from the sourced config:
#   EXE, TEST_NAME, SRC_FLAG, TGT_FLAG, DEFAULT_SRC, DEFAULT_TGT,
#   ONLINE_FLAGS, LOAD_FLAGS, DIGEST_FLAG, DIGEST_PREFIX, DIGEST_KERNELS,
#   DEFAULT_RANKS
# Optional:
#   EXTRA_ARGS, H5M_ARTIFACTS

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh"

run_bfb_harness() {
    local src="$1" tgt="$2"; shift 2
    local ranks=("$@")
    [ ${#ranks[@]} -gt 0 ] || ranks=( $DEFAULT_RANKS )

    setup_workdir
    cd "$RUN_WORKDIR" || die "cannot cd $RUN_WORKDIR"

    log "=== BfB regression: $TEST_NAME ==="
    log " EXE     : $EXE"
    log " SRC     : $src"
    log " TGT     : $tgt"
    log " RANKS   : ${ranks[*]}"
    log " WORKDIR : $RUN_WORKDIR"

    section "Phase A: serial compute online + write maps"
    run_logged step_A.log $MPIRUN -n 1 "$EXE" \
        $SRC_FLAG "$src" $TGT_FLAG "$tgt" \
        $ONLINE_FLAGS ${EXTRA_ARGS:-} \
        || die "serial map-write step exited non-zero"

    _stash_artifacts() {
        local n="$1"
        # Namespace the renamed copy by $TEST_NAME so concurrent runs of
        # different configs (which may write the same hardcoded basename)
        # cannot overwrite each other's reference artifacts in a shared
        # $WORKDIR.
        for art in ${H5M_ARTIFACTS:-}; do
            [ -s "$art" ] || die "expected artifact '$art' missing/empty after n=$n"
            mv "$art" "${TEST_NAME}_${art}.n${n}"
        done
    }

    section "Phase B: serial reload + baseline digests (n=1)"
    run_logged step_B.log $MPIRUN -n 1 "$EXE" \
        $SRC_FLAG "$src" $TGT_FLAG "$tgt" \
        $LOAD_FLAGS $DIGEST_FLAG ${EXTRA_ARGS:-} \
        || die "serial reload step exited non-zero"
    for k in $DIGEST_KERNELS; do
        local f="${DIGEST_PREFIX}_${k}_1.txt"
        [ -s "$f" ] || die "expected serial digest '$f' missing/empty"
    done
    _stash_artifacts 1

    section "Phase C: parallel reload + digests"
    for n in "${ranks[@]}"; do
        [ "$n" = "1" ] && continue
        run_logged "step_C_n${n}.log" $MPIRUN -n "$n" "$EXE" \
            $SRC_FLAG "$src" $TGT_FLAG "$tgt" \
            $LOAD_FLAGS $DIGEST_FLAG ${EXTRA_ARGS:-} \
            || die "parallel run at n=$n exited non-zero"
        _stash_artifacts "$n"
    done

    local dtol="${DIGEST_TOLERANCE:-}"
    if [ -n "$dtol" ]; then
        section "Phase D: tolerance digest compare (rel-tol=$dtol)"
    else
        section "Phase D: strict cmp of per-cell digests"
    fi
    local fail=0
    for k in $DIGEST_KERNELS; do
        for n in "${ranks[@]}"; do
            [ "$n" = "1" ] && continue
            if [ -n "$dtol" ]; then
                tol_digest_cmp "${DIGEST_PREFIX}_${k}_1.txt" \
                               "${DIGEST_PREFIX}_${k}_${n}.txt" \
                               "$dtol" "${k} n=${n}" \
                    || fail=$((fail+1))
            else
                strict_cmp "${DIGEST_PREFIX}_${k}_1.txt" \
                           "${DIGEST_PREFIX}_${k}_${n}.txt" \
                           "${k} n=${n}" \
                    || fail=$((fail+1))
            fi
        done
    done

    if [ -n "${H5M_ARTIFACTS:-}" ]; then
        section "Phase E: mbcmpfiles strict h5m compare (serial vs parallel)"
        for art in $H5M_ARTIFACTS; do
            for n in "${ranks[@]}"; do
                [ "$n" = "1" ] && continue
                mbcmp_h5m "${TEST_NAME}_${art}.n1" "${TEST_NAME}_${art}.n${n}" 0 "${art} n=${n}" \
                    || fail=$((fail+1))
            done
        done
    fi

    log ""
    if [ $fail -eq 0 ]; then
        if [ -n "$dtol" ]; then
            log "OK: $TEST_NAME is near-BfB (rel-tol=$dtol) across rank counts ${ranks[*]}"
        else
            log "OK: $TEST_NAME is BfB across rank counts ${ranks[*]}"
        fi
        return 0
    fi
    log "FAIL: $TEST_NAME has $fail BfB regression(s)"
    return 1
}
