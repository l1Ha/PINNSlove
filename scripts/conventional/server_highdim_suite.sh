#!/usr/bin/env bash
# Server-side orchestration for the high-dimensional (4D/6D) suites.
#
# Efficiency protocol (learned the hard way on this cluster):
#   * one process per independent workload, fanned out in parallel;
#   * OPENBLAS_NUM_THREADS = 1 for *sparse*-LU dominated jobs (SuperLU is
#     single-threaded; extra BLAS threads only cause oversubscription and
#     thrashing -- 12 jobs x 127 threads once turned a 3-minute solve into 8 h),
#     small values (2-4) only for *dense* linear-algebra dominated jobs;
#   * absolute paths everywhere (a `cd` inside a `&` group does not apply to
#     the sibling job);
#   * every job records wall time and peak RSS.
#
# Usage:  bash scripts/server_highdim_suite.sh
set -u
ROOT=/storage/home/lih/PINNSlove
PY=$HOME/aqd-venv/bin/python
LOGS=$ROOT/results/highdim_logs
mkdir -p "$LOGS"

run_job() {                    # name, threads, logfile, command...
    local name=$1 threads=$2 log=$3; shift 3
    local t0=$(date +%s.%N)
    env OPENBLAS_NUM_THREADS=$threads OMP_NUM_THREADS=$threads MKL_NUM_THREADS=$threads \
        "$@" > "$log" 2>&1
    local rc=$?
    local t1=$(date +%s.%N)
    local wall=$(awk "BEGIN{printf \"%.2f\", $t1 - $t0}")
    echo "$name rc=$rc wall=${wall}s" > "$log.time"
}

# 1) high-dimensional bound states (4D H2CO + 6D tetratomic): dense eigh -> 4 threads
run_job "highdim_bound" 4 "$LOGS/bound.log" \
    $PY "$ROOT/scripts/run_high_dim_bound.py" &
PID_BOUND=$!

# 2) 4D diatom-diatom scattering production: tiny dense solves -> 2 threads
run_job "scat4d_prod" 2 "$LOGS/scat4d.log" \
    $PY "$ROOT/scripts/run_diatom_diatom_4d.py" --production &
PID_SCAT=$!

wait $PID_BOUND $PID_SCAT

echo "=== timing ==="
cat "$LOGS"/*.time
echo "=== status ==="
tail -4 "$LOGS/bound.log"
tail -5 "$LOGS/scat4d.log"
