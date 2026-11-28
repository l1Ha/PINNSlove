#!/usr/bin/env bash
# =============================================================================
# 项目服务器流水线：代码同步 → 服务器实算（回归 + 全部基准）→ 结果回收。
#
# 用法:  bash scripts/server/run_pipeline.sh [SERVER]
#   SERVER 默认 192.168.6.24 (xb002, 192 核 + A100-40GB, venv ~/venvs/pinn)。
#
# 服务器端执行内容（全部实算，无任何硬编码结果，符合 AIPro.md）:
#   1. 8 级回归套件            scripts/validation/run_full_regression.py
#   2. 1D 散射三方法对比        scripts/pinn/run_scatter_1d.py
#   3. 2D CAP 基准             scripts/pinn/run_scatter_2d_cap_benchmark.py
#   4. 2D ELM 基准             scripts/pinn/run_scatter_2d_elm.py
#   5. 2D 载波基准             scripts/pinn/run_scatter_2d_carrier.py
#   6. 2D diabatic 基准        scripts/pinn/run_scatter_2d_diabatic.py
#   7. 收敛研究 / CAP 扫描      scripts/conventional/run_scatter2d_convergence.py
#                              scripts/conventional/run_cap_scan.py
#
# 结果回收: results/server_2d/（json + 日志），本地 python 结果保留在 results/。
# =============================================================================
set -euo pipefail

SERVER="${1:-192.168.6.24}"
LOCAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REMOTE_DIR='~/PINNSlove'
PY="$REMOTE_DIR/../../venvs/pinn/bin/python"   # ~/PINNSlove/../../venvs = ~/venvs
PY="~/venvs/pinn/bin/python"

echo "==> [1/4] 同步代码到 $SERVER:$REMOTE_DIR"
cd "$LOCAL_DIR"
tar czf - --exclude='__pycache__' --exclude='results' --exclude='.git' \
    --exclude='*.pyc' --exclude='peslib' --exclude='ref' --exclude='article' \
    --exclude='figures' --exclude='.zcode' . \
  | ssh "$SERVER" "mkdir -p $REMOTE_DIR && tar xzf - -C $REMOTE_DIR"

echo "==> [2/4] 服务器实算（顺序执行，日志逐一存档）"
ssh "$SERVER" "cd $REMOTE_DIR && mkdir -p results && bash -c '
set -x
$PY scripts/validation/run_full_regression.py      > results/regression.log        2>&1
$PY scripts/pinn/run_scatter_1d.py                 > results/scatter_1d.log        2>&1
$PY scripts/conventional/run_scatter2d_convergence.py > results/convergence.log    2>&1
$PY scripts/conventional/run_cap_scan.py           > results/cap_scan.log          2>&1
$PY scripts/pinn/run_scatter_2d_cap_benchmark.py   > results/cap_benchmark.log     2>&1
$PY scripts/pinn/run_scatter_2d_elm.py             > results/elm_benchmark.log     2>&1
$PY scripts/pinn/run_scatter_2d_carrier.py         > results/carrier_benchmark.log 2>&1
$PY scripts/pinn/run_scatter_2d_diabatic.py        > results/diabatic_benchmark.log 2>&1
echo ALL-DONE
'"

echo "==> [3/4] 回收结果到 results/server_2d/"
mkdir -p "$LOCAL_DIR/results/server_2d"
scp -q "$SERVER:$REMOTE_DIR/results/*.log" "$LOCAL_DIR/results/server_2d/" || true
scp -q "$SERVER:$REMOTE_DIR/results/*.json" "$LOCAL_DIR/results/server_2d/" || true

echo "==> [4/4] 摘要"
grep -E "PASS|FAIL|COMPLETED" "$LOCAL_DIR/results/server_2d/regression.log" | tail -12 || true
tail -3 "$LOCAL_DIR/results/server_2d/diabatic_benchmark.log" || true
echo "流水线完成。"
