#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════
# ABC_S 一体化管线: 构建 -> 生产运行(两波,断点续算) -> 提取速率 -> 对比图
#
# 用法:
#   bash run_pipeline.sh [SCRATCH_DIR] [THREADS_PER_RUN] [MODE]
#     SCRATCH_DIR  运行工作目录 (默认 <包>/runs; 建议放本地盘/高速盘)
#     THREADS      每套运行的 OpenBLAS/OMP 线程数 (默认 4)
#     MODE         all (默认) = NGP波+GP波+分析+出图
#                  ngp       = 只跑 NGP 波 + 分析 + 出图
#                  analyze   = 只重跑分析 + 出图 (不启动任何运行)
#
# 断点续算: 已达到 nnrg 行的运行自动跳过; 重跑本脚本即可续。
# 运行完后 figures/ 下生成三张对比图 (vs Kendrick 2018 / 2016 / GP-NGP 比值)。
# ═══════════════════════════════════════════════════════════════════════════
set -u
PKG="$(cd "$(dirname "$0")" && pwd)"
SCRATCH="${1:-$PKG/runs}"
TH="${2:-4}"
MODE="${3:-all}"
mkdir -p "$SCRATCH" "$PKG/figures" "$PKG/analysis"

# ── 1. 编译器与 BLAS: 服务器 ifort+MKL, 兜底 gfortran+OpenBLAS ──────────────
# Intel oneAPI 环境 (ifort 运行时库 + MKL 库路径):
for sv in /opt/intel/oneapi/setvars.sh "$HOME/intel/oneapi/setvars.sh"; do
  [ -f "$sv" ] && { source "$sv" > /dev/null 2>&1; break; }
done
if command -v ifort > /dev/null 2>&1; then
  FPM_FC="${FPM_FC:-ifort}"
  FPM_LDFLAGS="${FPM_LDFLAGS:-}"
  if [ -n "${MKLROOT:-}" ] || ls /opt/intel/oneapi/mkl >/dev/null 2>&1 || ls "$HOME/intel/oneapi/mkl" >/dev/null 2>&1; then
    # 激活 fpm.toml 中的 MKL 链接行 (幂等), blas 行注释掉
    sed -i 's|^link=\["blas","lapack"\]|#link=["blas","lapack"]|; s|^#link=\["mkl_intel_lp64","mkl_intel_thread","mkl_core","pthread","m","dl"\]|link=["mkl_intel_lp64","mkl_intel_thread","mkl_core","pthread","m","dl"]|' "$PKG/fpm.toml"
    FPM_LDFLAGS="-L${MKLROOT:-/opt/intel/oneapi/mkl}/lib/intel64 -L${MKLROOT:-/opt/intel/oneapi/mkl}/lib $FPM_LDFLAGS"
    echo "[$(date)] compiler: ifort + MKL (threaded)"
  else
    echo "[$(date)] compiler: ifort (MKL not found -> system blas/lapack)"
  fi
  BIN=""    # 包内预编译二进制为 gfortran 构建 -> ifort 环境强制重新构建
else
  FPM_FC="${FPM_FC:-gfortran}"
  FPM_LDFLAGS="${FPM_LDFLAGS:-}"
  # 确保 gfortran 用 blas/lapack 行 (幂等), MKL 行注释掉
  sed -i 's|^#link=\["blas","lapack"\]|link=["blas","lapack"]|; s|^link=\["mkl_intel_lp64","mkl_intel_thread","mkl_core","pthread","m","dl"\]|#link=["mkl_intel_lp64","mkl_intel_thread","mkl_core","pthread","m","dl"]|' "$PKG/fpm.toml"
  BIN=$(ls -t "$PKG"/build/*/app/abc 2>/dev/null | head -1)
  if [ -n "$BIN" ]; then
    "$BIN" < /dev/null > /dev/null 2>&1
    rc=$?
    if [ "$rc" -eq 126 ] || [ "$rc" -eq 127 ]; then
      echo "[$(date)] prebuilt binary not executable on this host (rc=$rc) -> rebuild"
      BIN=""
    else
      echo "[$(date)] using prebuilt binary: $BIN"
    fi
  fi
fi
if [ -z "$BIN" ]; then
  echo "[$(date)] fpm build --profile release ($FPM_FC, $FPM_LDFLAGS) ..."
  if ! (cd "$PKG" && FPM_FC="$FPM_FC" FPM_LDFLAGS="$FPM_LDFLAGS" fpm build --profile release); then
    echo "ERROR: fpm build failed (need fpm + $FPM_FC + BLAS/LAPACK)"; exit 1
  fi
  BIN=$(ls -t "$PKG"/build/*/app/abc 2>/dev/null | head -1)
fi
echo "[$(date)] binary: $BIN"

# ── 2. 运行组 ──────────────────────────────────────────────────────────────
WAVE_NGP="h3_J0_even_ngp h3_J0_odd_ngp h3_J1_ppos_ngp h3_J1_pneg_ngp h3_J2_ppos_ngp h3_J2_pneg_ngp dhd_ngp h3_J0_even_ngp_muk h3_J0_odd_ngp_muk dhd_ngp_muk"
WAVE_GP="h3_J0_even_gp h3_J0_odd_gp h3_J1_ppos_gp h3_J1_pneg_gp h3_J2_ppos_gp h3_J2_pneg_gp dhd_gp h3_J0_even_gp_muk h3_J0_odd_gp_muk dhd_gp_muk"

nnrg_of(){ sed -n 's/^[[:space:]]*nnrg[[:space:]]*=[[:space:]]*\([0-9]*\).*/\1/p' "$SCRATCH/$1/input.d" 2>/dev/null | head -1; }
rows_done(){ local f; for f in "$SCRATCH/$1/output/ReactionProbabilities.dat" "$SCRATCH/$1/output/ReactionProbabilities-complex.dat"; do
    [ -f "$f" ] && wc -l < "$f" && return; done; echo 0; }

launch_one(){
  local r="$1" d n
  d="$SCRATCH/$r"; mkdir -p "$d/output"
  n=$(rows_done "$r")
  local need; need=$(nnrg_of "$r"); need="${need:-300}"
  if [ "$n" -ge "$need" ]; then echo "  $r complete ($n/$need rows) - skip"; return 0; fi
  cp "$PKG/inputs/$r.d" "$d/input.d"
  (cd "$d" && OPENBLAS_NUM_THREADS=$TH OMP_NUM_THREADS=$TH MKL_NUM_THREADS=$TH \
      MKL_DYNAMIC=FALSE OMP_STACKSIZE=1000M \
      setsid nohup "$BIN" < input.d > output.log 2>&1 &)
  echo "  launched $r ($n/$need rows done, threads=$TH)"
  sleep 3
}

wait_wave(){
  local wave="$1"; shift
  local need r n allok
  while :; do
    sleep 300
    allok=1
    for r in "$@"; do
      need=$(nnrg_of "$r"); need="${need:-300}"
      n=$(rows_done "$r")
      if [ "$n" -lt "$need" ]; then allok=0; fi
    done
    date; for r in "$@"; do echo "   $r: $(rows_done $r)/$(nnrg_of $r) rows"; done
    if [ $allok -eq 1 ]; then echo "  WAVE $wave DONE"; return; fi
  done
}

analyze_and_plot(){
  echo "[$(date)] analysis -> rates.npz + summary"
  python3 "$PKG/scripts/analyze_rates.py" "$SCRATCH" || echo "WARN: analyze_rates failed"
  echo "[$(date)] comparison figures"
  python3 "$PKG/scripts/plot_compare.py" || echo "WARN: plot_compare failed"
  echo "[$(date)] done. Figures: $PKG/figures/compare_*.png"
}

if [ "$MODE" = "analyze" ]; then
  analyze_and_plot
  exit 0
fi

if [ "$MODE" = "all" ] || [ "$MODE" = "ngp" ]; then
  echo "[$(date)] wave 1: NGP (fast)"
  for r in $WAVE_NGP; do launch_one "$r"; done
  wait_wave NGP $WAVE_NGP
fi
if [ "$MODE" = "all" ]; then
  echo "[$(date)] wave 2: GP (slow, days per run)"
  for r in $WAVE_GP; do launch_one "$r"; done
  wait_wave GP $WAVE_GP
fi
analyze_and_plot
