#!/usr/bin/env bash
# ==============================================================================
# One-Click Deployment, Environment Self-Check & Scientific Validation
# ==============================================================================

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "================================================================================"
echo "      PINNSlove: PRODUCTION DEPLOYMENT & ENVIRONMENT INTEGRITY VERIFIER         "
echo "================================================================================"
echo "  Timestamp : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Host      : $(hostname) ($(uname -s) $(uname -m))"
echo "  CPUs      : $(nproc 2>/dev/null || echo 'Unknown') cores detected"
echo "  RAM       : $(free -h 2>/dev/null | awk '/Mem:/ {print $2}' || echo 'Unknown')"
echo "================================================================================"

# 1. Check Python environment
echo ""
echo "[1/4] Verifying Python3 and Scientific Libraries..."
python3 -c "
import numpy, scipy, matplotlib
print(f'  Python : {__import__(\"sys\").version.split()[0]}')
print(f'  NumPy  : {numpy.__version__}')
print(f'  SciPy  : {scipy.__version__}')
print(f'  Matplotlib: {matplotlib.__version__}')
" || { echo "ERROR: Missing required Python scientific libraries."; exit 1; }

# 2. Check BKMP2 Fortran potential library
echo ""
echo "[2/4] Verifying Authentic BKMP2 Potential Library (peslib/libbkmp2.so)..."
python3 -c "
import sys; sys.path.insert(0, '.')
from src import bkmp2
lib = bkmp2.get_lib(verbose=False)
v = bkmp2.V_h3_angstrom(0.74144, 30.0, 30.0, cm1=True)
print(f'  Loaded Fortran Library: {lib}')
print(f'  BKMP2 Isolated Asymptote: {v:.2f} cm^-1 (Exact: -38297.39 cm^-1)')
assert abs(v - (-38297.39)) < 0.1, 'BKMP2 potential check failed!'
" || { echo "ERROR: BKMP2 potential library validation failed."; exit 1; }

# 3. Run Scientific Audit
echo ""
echo "[3/4] Executing Scientific Audit..."
python3 scripts/validation/audit_sub_permill_all.py || { echo "ERROR: Scientific audit failed."; exit 1; }

# 4. Run Rigorous Regression Suite
echo ""
echo "[4/4] Executing Clean Regression Suite..."
python3 scripts/validation/run_full_regression.py || { echo "ERROR: Regression suite failed."; exit 1; }

echo ""
echo "================================================================================"
echo "  ENVIRONMENT & DEPLOYMENT VERIFICATION SUCCESSFUL (100% GREEN)"
echo "================================================================================"
