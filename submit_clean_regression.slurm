#!/usr/bin/env bash
#SBATCH --job-name=pinnslove_clean
#SBATCH --partition=air
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_pinnslove_%j.out
#SBATCH --error=slurm_pinnslove_%j.err

set -e

echo "================================================================================"
echo "    SLURM JOB LAUNCHED: PINNSLOVE CLEAN PIPELINE                               "
echo "================================================================================"
echo "  Job ID    : ${SLURM_JOB_ID}"
echo "  Job Name  : ${SLURM_JOB_NAME}"
echo "  Partition : ${SLURM_JOB_PARTITION}"
echo "  Node List : ${SLURM_JOB_NODELIST}"
echo "  CPUs/Task : ${SLURM_CPUS_PER_TASK}"
echo "  Start Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================================"

cd "${SLURM_SUBMIT_DIR:-/media/lh/新加卷/ProTest/PINNSlove}"

export PATH="${HOME}/.local/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/.local/lib:${LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1

# 1. Environment and BKMP2 library check
python3 scripts/check_environment.py

# 2. Strict sub-per-mille scientific audit
python3 scripts/audit_sub_permill_all.py

# 3. Full 9-level regression test suite
python3 scripts/run_full_regression.py

# 4. State-to-state multichannel scattering benchmarks
python3 scripts/run_state_to_state_benchmarks.py

# 5. Cryptographic manifest generation and verification
python3 scripts/generate_manifest.py
python3 scripts/verify_manifest.py

echo "================================================================================"
echo "  SLURM JOB ${SLURM_JOB_ID} COMPLETED SUCCESSFULLY AT: $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================================"
