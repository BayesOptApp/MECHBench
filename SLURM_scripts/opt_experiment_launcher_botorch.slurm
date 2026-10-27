#!/bin/bash
#SBATCH --job-name=mechbench-opt-gpu
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-10
#SBATCH --gres=gpu:1
# Add or edit #SBATCH --partition=... and --mail-user=... for your cluster.

set -euo pipefail

# ---- Edit for your cluster/account -----------------------------------------
PROJECT_DIR="${PROJECT_DIR:-$HOME/MECHBench}"
VENV_ACTIVATE="${VENV_ACTIVATE:-$HOME/venvs/mechbench/bin/activate}"
OPENRADIOSS_MAIN_PATH="${OPENRADIOSS_MAIN_PATH:-$HOME/OpenRadioss/OpenRadioss}"

# These module names are examples from the original ALICE setup. Replace them
# with modules available on your cluster, or remove the module commands if not
# needed. The virtual environment must contain this project's dependencies.
module purge
module load ALICE/default
module load OpenMPI/4.1.5-GCC-12.3.0
module load libGLU/9.0.3-GCCcore-12.3.0
module load Python/3.11.3-GCCcore-12.3.0
module load PyTorch/2.1.2-foss-2023a-CUDA-12.1.1

source "${VENV_ACTIVATE}"
cd "${PROJECT_DIR}"

# ---- Edit for this experiment ----------------------------------------------
ALGORITHM="BAxUS_botorch"       # pyCMA, BO_botorch, BAxUS_botorch, HEBO, turbo1
DIMENSION=20                    # Supported dimensions: 2, 5, 10, 20
PROBLEM_TYPE=2                  # 1, 2, or 3
OBJECTIVE_FUNCTION="penalized_mass"
MAX_ITER=600
DOE_MULT=3
NP=8                            # Keep aligned with #SBATCH --ntasks above
NT=1
H_LEVEL=1

RUN_ID="${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM array job}"

echo "MECHBench GPU: optimizer=${ALGORITHM}, dimension=${DIMENSION}, problem=${PROBLEM_TYPE}, run=${RUN_ID}, max_iter=${MAX_ITER}, np=${NP}, objective=${OBJECTIVE_FUNCTION}"
python3 "${PROJECT_DIR}/main_experiments.py" \
  --main_path "${OPENRADIOSS_MAIN_PATH}" \
  --problem_type "${PROBLEM_TYPE}" \
  --dimension "${DIMENSION}" \
  --objective_function "${OBJECTIVE_FUNCTION}" \
  --optimizer "${ALGORITHM}" \
  --run_id "${RUN_ID}" \
  --max_iter "${MAX_ITER}" \
  --doe_mult "${DOE_MULT}" \
  --np "${NP}" \
  --nt "${NT}" \
  --h_level "${H_LEVEL}"
