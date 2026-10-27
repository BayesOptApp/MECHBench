#!/bin/bash
#SBATCH --job-name=mechbench-sequential
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-01:00:00
# Add or edit #SBATCH --partition=... and --mail-user=... for your cluster.

set -euo pipefail

# ---- Edit for your cluster/account -----------------------------------------
PROJECT_DIR="${PROJECT_DIR:-$HOME/MECHBench}"
VENV_ACTIVATE="${VENV_ACTIVATE:-$HOME/venvs/mechbench/bin/activate}"
OPENRADIOSS_MAIN_PATH="${OPENRADIOSS_MAIN_PATH:-$HOME/OpenRadioss/OpenRadioss}"

# Purge modules
module purge
# These module names are examples from the original ALICE setup. Replace them
# with modules available on your cluster, or remove the module commands if not
# needed. The virtual environment must contain this project's dependencies.
module load ALICE/default
module load OpenMPI/4.1.5-GCC-12.3.0
module load libGLU/9.0.3-GCCcore-12.3.0
module load Python/3.11.3-GCCcore-12.3.0
module load PyTorch/2.1.2-foss-2023a-CUDA-12.1.1

source "${VENV_ACTIVATE}"
cd "${PROJECT_DIR}"

# ---- Edit for this experiment ----------------------------------------------
NUM_CORES=8                     # Keep within the allocated CPU count above
PROBLEM_TYPE=3
OBJECTIVE_FUNCTION="load_uniformity"
RANDOM_SEED=50
BIAS_NUM=0                      # Add this offset to each array task ID
DIMENSION=20                    # Must match the number of design variables
HEADER=1                        # 1 if the CSV has a header row, otherwise 0
MULTIPLIER=250
PICK_LAST=1                     # 1: take last D numeric columns; 0: first D

# --- config: point to your CSV ----------------------------------------------
DB="${PROJECT_DIR}/samples/${MULTIPLIER}D/sample_database_d_${DIMENSION}_mult_${MULTIPLIER}_seed_${RANDOM_SEED}.csv"
IDX=$((${SLURM_ARRAY_TASK_ID:?Submit this script as a SLURM array job} + BIAS_NUM))
# --- derive vector for this element from CSV ---------------------------------
export DB
export IDX
export DIMENSION
export HEADER
export PICK_LAST


# Extract exactly D columns starting at COL_OFFSET from row IDX (1-based)
vector_line=$(python3 - <<'PY'
import csv, os, sys, io

DB = os.environ["DB"]
D = int(os.environ["DIMENSION"])
HEADER = int(os.environ["HEADER"])
IDX = int(os.environ["IDX"])
PICK_LAST = int(os.environ["PICK_LAST"])

def sniff_delimiter(sample: str):
    # Try csv.Sniffer; fall back to common delimiters.
    try:
        sniffer = csv.Sniffer()
        dialect = sniffer.sniff(sample, delimiters=[",",";","\t","|"])
        return dialect.delimiter
    except Exception:
        # Heuristic: prefer semicolon if many ; present, else comma, else tab.
        counts = {d: sample.count(d) for d in [",",";","\t","|"]}
        return max(counts, key=counts.get) if any(counts.values()) else ","

def is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False

# Read a small sample to sniff delimiter
with open(DB, "rb") as fb:
    raw = fb.read(32768)
text_sample = raw.decode("utf-8-sig", errors="replace")  # strip BOM if present
delim = sniff_delimiter(text_sample)

# Now read full file with the chosen delimiter
with open(DB, newline="", encoding="utf-8-sig") as f:
    rdr = csv.reader(f, delimiter=delim)
    rows = [r for r in rdr if any(c.strip() for c in r)]  # drop blank lines

if not rows:
    sys.exit("CSV seems empty after removing blank lines.")

row_index = (IDX - 1) + HEADER
if row_index < 0 or row_index >= len(rows):
    sys.exit(f"Requested row {IDX} (data index {row_index}) is out of range; file has {len(rows)-HEADER} data rows.")

row = [c.strip() for c in rows[row_index]]

# If the row collapsed into one giant field (wrong delimiter), try whitespace split as a last resort
if len(row) == 1 and " " in row[0]:
    row = [tok for tok in row[0].split() if tok]

# Collect numeric fields only (ignores IDs / strings)
nums = [c for c in row if is_float(c)]

if len(nums) < D:
    # Helpful diagnostics
    msg = (
        f"Row has only {len(nums)} numeric fields but D={D}.\n"
        f"Detected delimiter='{delim}', row length={len(row)}.\n"
        f"Row preview: {row[:min(12,len(row))]}"
    )
    sys.exit(msg)

vals = nums[-D:] if PICK_LAST else nums[:D]
print(" ".join(vals), end="")
PY
)

if [[ -z "${vector_line}" ]]; then
  echo "Failed to extract vector for IDX=${IDX}. See log above." >&2
  exit 1
fi

# --- run your simulation ------------------------------------------------------
sim_id="${IDX}" # To be adjusted

echo "Running simulation with sim_id=${sim_id}, vector=${vector_line}"

mkdir -p results

# Start timing
job_start=$(date +%s)

# Run and keep only the final printed line (your float)
last_line=$(python3 -u "${PROJECT_DIR}/main_sequential.py" \
  --main_path "${OPENRADIOSS_MAIN_PATH}" \
  --sim_id "${sim_id}" \
  --vector ${vector_line} \
  --problem_type "${PROBLEM_TYPE}" \
  --objective_function "${OBJECTIVE_FUNCTION}" \
  --np "${NUM_CORES}" --nt 1 --h_level 1 --gmsh_verbosity 0 \
  | tail -n 1)

# End timing
job_end=$(date +%s)
job_runtime=$((job_end - job_start))

echo "Job finished at: $(date)"
echo "Total SLURM runtime: ${job_runtime} seconds"

# Per-element result file (just the float)
RESULTS_DIR="${PROJECT_DIR}/results/individual_points/problem_${PROBLEM_TYPE}"
mkdir -p "${RESULTS_DIR}"
echo "${last_line}" > "${RESULTS_DIR}/result_dim_${DIMENSION}_seed_${RANDOM_SEED}_${sim_id}.txt"
echo "${job_runtime}" > "${RESULTS_DIR}/result_dim_${DIMENSION}_seed_${RANDOM_SEED}_${sim_id}_runtime.txt"

# Append to global summary safely (requires 'flock' from util-linux)
mkdir -p results
{
  exec 200>"${RESULTS_DIR}/summary_dim_${DIMENSION}_seed_${RANDOM_SEED}.lock"
  flock 200
  printf "%s\t%s\n" "${sim_id}" "${last_line}" >> "${RESULTS_DIR}/summary_dim_${DIMENSION}_seed_${RANDOM_SEED}.txt"
  printf "%s\t%s\n" "${sim_id}" "${job_runtime}" >> "${RESULTS_DIR}/summary_dim_${DIMENSION}_seed_${RANDOM_SEED}_runtime.txt"
}
