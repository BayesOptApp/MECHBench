#!/bin/bash
#SBATCH --job-name=problem_3_rs_50_dim_20    
#SBATCH --mail-user=olarterodriguezi@vuw.leidenuniv.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/olarterodriguezi/logs/experiments_mult_250_dim_20_problem_3_rs_50_%A_%a.out
#SBATCH --error=/home/olarterodriguezi/logs/experiments_mult_250_dim_20_problem_3_rs_50_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=9
#SBATCH --nodes=1
#SBATCH --array=1-1000
#SBATCH --mem-per-cpu=2048M
#SBATCH --partition=cpu-zen4
#SBATCH --time=0-00:20:00



# ###--array=1
#djslkadSBATCH --array=1-600
# partition=cpu-zen4
# partition=cpu-skylake

# Purge modules
module purge
# load modules (assuming you start from the default environment)
# we explicitly call the modules to improve reproducibility
# in case the default settings change

#FOSS
module load ALICE/default
#module load slurm

# ACTIVATE OPENMPI
module load OpenMPI/4.1.5-GCC-12.3.0

#Activate libGLU (This is a helper if ALICE/default does not load it)
module load libGLU/9.0.3-GCCcore-12.3.0

#PYTHON
module load Python/3.11.3-GCCcore-12.3.0

#PYTORCH
module load PyTorch/2.1.2-foss-2023a-CUDA-12.1.1

# Activate the Python virtual environment
source /home/olarterodriguezi/BO_torch/bo-torch-run/bin/activate

# Change source directory to the home directory
cd /home/olarterodriguezi/MECHBench

# Add executable permissions to the open_radioss binaries
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/starter_linux64_gf
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/engine_linux64_gf
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/engine_linux64_gf_ompi
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/th_to_csv_linux64_gf
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/anim_to_vtk_linux64_gf


#OBJECTIVE_FUNCTION="penalized_sea"
#OBJECTIVE_FUNCTION="penalized_mass"
#OBJECTIVE_FUNCTION="load_uniformity"
#OBJECTIVE_FUNCTION="mass intrusion penalized_mass"
#OBJECTIVE_FUNCTION="specific_energy_absorbed intrusion penalized_sea"

# Specify the path to the config file
NUM_CORES=8
PROBLEM_TYPE=3
OBJECTIVE_FUNCTION="load_uniformity"
#OBJECTIVE_FUNCTION="mass"
#OBJECTIVE_FUNCTION="penalized_mass"
#OBJECTIVE_FUNCTION="specific_energy_absorbed"
RANDOM_SEED=50
BIAS_NUM=3000
DIMENSION=20
HEADER=1  # Set to 1 if the CSV has a header row, 0 otherwise
COL_OFFSET=1
MULTIPLIER=250
IDX=$(($SLURM_ARRAY_TASK_ID + $BIAS_NUM))
#IDX=1

set -euo pipefail

# --- config: point to your CSV ----------------------------------------------
DB="/home/olarterodriguezi/MECHBench/samples/${MULTIPLIER}D/sample_database_d_${DIMENSION}_mult_${MULTIPLIER}_seed_${RANDOM_SEED}.csv"   # path to the  CSV
#   /home/olarterodriguezi/MECHBench/samples/30D/sample_database_d_5_mult_30_seed_43.csv
# --- derive vector for this element from CSV ---------------------------------
export DB
export IDX
export DIMENSION
export HEADER
export PICK_LAST=1  # Set to 1 to pick the last D columns, 0 to pick the first D columns


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
last_line=$(python3 -u main_sequential.py \
  --sim_id "${sim_id}" \
  --vector ${vector_line} \
  --problem_type ${PROBLEM_TYPE} \
  --objective_function ${OBJECTIVE_FUNCTION} \
  --np ${NUM_CORES} --nt 1 --h_level 1 --gmsh_verbosity 0 \
  | tail -n 1)

# End timing
job_end=$(date +%s)
job_runtime=$((job_end - job_start))

echo "Job finished at: $(date)"
echo "Total SLURM runtime: ${job_runtime} seconds"

# Per-element result file (just the float)
echo "${last_line}" > "/home/olarterodriguezi/MECHBench/results/individual_points/problem_${PROBLEM_TYPE}/result_dim_${DIMENSION}_seed_${RANDOM_SEED}_${sim_id}.txt"
echo "${job_runtime}" > "/home/olarterodriguezi/MECHBench/results/individual_points/problem_${PROBLEM_TYPE}/result_dim_${DIMENSION}_seed_${RANDOM_SEED}_${sim_id}_runtime.txt"

# Append to global summary safely (requires 'flock' from util-linux)
mkdir -p results
{
  exec 200>/home/olarterodriguezi/MECHBench/results/individual_points/problem_${PROBLEM_TYPE}/summary_dim_${DIMENSION}_seed_${RANDOM_SEED}.lock
  flock 200
  printf "%s\t%s\n" "${sim_id}" "${last_line}" >> /home/olarterodriguezi/MECHBench/results/individual_points/problem_${PROBLEM_TYPE}/summary_dim_${DIMENSION}_seed_${RANDOM_SEED}.txt
  printf "%s\t%s\n" "${sim_id}" "${job_runtime}" >> /home/olarterodriguezi/MECHBench/results/individual_points/problem_${PROBLEM_TYPE}/summary_dim_${DIMENSION}_seed_${RANDOM_SEED}_runtime.txt
}
