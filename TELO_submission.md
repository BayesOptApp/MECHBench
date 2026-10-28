# TELO submission environment

This guide explains how to install MECHBench and the optional dependency groups
declared in `pyproject.toml`.

## Install MECHBench

From the repository root, create and activate a Python environment, then install
the project and its core dependencies:

```bash
python -m venv .venv
```

Activate it before continuing:

```bash
# Linux/macOS
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install MECHBench in editable mode:

```bash
python -m pip install -e .
```

## Install optional dependencies

Optional groups are declared under `[project.optional-dependencies]` in
`pyproject.toml`. Install a group by appending its name in square brackets to
the project install command:

```bash
# HDBO and related optimizer dependencies
python -m pip install -e ".[hdbo]"

# ELA (Exploratory Landscape Analysis) dependencies
python -m pip install -e ".[ela]"

# Both optional groups
python -m pip install -e ".[hdbo,ela]"
```

The `hdbo` group includes dependencies for the HDBO and related optimizers,
including BoTorch, GPyTorch, PyTorch, HEBO, and CMA-ES. The `ela` group installs
`pflacco`.

If MECHBench is already installed in the active environment, rerunning the
appropriate command installs the selected optional dependencies as well.

## Reproduce the BBOB experiments

The repository includes two BBOB workflows that can be run from the repository
root. They use IOHexperimenter test functions directly and do not require
OpenRadioss:

- [`main_experiments_bbob.py`](main_experiments_bbob.py) runs one black-box
  optimization experiment and writes IOH logs.
- [`ELA_Features/Sampling_ELA_with_BBOB/ela_computation_bbob.py`](ELA_Features/Sampling_ELA_with_BBOB/ela_computation_bbob.py)
  evaluates reproducible designs and computes ELA feature rows.

Install the corresponding optional dependencies before running them:

```bash
# BBOB optimizer experiments
python -m pip install -e ".[hdbo]"

# BBOB ELA computation
python -m pip install -e ".[ela]"

# Both workflows in the same environment
python -m pip install -e ".[hdbo,ela]"
```

### Compute BBOB ELA features

From the repository root, run the default feature experiment:

```bash
python ELA_Features/Sampling_ELA_with_BBOB/ela_computation_bbob.py
```

The default run covers function IDs 1–24, instances 1–5, dimensions 5, 10,
and 20, and sampling seeds 1–10. It uses Latin-hypercube designs with 250
points per dimension and computes both raw-objective and min-max-normalized
feature rows, for 7,200 rows total. Results are saved to
`ELA_Features/Sampling_ELA_with_BBOB/output_ela/ela_features_bbob.csv`.
Completed dimension/seed batches are checkpointed to the CSV while the run is
in progress.

For a quick smoke test before the full run:

```bash
python ELA_Features/Sampling_ELA_with_BBOB/ela_computation_bbob.py \
  --fids 1 8 \
  --iids 1 \
  --dimensions 5 \
  --seeds 43 \
  --multiplier 10 \
  --normalization both \
  --output ELA_Features/Sampling_ELA_with_BBOB/output_ela/ela_features_smoke_test.csv
```

This smaller run produces four rows (two functions × one instance × one
dimension × one seed × two normalization modes). Override `--fids`, `--iids`,
`--dimensions`, `--seeds`, `--multiplier`, `--sampling`, or
`--normalization` to select another experiment. The script accepts sampling
methods `lhs` and `random`; it saves raw and min-max features by default.
See the [BBOB ELA README](ELA_Features/Sampling_ELA_with_BBOB/README.md) for
the full parameter and output description.

### Run BBOB optimization experiments

`main_experiments_bbob.py` runs one optimizer/function/instance/dimension
combination per invocation. For example, run CMA-ES (`pyCMA`) on BBOB
function 1, instance 1, in dimension 5:

```bash
python main_experiments_bbob.py \
  --optimizer pyCMA \
  --function_id 1 \
  --instance 1 \
  --dimension 5 \
  --run_id 0 \
  --max_iter 100 \
  --doe_mult 3 \
  --seed 42 \
  --ioh_logger_folder logdata_optimization
```

The script supports optimizers `pyCMA`, `BO_botorch`, `BAxUS_botorch`, `HEBO`,
and `turbo1`; BBOB function IDs 1–24; dimensions 2, 5, 10, and 20; and
positive instance IDs. `--run_id` is required. Each run uses optimizer seed
`--seed + --run_id`; use a distinct run ID for each repetition. The initial
design size is `--doe_mult * --dimension` and must not exceed `--max_iter`,
which is the total evaluation budget. Set the parameters to the values used
in the experiment you intend to reproduce.

For example, a small functional test can use a short budget:

```bash
python main_experiments_bbob.py \
  --optimizer pyCMA \
  --function_id 1 \
  --instance 1 \
  --dimension 2 \
  --run_id 0 \
  --max_iter 10 \
  --doe_mult 2 \
  --ioh_logger_folder logdata_optimization/smoke_test
```

To run independent repetitions locally, use distinct run IDs. This example
runs ten repetitions with the same function, instance, dimension, optimizer,
and budget:

```bash
for run_id in $(seq 0 9); do
  python main_experiments_bbob.py \
    --optimizer pyCMA \
    --function_id 1 \
    --instance 1 \
    --dimension 5 \
    --run_id "$run_id" \
    --max_iter 100 \
    --doe_mult 3 \
    --seed 42 \
    --ioh_logger_folder logdata_optimization
done
```

IOH logs are written beneath the selected logger root (by default
`logdata_optimization` in the current working directory). Each run folder
records the function, instance, dimension, run ID, and optimizer; the logger
stores positions and the raw best-so-far values. Give each concurrent run a
unique combination of run ID and logger folder to avoid output collisions.
`main_experiments_bbob.py` imports `wrapopt` from the repository's
`ML_Benchmarks` package, so run it from the repository root in the environment
where the required optimizer dependencies are installed. This BBOB workflow
does not launch the MECHBench finite-element models.

## PyTorch and cluster environments

For GPU runs, install a PyTorch build compatible with the cluster's CUDA
version. On clusters that provide PyTorch as a module, follow the cluster's
environment instructions and confirm the activated Python environment can
import `torch`; a module and a virtual environment do not always share packages.

The repository's `requirements.txt` lists core dependencies. To install optional
groups, use the `pyproject.toml` commands above rather than assuming
`pip install -r requirements.txt` includes them.
