# MECHBench SLURM templates

These scripts are starting points for running MECHBench on a SLURM cluster:

- `opt_experiment_launcher.slurm`: CPU optimization runs using `main_experiments.py`.
- `opt_experiment_launcher_botorch.slurm`: GPU-requesting optimization runs using `main_experiments.py`.
- `sequential_launcher.slurm`: evaluates rows from a CSV sample database using `main_sequential.py`.

The scripts are templates, not cluster-independent job specifications. Update
the cluster settings and paths below before submitting them. Submit from the
repository root, for example:

```bash
sbatch SLURM_scripts/opt_experiment_launcher.slurm
```

Each optimization array task uses its `SLURM_ARRAY_TASK_ID` as the run ID. The
sequential launcher uses it as the sample row number, plus `BIAS_NUM`.

## Fields to edit

| Field | Where | What to set |
| --- | --- | --- |
| `#SBATCH --job-name` | All scripts | A recognizable name for the job. |
| `#SBATCH --output`, `--error` | All scripts | Log file paths/patterns. The defaults write to the submission directory; any parent directory in a custom path must already exist when submitting. |
| `#SBATCH --array` | All scripts | The task IDs to run. Use the range/count appropriate to the experiment and data. For sequential runs, ensure each resulting `IDX` is a valid CSV data row. |
| `#SBATCH --ntasks`, `--cpus-per-task`, `--mem-per-cpu`, `--time` | All scripts | Resource requests supported by the cluster. Keep `NP` (optimization) or `NUM_CORES` (sequential) within the allocated CPU capacity. |
| `#SBATCH --partition` | All scripts (optional) | Uncomment/add a partition directive using a partition name on your cluster. |
| `#SBATCH --gres=gpu:1` | GPU optimization script | Change or remove the GPU request to match your cluster's GPU syntax and the experiment's needs. |
| `#SBATCH --mail-user`, `--mail-type` | All scripts (optional) | Add your email address and desired notifications if the cluster supports job email. |
| `PROJECT_DIR` | All scripts | Absolute path to the MECHBench checkout on the cluster. |
| `VENV_ACTIVATE` | All scripts | Path to the Python environment's `bin/activate` script. |
| `OPENRADIOSS_MAIN_PATH` | All scripts | OpenRadioss installation path used by MECHBench. |
| `module load ...` | All scripts | Replace the example ALICE module names with modules available on the target cluster, or remove these lines if the environment is already configured. |

### Optimization experiment fields

Edit the experiment values in either optimization script:

- `ALGORITHM`: `pyCMA`, `BO_botorch`, `BAxUS_botorch`, `HEBO`, or `turbo1`.
- `DIMENSION`: supported values are `2`, `5`, `10`, and `20`.
- `PROBLEM_TYPE`: `1`, `2`, or `3`.
- `OBJECTIVE_FUNCTION`: one of `load_uniformity`, `intrusion`, `absorbed_energy`,
  `mass`, `specific-energy-absorbed`, `penalized_sea`, or `penalized_mass`.
- `MAX_ITER`, `DOE_MULT`, `NP`, `NT`, and `H_LEVEL`: optimizer and simulation
  settings. Match `NP` and the SLURM task allocation, and adjust `NT` and
  `--cpus-per-task` consistently if using threaded simulations.

The GPU script requests one GPU but does not force the selected optimizer or
the Python code to use it; choose an appropriate algorithm and environment.

### Sequential sample fields

Edit these values in `sequential_launcher.slurm`:

- `PROBLEM_TYPE` and `OBJECTIVE_FUNCTION`: benchmark problem and output to evaluate.
- `DIMENSION`, `RANDOM_SEED`, and `MULTIPLIER`: used to construct the default
  CSV filename under `${PROJECT_DIR}/samples/${MULTIPLIER}D/`.
- `DB`: override the CSV path if the file does not use the default naming pattern.
- `HEADER`: `1` when the CSV has a header row, otherwise `0`.
- `BIAS_NUM`: offset added to each array task ID to get the 1-based data row.
- `PICK_LAST`: `1` selects the last `DIMENSION` numeric CSV values; `0` selects
  the first `DIMENSION` values.
- `NUM_CORES`: number of simulation processes, within the SLURM allocation.

The sequential script writes per-sample results and summary files under
`${PROJECT_DIR}/results/individual_points/problem_${PROBLEM_TYPE}/`. Avoid
submitting overlapping jobs that write the same sample IDs unless that is
intended.
