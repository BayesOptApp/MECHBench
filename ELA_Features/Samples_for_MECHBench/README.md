# Sample databases for MECHBench ELA

This folder contains [sample_generator.py](sample_generator.py) and pre-generated
sample databases for evaluating MECHBench at a set of design points. Each CSV
contains `dimension * multiplier` points in the MECHBench search range
`[-5, 5]`, along with an `id` column. The checked-in databases use multiplier
250 and seeds 43 through 53.

## Generate a sample database

The generator uses NumPy, pandas, and SciPy (specifically SciPy's QMC samplers).
Install the project's core dependencies, then run from the repository root:

```bash
python ELA_Features/Samples_for_MECHBench/sample_generator.py \
  --dim 20 \
  --multiplier 250 \
  --seed 43 \
  --sample_type lhs \
  --optimization random-cd \
  --output ELA_Features/Samples_for_MECHBench/250D/sample_database_d_20_mult_250_seed_43.csv
```

The generator supports `lhs`, `random`, and `sobol` sample types, configurable
bounds and seeds, and an output path. By default it creates a Latin
hypercube (`lhs`) with SciPy's `random-cd` optimization. Its total sample count
is `dim * multiplier`. For example, `--dim 20 --multiplier 250` creates 5,000
design points. The resulting CSV is whitespace-separated, has a header, and
uses 1-based row IDs.

## Run evaluations in parallel

These databases can be evaluated as independent simulations with SLURM job
arrays. Each array task selects one sample row and calls
[`main_sequential.py`](../../main_sequential.py) to evaluate that design vector.
This is a good fit for parallelization because each point's finite-element
simulation is independent: tasks can run concurrently without waiting for
previous objective values. Pre-generating the points also makes the design
fixed and repeatable across runs; choose a different seed to generate another
design.

Use the [sequential SLURM template](../../SLURM_scripts/sequential_launcher.slurm)
and configure these fields for the chosen CSV:

```bash
DIMENSION=20
MULTIPLIER=250
RANDOM_SEED=43
HEADER=1
BIAS_NUM=0
PICK_LAST=1
```

Set `DB` to the corresponding file in this directory. For example:

```bash
DB="${PROJECT_DIR}/ELA_Features/Samples_for_MECHBench/250D/sample_database_d_${DIMENSION}_mult_${MULTIPLIER}_seed_${RANDOM_SEED}.csv"
```

The template's default `DB` points to `${PROJECT_DIR}/samples/`, not this
directory, so override that line before submitting. `HEADER=1` skips the CSV
header; `PICK_LAST=1` selects the last `DIMENSION` numeric values from each
row, which excludes the numeric `id` column. The task ID plus `BIAS_NUM` gives
the 1-based sample ID and the `sim_id` passed to `main_sequential.py`.

For the full database, set the array range to the number of rows, e.g.
`#SBATCH --array=1-5000` for dimension 20 and multiplier 250. Smaller array
ranges can be used to test a subset. Keep `NUM_CORES` within the allocated
SLURM CPU resources. The launcher writes individual results by sample ID and
uses `flock` when appending summary files. Use distinct sample IDs when
submitting concurrent jobs to avoid overwriting individual result files.

## Available sample files

All listed databases are under [`250D/`](250D/); the `250D` directory denotes
the multiplier, while the dimension and seed are in each filename.

### Dimension 5

- [Seed 43](250D/sample_database_d_5_mult_250_seed_43.csv)
- [Seed 44](250D/sample_database_d_5_mult_250_seed_44.csv)
- [Seed 45](250D/sample_database_d_5_mult_250_seed_45.csv)
- [Seed 46](250D/sample_database_d_5_mult_250_seed_46.csv)
- [Seed 47](250D/sample_database_d_5_mult_250_seed_47.csv)
- [Seed 48](250D/sample_database_d_5_mult_250_seed_48.csv)
- [Seed 49](250D/sample_database_d_5_mult_250_seed_49.csv)
- [Seed 50](250D/sample_database_d_5_mult_250_seed_50.csv)
- [Seed 51](250D/sample_database_d_5_mult_250_seed_51.csv)
- [Seed 52](250D/sample_database_d_5_mult_250_seed_52.csv)
- [Seed 53](250D/sample_database_d_5_mult_250_seed_53.csv)

### Dimension 10

- [Seed 43](250D/sample_database_d_10_mult_250_seed_43.csv)
- [Seed 44](250D/sample_database_d_10_mult_250_seed_44.csv)
- [Seed 45](250D/sample_database_d_10_mult_250_seed_45.csv)
- [Seed 46](250D/sample_database_d_10_mult_250_seed_46.csv)
- [Seed 47](250D/sample_database_d_10_mult_250_seed_47.csv)
- [Seed 48](250D/sample_database_d_10_mult_250_seed_48.csv)
- [Seed 49](250D/sample_database_d_10_mult_250_seed_49.csv)
- [Seed 50](250D/sample_database_d_10_mult_250_seed_50.csv)
- [Seed 51](250D/sample_database_d_10_mult_250_seed_51.csv)
- [Seed 52](250D/sample_database_d_10_mult_250_seed_52.csv)
- [Seed 53](250D/sample_database_d_10_mult_250_seed_53.csv)

### Dimension 20

- [Seed 43](250D/sample_database_d_20_mult_250_seed_43.csv)
- [Seed 44](250D/sample_database_d_20_mult_250_seed_44.csv)
- [Seed 45](250D/sample_database_d_20_mult_250_seed_45.csv)
- [Seed 46](250D/sample_database_d_20_mult_250_seed_46.csv)
- [Seed 47](250D/sample_database_d_20_mult_250_seed_47.csv)
- [Seed 48](250D/sample_database_d_20_mult_250_seed_48.csv)
- [Seed 49](250D/sample_database_d_20_mult_250_seed_49.csv)
- [Seed 50](250D/sample_database_d_20_mult_250_seed_50.csv)
- [Seed 51](250D/sample_database_d_20_mult_250_seed_51.csv)
- [Seed 52](250D/sample_database_d_20_mult_250_seed_52.csv)
- [Seed 53](250D/sample_database_d_20_mult_250_seed_53.csv)
