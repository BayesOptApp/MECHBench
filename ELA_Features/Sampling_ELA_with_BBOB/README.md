# Generating ELA features from BBOB

[ela_computation_bbob.py](ela_computation_bbob.py) generates reproducible
sample designs in the BBOB domain, evaluates them with IOH, computes ELA
features with `pflacco`, and writes the feature rows to a CSV. Its default
configuration is intended to reproduce the feature-generation setup used for
the paper: the same seven ELA feature groups (including PCA), computed for raw
and min-max-normalized objective values. The input coordinates are not
normalized.

The script uses SciPy for sampling rather than `pflacco`; ELA feature
calculation itself requires `pflacco`.

## Install dependencies

From the repository root, install MECHBench with the ELA extra. This will install `scikit-learn` to use the `min-max` scaler from
such library.

```bash
python -m pip install -e ".[ela]"
```

The script uses IOH, NumPy, pandas, SciPy, scikit-learn, and pflacco. Use a
Python environment compatible with the installed package versions.

## Reproduce the full feature dataset

Run the script from the repository root:

```bash
python ELA_Features/Sampling_ELA_with_BBOB/ela_computation_bbob.py
```

The defaults are:

| Setting | Default |
| --- | --- |
| BBOB function IDs (`--fids`) | 1–24 |
| BBOB instance IDs (`--iids`) | 1–5 |
| Dimensions (`--dimensions`) | 5, 10, 20 |
| Sampling seeds (`--seeds`) | 43 |
| Sampling method (`--sampling`) | Latin hypercube (`lhs`) |
| Samples per dimension (`--multiplier`) | 250 |
| Objective normalization (`--normalization`) | Both `raw` and `minmax` |
| Output | `ELA_Features/Sampling_ELA_with_BBOB/output_ela/ela_features_bbob.csv` |

With these settings, the script produces 720 feature rows: 24 functions × 5
instances × 3 dimensions × 1 sampling seed × 2 normalization modes. For each
dimension and sampling seed, it generates a design of `multiplier * dimension`
points in `[-5, 5]^dimension`. That same design is evaluated on each selected
function and instance for the dimension/seed combination. The sampling seed
controls the design and is independent of the BBOB instance ID.

To specify the defaults explicitly and choose an output path:

```bash
python ELA_Features/Sampling_ELA_with_BBOB/ela_computation_bbob.py \
  --fids 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 \
  --iids 1 2 3 4 5 \
  --dimensions 5 10 20 \
  --seeds 43 \
  --multiplier 250 \
  --sampling lhs \
  --normalization both \
  --output ELA_Features/Sampling_ELA_with_BBOB/output_ela/ela_features_bbob.csv
```

The output file is created if needed and rewritten as the run proceeds. The
script saves after each dimension/seed batch, so completed batches are retained
if a later batch fails. Re-running the command starts the computation again
and overwrites that output file.

## Test with a smaller run

Before computing the full dataset, run a small subset:

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

This computes 4 rows (2 functions × 1 instance × 1 dimension × 1 seed × 2
normalization modes) using 50 sample points per function evaluation.

## Normalization formula

For the `minmax` mode, the script fits a separate scaler to the objective
values from each function/instance/design evaluation. For objective value
$y_i$, the transformation is:

$$
y_i^{\mathrm{scaled}} =
\frac{y_i - y_{\min}}{y_{\max} - y_{\min}},
\qquad
y_{\min} = \min_j y_j,\quad y_{\max} = \max_j y_j.
$$

Thus non-constant objective values map to the interval \([0, 1]\). The `raw`
mode uses the original objective values. When all values are equal,
scikit-learn's `MinMaxScaler` maps the constant values to zero. Only the
objective values are scaled; the sample coordinates remain in `[-5, 5]`.

## Feature groups and output

For each evaluation the script computes seven feature groups:

- ELA level
- ELA meta
- ELA distribution
- Information content
- PCA
- Dispersion
- Nearest-better clustering (NBC)

The CSV starts with metadata columns (`fid`, `iid`, `dim`, `seed_number`,
`n_samples`, `sampling`, and `normalization`), followed by the feature values.
Undefined feature values may appear as blank CSV cells; a feature-group
computation failure stops the run with an error rather than exporting an
incomplete feature row.

The repository includes an [example generated dataset](output_ela/ela_features_bbob.csv).
