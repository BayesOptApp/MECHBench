"""Generate BBOB samples with SciPy, evaluate with IOHexperimenter, export ELA.

Run:     python ela_computation_bbob.py
Small:   python ela_computation_bbob.py --fids 1 8 --iids 1 --dimensions 5

One CSV row represents one (fid, iid, dim, seed_number, normalization).
Both raw and min-max-normalized objective values are used by default, as in
ela_computation_raw.py. Input coordinates are never normalized. All seven
feature groups from that script (including PCA) are retained. Sampling uses
SciPy rather than pflacco; feature extraction still requires pflacco.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import ioh
import numpy as np
import pandas as pd
from scipy.stats import qmc
from sklearn.preprocessing import MinMaxScaler
from pflacco.classical_ela_features import (
    calculate_ela_level,
    calculate_ela_meta,
    calculate_ela_distribution,
    calculate_information_content,
    calculate_pca,
    calculate_dispersion,
    calculate_nbc,
)

# Edit these defaults here or override them on the command line.
FUNCTION_IDS = list(range(1, 25))
#FUNCTION_IDS = list(range(1, 2))
INSTANCE_IDS = list(range(1, 6))
#INSTANCE_IDS = list(range(1, 2))
DIMENSIONS = [5, 10, 20]
#DIMENSIONS = [5]
SEEDS = list(range(1, 11))  # Random seeds for reproducibility.
#SEEDS = list(range(1, 3))  # Random seeds for reproducibility.
SAMPLE_MULTIPLIER = 250  # Number of sample points = multiplier * dimension.
OUTPUT_CSV = Path(__file__).resolve().parent / 'output_ela' / 'ela_features_bbob.csv'


def generate_sample(dimension: int, n_samples: int, seed: int,
                    method: str = 'lhs') -> np.ndarray:
    """Generate a reproducible design on the BBOB study domain [-5, 5]^d.

    The same design is reused across functions and instances for a fixed
    (dimension, sample size, seed). A seed identifies a sampling repetition;
    it is independent of the BBOB instance ID.
    """
    if dimension < 2 or n_samples < 2 or seed < 0:
        raise ValueError('Require dimension >= 2, n_samples >= 2, seed >= 0.')
    if method == 'lhs':
        unit_sample = qmc.LatinHypercube(d=dimension, seed=seed).random(n_samples)
    elif method == 'random':
        unit_sample = np.random.default_rng(seed).random((n_samples, dimension))
    else:
        raise ValueError(f'Unknown sampling method: {method}')
    return qmc.scale(unit_sample, [-5.0] * dimension, [5.0] * dimension)


def evaluate_bbob(X: np.ndarray, fid: int, iid: int) -> np.ndarray:
    """Evaluate every sample point on the requested noiseless BBOB instance."""
    problem = ioh.get_problem(
        fid=fid, instance=iid, dimension=X.shape[1],
        problem_class=ioh.ProblemClass.BBOB,
    )
    # Explicit pointwise evaluation also works on IOH versions without batching.
    y = np.asarray([problem(x.tolist()) for x in X], dtype=float)
    if not np.isfinite(y).all():
        raise ValueError(f'Non-finite objective values for fid={fid}, iid={iid}.')
    return y


def compute_ela_features(X: np.ndarray, y: np.ndarray, seed: int = 43,
                         normalize: bool = False) -> dict[str, Any]:
    """Compute every feature returned by the original seven feature groups.

    Fail explicitly if a group cannot be computed; incomplete feature rows
    are not silently exported. Constant objectives and some small designs
    may produce undefined feature values (NaN), represented as blank CSV cells.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    if X.ndim != 2 or len(X) != len(y):
        raise ValueError('X must be an (n, d) array and y must have length n.')
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError('X and y must contain only finite values.')
    
    # This is for normalization of the objective values, not the input coordinates.
    # This follows the recommendations from previous analysis on ELA features
    # Thus this allows the feature values to be comparable across different functions and instances.
    if normalize:
        y = MinMaxScaler().fit_transform(y.reshape(-1, 1)).ravel()
    groups = (
        ('ela_level', lambda: calculate_ela_level(
            X, y, ela_level_quantiles=[0.1, 0.25, 0.5])),
        ('ela_meta', lambda: calculate_ela_meta(X, y)),
        ('ela_distribution', lambda: calculate_ela_distribution(X, y)),
        ('information_content', lambda: calculate_information_content(X, y, seed=seed)),
        ('pca', lambda: calculate_pca(X, y)),
        ('dispersion', lambda: calculate_dispersion(X, y)),
        ('nbc', lambda: calculate_nbc(X, y)),
    )
    features = {}
    for name, calculate in groups:
        try:
            features.update(calculate())
        except Exception as exc:
            raise RuntimeError(f'Failed to compute {name}: {exc}') from exc
    return features


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fids', nargs='+', type=int, default=FUNCTION_IDS)
    parser.add_argument('--iids', nargs='+', type=int, default=INSTANCE_IDS)
    parser.add_argument('--dimensions', nargs='+', type=int, default=DIMENSIONS)
    parser.add_argument('--seeds', nargs='+', type=int, default=SEEDS)
    parser.add_argument('--multiplier', type=int, default=SAMPLE_MULTIPLIER)
    parser.add_argument('--n-samples', type=int, default=None,
                        help='Fixed sample count, overriding multiplier * dimension.')
    parser.add_argument('--sampling', choices=['lhs', 'random'], default='lhs')
    parser.add_argument('--normalization', choices=['both', 'raw', 'minmax'], default='both')
    parser.add_argument('--output', type=Path, default=OUTPUT_CSV)
    args = parser.parse_args()
    if any(fid < 1 or fid > 24 for fid in args.fids):
        parser.error('BBOB function IDs must be between 1 and 24.')
    if any(iid < 1 for iid in args.iids):
        parser.error('Instance IDs must be positive.')
    if any(dim < 2 for dim in args.dimensions):
        parser.error('Dimensions must be >= 2.')
    if any(seed < 0 for seed in args.seeds):
        parser.error('Seeds must be nonnegative.')
    if args.multiplier < 1 or (args.n_samples is not None and args.n_samples < 2):
        parser.error('Use multiplier >= 1 and n-samples >= 2.')
    for field in ('fids', 'iids', 'dimensions', 'seeds'):
        setattr(args, field, list(dict.fromkeys(getattr(args, field))))
    return args


def main() -> None:
    args = parse_args()
    modes = ['raw', 'minmax'] if args.normalization == 'both' else [args.normalization]
    rows = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Checkpoint after each (dimension, seed) to retain completed batches.
    for dim in args.dimensions:
        n_samples = args.n_samples if args.n_samples is not None else args.multiplier * dim
        for repetition, seed in enumerate(args.seeds):
            X = generate_sample(dim, n_samples, seed, args.sampling)
            for fid in args.fids:
                for iid in args.iids:
                    y = evaluate_bbob(X, fid, iid)
                    for mode in modes:
                        metadata = dict(fid=fid, iid=iid, dim=dim, rep=repetition + 1,
                                        seed_number=seed, n_samples=n_samples,
                                        sampling=args.sampling, normalization=mode,
                                        )
                        try:
                            features = compute_ela_features(X, y, seed, mode == 'minmax')
                        except Exception as exc:
                            raise RuntimeError(f'ELA computation failed for {metadata}') from exc
                        rows.append({**features, **metadata})
                    print(f'Computed fid={fid}, iid={iid}, dim={dim}, seed={seed}, n={n_samples}',
                          flush=True)
            pd.DataFrame(rows).to_csv(args.output, index=False)
            print(f'Saved {len(rows)} rows to {args.output.resolve()}', flush=True)


if __name__ == '__main__':
    main()
