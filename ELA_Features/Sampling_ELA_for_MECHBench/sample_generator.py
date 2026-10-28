"""Generate dim * multiplier samples using SciPy, without pflacco.

Dependencies: numpy, pandas, scipy.
Example:
    python sample_generator_scipy.py --dim 20 --multiplier 250 --seed 53
"""

import argparse
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy.stats.qmc import LatinHypercube, Sobol, scale


Bound = Union[float, Sequence[float]]


def create_sample_database(
    dim: int,
    multiplier: int,
    lower_bound: Bound = -5.0,
    upper_bound: Bound = 5.0,
    seed: Optional[int] = 43,
    sample_type: str = "lhs",
    optimization: Optional[str] = "random-cd",
) -> pd.DataFrame:
    """Return dim * multiplier samples, with columns x0, ..., x(dim-1).

    Bounds may be scalars or sequences of length dim. The index starts at 1.
    LHS uses scrambled SciPy LatinHypercube with random-cd optimization by
    default. Set optimization=None to disable optimization.
    Random sampling uses NumPy's default_rng. Sobol uses a scrambled SciPy
    sequence; its balance properties require a power-of-two sample count.
    The optimization argument applies only to LHS.
    """
    for name, value in (("dim", dim), ("multiplier", multiplier)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value <= 0:
            raise ValueError(f"{name} must be a positive integer.")

    if sample_type not in {"lhs", "random", "sobol"}:
        raise ValueError("sample_type must be 'lhs', 'random', or 'sobol'.")
    if optimization not in {None, "random-cd", "lloyd"}:
        raise ValueError("optimization must be None, 'random-cd', or 'lloyd'.")

    lower = np.broadcast_to(np.asarray(lower_bound, dtype=float), (dim,))
    upper = np.broadcast_to(np.asarray(upper_bound, dtype=float), (dim,))
    if not (np.isfinite(lower).all() and np.isfinite(upper).all()):
        raise ValueError("Bounds must be finite.")
    if np.any(lower >= upper):
        raise ValueError("Each lower bound must be strictly below its upper bound.")

    n_samples = int(dim * multiplier)
    if sample_type == "lhs":
        sampler = LatinHypercube(
            d=dim, seed=seed, optimization=optimization, scramble=True
        )
        samples = sampler.random(n=n_samples)
    elif sample_type == "sobol":
        sampler = Sobol(d=dim, seed=seed, scramble=True)
        # Keep the exact requested count. SciPy warns for non-power-of-two n.
        samples = sampler.random(n=n_samples)
    else:
        samples = np.random.default_rng(seed).random((n_samples, dim))

    samples = scale(samples, l_bounds=lower, u_bounds=upper)
    return pd.DataFrame(
        samples,
        columns=[f"x{i}" for i in range(dim)],
        index=pd.RangeIndex(1, n_samples + 1, name="id"),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dim", type=int, default=2, help="Number of dimensions")
    parser.add_argument("--multiplier", type=int, default=50,
                        help="Samples per dimension; total = dim * multiplier")
    parser.add_argument("--seed", type=int, default=43, help="Random seed for reproducibility")
    parser.add_argument("--lower_bound", type=float, default=-5.0, help="Lower bound for each dimension")
    parser.add_argument("--upper_bound", type=float, default=5.0, help="Upper bound for each dimension")
    parser.add_argument("--sample_type", choices=["lhs", "random", "sobol"], default="lhs", help="Type of sampling")
    parser.add_argument("--optimization", choices=["random-cd", "lloyd", "none"],
                        default="random-cd", help="LHS optimization method")
    parser.add_argument("--output", default=None, help="Optional output CSV path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = create_sample_database(
        dim=args.dim,
        multiplier=args.multiplier,
        lower_bound=args.lower_bound,
        upper_bound=args.upper_bound,
        seed=args.seed,
        sample_type=args.sample_type,
        optimization=None if args.optimization == "none" else args.optimization,
    )
    output = args.output or (
        f"sample_database_d_{args.dim}_mult_{args.multiplier}_seed_{args.seed}.csv"
    )
    df.to_csv(output, index=True, index_label="id", sep=" ", float_format="%.8f")
    print(f"Generated {len(df)} {args.sample_type} samples in {args.dim} dimensions: {output}")


if __name__ == "__main__":
    main()
