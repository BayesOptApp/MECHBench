"""Compute exploratory landscape analysis features from MECHBench evaluations.

Sample rows and evaluation rows are joined by their 1-based simulation ID.
See the adjacent README for the expected input formats and examples.
"""

from __future__ import annotations

import argparse
import ast
import csv
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from pflacco.classical_ela_features import (
    calculate_dispersion,
    calculate_ela_distribution,
    calculate_ela_level,
    calculate_ela_meta,
    calculate_information_content,
    calculate_nbc,
    calculate_pca,
)
from sklearn.preprocessing import MinMaxScaler


DEFAULT_SAMPLE_FILE = (
    Path(__file__).resolve().parent
    / "250D"
    / "sample_database_d_20_mult_250_seed_43.csv"
)
DEFAULT_EVALUATION_DIR = Path("results") / "individual_points"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output_ela" / "ela_features_mechbench.csv"
SAMPLE_NAME_PATTERN = re.compile(
    r"sample_database_d_(?P<dimension>\d+)_mult_(?P<multiplier>\d+)_seed_(?P<seed>\d+)",
    re.IGNORECASE,
)
SUMMARY_NAME_PATTERN = re.compile(
    r"summary_dim_(?P<dimension>\d+)_seed_(?P<seed>\d+)",
    re.IGNORECASE,
)


def _read_delimited_table(path: Path) -> pd.DataFrame:
    """Read a CSV/TSV/whitespace table, detecting common delimiters."""
    try:
        return pd.read_csv(path, sep=None, engine="python", comment="#")
    except (pd.errors.ParserError, UnicodeDecodeError):
        return pd.read_csv(path, sep=r"\s+", engine="python", comment="#")


def _numeric_id(value: Any, source: Path) -> int:
    try:
        numeric_value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid simulation ID {value!r} in {source}.") from exc
    if not np.isfinite(numeric_value) or not numeric_value.is_integer():
        raise ValueError(f"Simulation IDs must be finite integers; got {value!r} in {source}.")
    return int(numeric_value)


def load_samples(path: Path) -> Tuple[List[int], np.ndarray]:
    """Load design points from a delimited file with an ID column or index."""
    table = _read_delimited_table(path)
    if table.empty:
        raise ValueError(f"No sample rows found in {path}.")

    id_column = next(
        (column for column in table.columns if str(column).strip().lower() in {"id", "sim_id"}),
        None,
    )
    if id_column is None and str(table.columns[0]).strip().lower().startswith("unnamed"):
        id_column = table.columns[0]

    if id_column is None:
        ids = list(range(1, len(table) + 1))
    else:
        ids = [_numeric_id(value, path) for value in table.pop(id_column)]

    try:
        points = table.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Sample file must contain only numeric design columns: {path}") from exc

    if points.ndim != 2 or points.shape[1] == 0:
        raise ValueError(f"No design-variable columns found in {path}.")
    if not np.isfinite(points).all():
        raise ValueError(f"Sample points in {path} must be finite.")
    if len(set(ids)) != len(ids):
        raise ValueError(f"Sample IDs must be unique in {path}.")
    return ids, points


def _evaluation_value(value: str, source: Path, line_number: int) -> float:
    """Read a scalar objective or the final value in a legacy list value."""
    text = value.strip()
    if not text:
        raise ValueError(f"Missing objective value at {source}:{line_number}.")
    if text.startswith("["):
        try:
            values = ast.literal_eval(text)
        except (SyntaxError, ValueError) as exc:
            raise ValueError(f"Invalid objective list at {source}:{line_number}.") from exc
        if not isinstance(values, (list, tuple)) or not values:
            raise ValueError(f"Expected a non-empty objective list at {source}:{line_number}.")
        text = str(values[-1])
    try:
        result = float(text)
    except ValueError as exc:
        raise ValueError(f"Invalid objective value {value!r} at {source}:{line_number}.") from exc
    if not np.isfinite(result):
        raise ValueError(f"Objective values must be finite at {source}:{line_number}.")
    return result


def load_evaluations(path: Path) -> Dict[int, float]:
    """Load MECHBench summary rows (simulation ID, objective value).

    Also accepts two-column CSV/TSV files with ``sim_id``/``id`` and
    ``objective``/``objective_value`` columns, or positional columns.
    """
    evaluations: Dict[int, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter="\t"), start=1):
            if not row or not any(field.strip() for field in row):
                continue
            if len(row) < 2:
                # Permit comma-separated summaries as well.
                row = next(csv.reader([row[0]], delimiter=","))
            if len(row) < 2:
                raise ValueError(f"Expected an ID and objective value at {path}:{line_number}.")

            id_text = row[0].strip()
            value_text = row[1].strip()
            if line_number == 1 and id_text.lower() in {"id", "sim_id", "evaluation_id"}:
                continue
            simulation_id = _numeric_id(id_text, path)
            if simulation_id in evaluations:
                raise ValueError(f"Duplicate simulation ID {simulation_id} in {path}.")
            evaluations[simulation_id] = _evaluation_value(value_text, path, line_number)

    if not evaluations:
        raise ValueError(f"No evaluation rows found in {path}.")
    return evaluations


def discover_summaries(path: Path) -> List[Path]:
    """Return one summary file or recursively discover MECHBench summaries."""
    if path.is_file():
        return [path]
    if not path.is_dir():
        raise FileNotFoundError(f"Evaluation path does not exist: {path}")
    summaries = sorted(path.rglob("summary_dim_*_seed_*.txt"))
    if not summaries:
        raise FileNotFoundError(
            f"No MECHBench summary files found under {path}. "
            "Expected summary_dim_<dimension>_seed_<seed>.txt."
        )
    return summaries


def compute_ela_features(
    X: np.ndarray,
    y: np.ndarray,
    seed: int = 43,
    normalize: bool = False,
) -> Dict[str, Any]:
    """Compute seven pflacco feature groups for one design/evaluation pair."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    if X.ndim != 2 or X.shape[0] != y.size:
        raise ValueError("X must have shape (n_samples, dimension) and y length n_samples.")
    if X.shape[0] < 2:
        raise ValueError("At least two point/evaluation pairs are required.")
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("X and y must contain only finite values.")
    if normalize:
        y = MinMaxScaler().fit_transform(y.reshape(-1, 1)).ravel()

    groups = (
        ("ela_level", lambda: calculate_ela_level(X, y, ela_level_quantiles=[0.1, 0.25, 0.5])),
        ("ela_meta", lambda: calculate_ela_meta(X, y)),
        ("ela_distribution", lambda: calculate_ela_distribution(X, y)),
        ("information_content", lambda: calculate_information_content(X, y, seed=seed)),
        ("pca", lambda: calculate_pca(X, y)),
        ("dispersion", lambda: calculate_dispersion(X, y)),
        ("nbc", lambda: calculate_nbc(X, y)),
    )
    features: Dict[str, Any] = {}
    for group_name, calculate in groups:
        try:
            features.update(calculate())
        except Exception as exc:
            raise RuntimeError(f"Failed to compute {group_name}: {exc}") from exc
    return features


def _metadata(sample_path: Path, summary_path: Path, dimension: int) -> Dict[str, Any]:
    sample_match = SAMPLE_NAME_PATTERN.search(sample_path.stem)
    summary_match = SUMMARY_NAME_PATTERN.search(summary_path.stem)
    problem_match = re.search(r"problem_(\d+)", str(summary_path.parent), re.IGNORECASE)

    if summary_match:
        summary_dimension = int(summary_match.group("dimension"))
        if summary_dimension != dimension:
            raise ValueError(
                f"Sample dimension {dimension} does not match {summary_path.name} "
                f"(dimension {summary_dimension})."
            )
    return {
        "fid": int(problem_match.group(1)) if problem_match else None,
        "dim": dimension,
        "seed_number": int(summary_match.group("seed")) if summary_match else (
            int(sample_match.group("seed")) if sample_match else None
        ),
        "sample_file": str(sample_path.resolve()),
        "evaluation_file": str(summary_path.resolve()),
    }


def _matching_pairs(
    sample_ids: Iterable[int],
    points: np.ndarray,
    evaluations: Dict[int, float],
    source: Path,
) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    point_by_id = dict(zip(sample_ids, points))
    common_ids = sorted(point_by_id.keys() & evaluations.keys())
    if len(common_ids) < 2:
        raise ValueError(
            f"Only {len(common_ids)} matching point/evaluation IDs in {source}; "
            "at least two are required."
        )
    unmatched_points = sorted(point_by_id.keys() - evaluations.keys())
    unmatched_evaluations = sorted(evaluations.keys() - point_by_id.keys())
    if unmatched_points or unmatched_evaluations:
        print(
            f"Warning: {source.name}: using {len(common_ids)} matching IDs; "
            f"{len(unmatched_points)} sample IDs and {len(unmatched_evaluations)} "
            "evaluation IDs have no counterpart.",
            flush=True,
        )
    X = np.asarray([point_by_id[simulation_id] for simulation_id in common_ids], dtype=float)
    y = np.asarray([evaluations[simulation_id] for simulation_id in common_ids], dtype=float)
    return X, y, common_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=DEFAULT_SAMPLE_FILE,
                        help="Sample CSV/TSV with design variables and optional id/sim_id column.")
    parser.add_argument("--evaluations", type=Path, default=DEFAULT_EVALUATION_DIR,
                        help="MECHBench summary file or directory containing summary files.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="Output CSV path.")
    parser.add_argument("--normalization", choices=("both", "raw", "minmax"), default="both",
                        help="Compute features from raw objectives, min-max values, or both.")
    parser.add_argument("--seed", type=int, default=None,
                        help="Seed for pflacco information-content features (defaults to sample seed or 43).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_ids, points = load_samples(args.samples)
    summaries = discover_summaries(args.evaluations)
    sample_match = SAMPLE_NAME_PATTERN.search(args.samples.stem)
    if args.evaluations.is_dir():
        matching_summaries = []
        for summary_path in summaries:
            summary_match = SUMMARY_NAME_PATTERN.search(summary_path.stem)
            if summary_match:
                if int(summary_match.group("dimension")) != points.shape[1]:
                    continue
                if (
                    sample_match
                    and int(summary_match.group("seed")) != int(sample_match.group("seed"))
                ):
                    continue
            matching_summaries.append(summary_path)
        summaries = matching_summaries
        if not summaries:
            raise FileNotFoundError(
                f"No evaluation summaries match dimension {points.shape[1]}"
                + (
                    f" and seed {sample_match.group('seed')}"
                    if sample_match
                    else ""
                )
                + f" under {args.evaluations}."
            )
    modes = ("raw", "minmax") if args.normalization == "both" else (args.normalization,)
    rows: List[Dict[str, Any]] = []

    for summary_path in summaries:
        metadata = _metadata(args.samples, summary_path, points.shape[1])
        if args.seed is not None:
            metadata["seed_number"] = args.seed
        seed = metadata["seed_number"] if metadata["seed_number"] is not None else 43
        evaluations = load_evaluations(summary_path)
        X, y, ids = _matching_pairs(sample_ids, points, evaluations, summary_path)
        for mode in modes:
            features = compute_ela_features(X, y, seed=seed, normalize=mode == "minmax")
            rows.append({
                **metadata,
                "n_samples": len(ids),
                "normalization": mode,
                **features,
            })
        print(
            f"Computed ELA for {summary_path} using {len(ids)} matched points.",
            flush=True,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"Saved {len(rows)} ELA rows to {args.output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
