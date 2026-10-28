"""Optimize an IOHexperimenter BBOB function using the existing wrapopt algorithms.

Example:
    python main_experiments.py --optimizer pyCMA --function_id 1 \
        --instance 3 --dimension 5 --run_id 0 --max_iter 100

Requires ioh and the same ML_Benchmarks/optimizer dependencies as the original
script. OpenRadioss, MECHBench, and src.sob are no longer needed.
"""

# GENERAL LIBRARIES
import argparse
import os
import platform
from pathlib import Path

# ---- CONSTANTS ----
OPTIMIZERS = ['pyCMA', 'BO_botorch', 'BAxUS_botorch', 'HEBO', 'turbo1']


# THIS IS A HELPER FUNCTION TO PARSE ARGUMENTS. IT IS NOT PART OF THE MAIN LOGIC.
def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Optimize a BBOB function with IOHexperimenter and algorithms from the HDBO Benchmarking Study from Santoni et al. (2023).')
    parser.add_argument('--optimizer', choices=OPTIMIZERS, default='pyCMA',
                        help='Optimizer to use (default: pyCMA).')
    parser.add_argument('--function_id', '--problem_type', dest='function_id',
                        type=int, choices=[*range(1, 25)], default=1,
                        help='BBOB function ID, from 1 to 24 (default: 1, Sphere).')
    parser.add_argument('--instance', '--instance_id', dest='instance',
                        type=positive_int, default=1,
                        help='BBOB instance ID (default: 1), independent of run_id.')
    parser.add_argument('--run_id', type=int, required=True,
                        help='Independent optimizer run ID; also contributes to the seed.')
    parser.add_argument('--dimension', type=int, choices=[2, 5, 10, 20], default=5,
                        help='Problem dimension (default: 5).')
    parser.add_argument('--max_iter', type=positive_int, default=100,
                        help='Total evaluation budget passed to wrapopt (default: 100).')
    parser.add_argument('--doe_mult', type=positive_int, default=3,
                        help='DoE size = doe_mult * dimension (default: 3).')
    parser.add_argument('--seed', type=int, default=42,
                        help='Base seed; optimizer seed = seed + run_id (default: 42).')
    parser.add_argument('--ioh_logger_folder', type=str, default=None,
                        help='IOH log root (default: logdata_optimization beside this script).')
    
    args = parser.parse_args()
    if args.seed + args.run_id < 0:
        parser.error('seed + run_id must be nonnegative')
    if args.doe_mult * args.dimension > args.max_iter:
        parser.error('DoE size (doe_mult * dimension) must not exceed max_iter')
    return args


def main() -> None:
    args = parse_args()

    import ioh
    try:
        from ML_Benchmarks import wrapopt
    except ImportError as e:
        raise ImportError("Could not import wrapopt from ML_Benchmarks. Ensure that the ML_Benchmarks directory is in the PYTHONPATH.") from e


    # Resolve user paths before changing the working directory.
    original_cwd = Path.cwd()

    if args.ioh_logger_folder is not None:
        ioh_logger_folder = Path(args.ioh_logger_folder).expanduser().resolve()
    else:
        ioh_logger_folder = original_cwd.resolve() / 'logdata_optimization'
    ioh_logger_folder.mkdir(parents=True, exist_ok=True)

    # Use the native BBOB problem so IOH retains its function/instance metadata,
    # bounds, known optimum, and evaluation counter.
    ioh_problem = ioh.get_problem(
        fid=args.function_id,
        instance=args.instance,
        dimension=args.dimension,
        problem_class=ioh.ProblemClass.BBOB,
    )
    doe_size = args.doe_mult * args.dimension
    random_seed = args.seed + args.run_id

    logger = ioh.logger.Analyzer(
        root=ioh_logger_folder.resolve().as_posix(),
        triggers=[ioh.logger.trigger.ALWAYS],
        additional_properties=[ioh.logger.property.RAWYBEST],
        folder_name=(f'BBOB_f{args.function_id}_instance_{args.instance}'
                     f'_d{args.dimension}_run_{args.run_id}_{args.optimizer}'),
        algorithm_name=args.optimizer,
        algorithm_info=f'DOE_SIZE:{doe_size};SEED:{random_seed}',
        store_positions=True,
    )
    try:
        ioh_problem.attach_logger(logger)
        # Preserve the optimizer interface, DoE rule, and total-budget logic.
        optimizer = wrapopt(
            args.optimizer, ioh_problem,
            ml_total_budget=args.max_iter,
            ml_dim=args.dimension,
            ml_DoE_size=doe_size,
            random_seed=random_seed,
        )
        optimizer.run()

    finally:
        # Flush the run even when optimization fails, then release the logger.
        try:
            ioh_problem.reset()
        finally:
            try:
                ioh_problem.detach_logger()
            finally:
                logger.close()
                os.chdir(original_cwd)


if __name__ == '__main__':
    main()
