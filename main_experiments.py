# Import the wrapopt function
from ML_Benchmarks import wrapopt

# Import IOH
import ioh

# Install dependencies to control MECHBENCH
from src import sob
import numpy as np
import pandas as pd
import time

# Import the platform
import platform

# Import system libraries
import os, sys

# Import argparse
import argparse


'''
Added the following block to identify if the current system 
'''

r'''
Once the optimization problem instance has been generated, 
the model is determined (mesh and fem data loaded) only when the variable array has been input.
'''

linux_system = not (platform.system() == 'Windows')

if linux_system:
    #orss_main_path = "/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/"
    orss_main_path = "/home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/"
    
else:
    orss_main_path = "C:/Users/iolar/Documents/OpenRadioss/OpenRadioss"

runnerOptions = {"open_radioss_main_path":orss_main_path,
                 "write_vtk":False,
                 "np":32, # Number of processes to run the simulation
                 "nt":1,
                 "h_level":1,
                 "gmsh_verbosity":0,
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parse arguments for the simulation.")
    parser.add_argument("--orss_main_path", type=str, default=orss_main_path,
                        help="Path to the OpenRadioss main directory.")
    parser.add_argument("--np", type=int, default=4,
                        help="Number of processes to run the simulation.")
    parser.add_argument("--nt", type=int, default=1,
                        help="Number of threads to use in the simulation.")
    parser.add_argument("--h_level", type=int, default=1,
                        help="Hierarchy level for the simulation.")
    
    parser.add_argument("--write_vtk", action='store_true',
                        help="Flag to write VTK files during the simulation.")
    
    parser.add_argument("--gmsh_verbosity", type=int, default=0,
                        help="Verbosity level for GMSH operations.")


    parser.add_argument("--optimizer", type=str, choices=['pyCMA','BO_botorch','BAxUS_botorch','HEBO',"turbo1"],
                        default='pyCMA', help="Optimizer to use to optimize a defined problem.")
    parser.add_argument("--problem_type", type=int, choices=[1,2,3], default=3,
                        help="Type of problem to solve. Pick a number from 1 to 3.")
    parser.add_argument("--run_id", type=int, required=True,
                        help="Run ID for the simulation.")
    parser.add_argument("--objective_function", type=str, default="load_uniformity",
                        choices=['load_uniformity', 'intrusion', 'absorbed_energy',
                                 'mass', 'specific-energy-absorbed',"penalized_sea","penalized_mass"],
                        help="Objective function to be used in the simulation. Default is 'load_uniformity'.")
    parser.add_argument("--dimension", type=int, choices=[2,5,10,20], default=5,
                        help="Dimension of the problem. Default is 5.")
    parser.add_argument("--max_iter", type=int, default=100,
                        help="Maximum number of iterations for the optimizer. Default is 100.")
    parser.add_argument("--doe_mult", type=int, default=3,
                        help="Multiplicative factor for the Design of Experiments (DoE) size.")

    # Add more arguments as needed for your simulation


    return parser.parse_args()

# def main():

#     ### NOTE: This is an addition to the original code to save the results in a specific folder

#     #os.chdir("/scratchdata/olarterodriguezi/") # Change to the directory where the results will be saved
#     sim_id = 257 # Attribute to define the simulation id and connected results folder name
#     vector = [2.5,-1] # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
#     f = sob.get_problem(3,2,runnerOptions,"intrusion",sequential_id_numbering=False)
#     # Start timing
#     start_time = time.perf_counter()
#     obj_value = f(vector,sim_id)
#     end_time = time.perf_counter()
#     print(obj_value)
#     print(f"Time taken: {end_time - start_time} seconds")

if __name__ == '__main__':
    args = parse_args()

    # Generate the scratchdata directory if it does not exist
    scratchdata_dir = f"/scratchdata/olarterodriguezi/problem_{args.problem_type}/{args.optimizer}/{args.dimension}/{args.run_id}"
    if not os.path.exists(scratchdata_dir):
        os.makedirs(scratchdata_dir)
    # Move the main directory to the scratchdata directory
    os.chdir(scratchdata_dir)  # Change to the directory where the results will be saved
    #os.makedirs(f"/home/olarterodriguezi/MECHBench/results/optimization/MECHBENCH_problem_{args.problem_type}_{args.objective_function}/run_{args.run_id}_{args.optimizer}/ARF/", exist_ok=True)
    #os.chdir(f"/home/olarterodriguezi/MECHBench/results/optimization/MECHBENCH_problem_{args.problem_type}_{args.objective_function}/run_{args.run_id}_{args.optimizer}/ARF/")  # Change to the directory where the results will be saved    
    # Update runnerOptions with parsed arguments
    runnerOptions.update({
        "open_radioss_main_path": args.orss_main_path,
        "np": args.np,
        "nt": args.nt,
        "h_level": args.h_level,
        "write_vtk": args.write_vtk,
        "gmsh_verbosity": args.gmsh_verbosity
    })

    problem_type = args.problem_type
    objective_function = args.objective_function

    # Get the problem instance
    f = sob.get_problem(problem_type, args.dimension, runnerOptions, objective_function, sequential_id_numbering=False)

    # DoE Size
    doe_size = args.doe_mult * args.dimension


    n_iter:int = 0

    

    # Make a wrapper function

    def func_wrapper(vector)->float:
        global n_iter
        n_iter += 1
        return f(vector, n_iter)

    # Use the `get_problem` from IOH
    ioh_problem = ioh.wrap_problem(function=func_wrapper,
                                   name=f"MECHBENCH_problem_{problem_type}_{objective_function}",
                                   dimension=args.dimension,
                                   instance=0,
                                   lb=-5.0,
                                   ub=5.0)
    
    # Set up the logger
    logger = ioh.logger.Analyzer(
        root="/home/olarterodriguezi/MECHBench/results/optimization",
        triggers=[ioh.logger.trigger.ALWAYS],
        additional_properties=[ioh.logger.property.RAWYBEST], # Log the RAW_Y_BEST
        folder_name=f"MECHBENCH_problem_{problem_type}_{objective_function}/run_{args.run_id}_{args.optimizer}",
        algorithm_name=args.optimizer,
        algorithm_info=f"DOE_SIZE:{doe_size}",
        store_positions=True
    )

    ioh_problem.attach_logger(logger)

    # Use the optimizer
    optimizer = wrapopt(args.optimizer, ioh_problem, ml_total_budget=args.max_iter,
                        ml_dim=args.dimension, ml_DoE_size=doe_size, random_seed=args.run_id+42)

    # Run optimization
    optimizer.run()

    ioh_problem.reset()
    ioh_problem.detach_logger()

    logger.close()

    #start_time = time.perf_counter()
    #obj_value = f(vector,sim_id)
    #end_time = time.perf_counter()
    
    #print(f"Time taken: {end_time - start_time} seconds")

    #print(obj_value)
    

    
