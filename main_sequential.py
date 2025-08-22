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
Once the optimization problem instance has been generate, 
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
    parser.add_argument("--np", type=int, default=8,
                        help="Number of processes to run the simulation.")
    parser.add_argument("--nt", type=int, default=1,
                        help="Number of threads to use in the simulation.")
    parser.add_argument("--h_level", type=int, default=1,
                        help="Hierarchy level for the simulation.")
    parser.add_argument("--write_vtk", action='store_true',
                        help="Flag to write VTK files during the simulation.")
    
    parser.add_argument("--gmsh_verbosity", type=int, default=0,
                        help="Verbosity level for GMSH operations.")
    
    parser.add_argument("--sim_id", type=int, required=True,
                        help="Simulation ID to define the results folder name.")
    parser.add_argument("--vector", nargs='+', type=float, required=True,
                        help="Vector where the objective function is evaluated. It should have as many components as the second input argument in get_problem.")
    
    parser.add_argument("--problem_type", type=int, choices=[1,2,3], default=3,
                        help="Type of problem to solve. Pick a number from 1 to 3.")
    
    parser.add_argument("--objective_function", type=str, default="load_uniformity",
                        help="Objective function to be used in the simulation. Default is 'load_uniformity'.")

    # Add more arguments as needed for your simulation


    return parser.parse_args()

def main():

    ### NOTE: This is an addition to the original code to save the results in a specific folder

    #os.chdir("/scratchdata/olarterodriguezi/") # Change to the directory where the results will be saved
    sim_id = 257 # Attribute to define the simulation id and connected results folder name
    vector = [2.5,-1] # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    f = sob.get_problem(3,2,runnerOptions,"intrusion",sequential_id_numbering=False)
    # Start timing
    start_time = time.perf_counter()
    obj_value = f(vector,sim_id)
    end_time = time.perf_counter()
    print(obj_value)
    print(f"Time taken: {end_time - start_time} seconds")

if __name__ == '__main__':
    args = parse_args()

    # Generate the scratchdata directory if it does not exist
    scratchdata_dir = f"/scratchdata/olarterodriguezi/{args.sim_id}/"
    if not os.path.exists(scratchdata_dir):
        os.makedirs(scratchdata_dir)
    # Move the main directory to the scratchdata directory
    os.chdir(f"/scratchdata/olarterodriguezi/{args.sim_id}/")  # Change to the directory where the results will be saved

    # Update runnerOptions with parsed arguments
    runnerOptions.update({
        "open_radioss_main_path": args.orss_main_path,
        "np": args.np,
        "nt": args.nt,
        "h_level": args.h_level,
        "write_vtk": args.write_vtk,
        "gmsh_verbosity": args.gmsh_verbosity
    })

    sim_id = args.sim_id
    vector = args.vector

    # check if the vector has a range between -5 and 5
    if not all(-5.00 <= v <= 5.00 for v in vector):
        raise ValueError("All components of the vector must be between -5.00 and 5.00.")

    dim = len(vector)  # Dimension is determined by the length of the vector
    problem_type = args.problem_type
    objective_function = args.objective_function

    # Get the problem instance
    f = sob.get_problem(problem_type, dim, runnerOptions, objective_function, sequential_id_numbering=False)

    start_time = time.perf_counter()
    obj_value = f(vector,sim_id)
    end_time = time.perf_counter()
    
    print(f"Time taken: {end_time - start_time} seconds")

    print(obj_value)
    

    
