from src import sob
import numpy as np
import pandas as pd



'''
Added the following block to identify if the current system 
'''

#runnerOptions = {"open_radioss_main_path":"/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/",
#                 "write_vtk":False,
#                 "np":4, # Number of processes to run the simulation
#                 "nt":1,
#                 "h_level":1,
#                 "gmsh_verbosity":0,
#}   

runnerOptions = {"np":1, # Number of processes to run the simulation
                 "nt":1,
                 "h_level":1,
                 "gmsh_verbosity":0,
}

r'''
Once the optimization problem instance has been generate, 
the model is determined (mesh and fem data loaded) only when the variable array has been input.
'''


def main():
    sim_id = 237 # Attribute to define the simulation id and connected results folder name
    dim = 5#vector = [np.zeros((20,)).tolist()] # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    problem_id = 3 # 1: star box, 2: three point bending, 3: crash tube

    vector = np.random.uniform(0,0,(dim,)).tolist()  # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    print(f"Evaluating vector: {vector}")
    f = sob.get_problem(problem_id,dim, runner_options=runnerOptions)
    obj_value = f(vector,sim_id)
    print(obj_value)
    

if __name__ == '__main__':
    main()
