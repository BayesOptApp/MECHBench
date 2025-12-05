from src import sob
import numpy as np
import pandas as pd

# Import the platform
import platform


'''
Added the following block to identify if the current system 
'''

r'''
Once the optimization problem instance has been generate, 
the model is determined (mesh and fem data loaded) only when the variable array has been input.
'''

linux_system = not (platform.system() == 'Windows')

if linux_system:
    orss_main_path = "/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/"
    
else:
    orss_main_path = "C:/Users/iolar/Documents/OpenRadioss/OpenRadioss"

runnerOptions = {"open_radioss_main_path":orss_main_path,
                 "write_vtk":True,
                 "np":8,
                 "nt":1,
                 "h_level":1,
                 "gmsh_verbosity":0,
}

def main():
    sim_id = 275 # Attribute to define the simulation id and connected results folder name
    vector = np.random.uniform(-5,5,(5,)).tolist()  # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    # vector = np.zeros((30,)).tolist()] # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    print(f"Evaluating vector: {vector}")
    f = sob.get_problem(2,5,runnerOptions,["intrusion","mass","mean_impact_force","max_impact_force","load_uniformity"],sequential_id_numbering=False)
    obj_value = f(vector,sim_id)
    print(obj_value)
    

if __name__ == '__main__':
    main()
