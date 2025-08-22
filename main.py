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
                 "np":1,
                 "nt":4,
                 "h_level":1,
                 "gmsh_verbosity":0,
}

def main():
    sim_id = 255 # Attribute to define the simulation id and connected results folder name
    vector = [3.42834656, 2.15259375, -1.49093763, 1.69472419, -0.29146718, 1.32146257, 4.98257611, -2.80885004, -3.23567450, -4.48471372]
    # vector = np.zeros((30,)).tolist()] # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    f = sob.get_problem(3,10,runnerOptions,"load_uniformity",sequential_id_numbering=False)
    obj_value = f(vector,sim_id)
    print(obj_value)
    

if __name__ == '__main__':
    main()
