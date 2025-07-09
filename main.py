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
    orss_main_path = "C:/Users/iolar/Downloads/OpenRadioss_win64_2/OpenRadioss/"

runnerOptions = {"open_radioss_main_path":orss_main_path,
                 "write_vtk":True,
                 "np":4,
                 "nt":1,
                 "h_level":1,
                 "gmsh_verbosity":0,
}

def main():
    a = sob.get_problem(1,2,"intrusion",runnerOptions,sequential_id_numbering=False)
    vector = [2.5,-1]
    intrusion1 = a(vector,253)
    print(intrusion1)
    

if __name__ == '__main__':
    main()
