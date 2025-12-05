import numpy as np
import pandas as pd

# Import the platform
import platform

from src.sob.sampler import Sampler
from src.sob.observer import Observer


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
                 "np":4,
                 "nt":1,
                 "h_level":1,
                 "gmsh_verbosity":0,
}

def main():
    sampler_instance = Sampler(model_number=2,
                                   dimension=5,
                                   runner_options=runnerOptions,
                                   output_data=["intrusion","mass","mean_impact_force","max_impact_force","load_uniformity"],
                                   root_folder=None,
                                   )
    
    observer_instance = Observer(root="results/run1",folder_name="test_run")

    sampler_instance.attach_observer(observer_instance)

    for i in range(10):
        vector = np.random.uniform(-5,5,(5,)).tolist()  # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
        print(f"Evaluating vector: {vector}")
        obj_value = sampler_instance(vector)
        print(obj_value)
    
    sampler_instance.detach_observer()

    

if __name__ == '__main__':
    main()
