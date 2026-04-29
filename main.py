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
    orss_main_path = "/home/ivanolar/Documents/MECHBench_v2/MECHBench/OpenRadioss_linux64/OpenRadioss"
    
else:
    orss_main_path = "C:/Users/iolar/Documents/OpenRadioss/OpenRadioss"

#runnerOptions = {"open_radioss_main_path":orss_main_path}

runnerOptions = {"np":1}

def main():
    sim_id = 236 # Attribute to define the simulation id and connected results folder name
    dim = 5
    vector = np.random.uniform(-5,5,(dim,)).tolist()  # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    #vector = [np.zeros((20,)).tolist()] # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
    #vector = [3.129063083151652, 1.276036348453766, 4.649091432673452, -1.0578641202273658, -2.223620398000797, 4.972163051323311, -1.5671976458618349, -2.2839067338083954, -0.1137787419587058, 0.7028705602294962, -1.935981339832131, -4.986845327202568, -4.979443528081319, 4.248263867987974, 2.1888934443123524, -3.850223213447326, 3.6327155165559946, 2.3588020285909357, 2.9293771388783707, -1.263669089834881]
    print(f"Evaluating vector: {vector}")
    f = sob.get_problem(2,dim,runnerOptions,["load_uniformity","intrusion"])
    obj_value = f(vector,sim_id)
    print(obj_value)
    

if __name__ == '__main__':
    main()
