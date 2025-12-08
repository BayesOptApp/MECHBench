import numpy as np
import pandas as pd

# Import the platform
import platform

from src.sob.problems import IOH_Single_Objective
import ioh.iohcpp.logger as IOHLogger


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

def mult1(x:list, *output_data,**kwargs)->float:
    r"""
    This is a function definition to evaluate the SEA from the absorbed energy and mass values.

    Args:
        x (list): The input variable list (not used in this function but required for compatibility).
        output_data (tuple): A tuple containing the output data from the physical model.
        kwargs: Additional keyword arguments that should include 'absorbed_energy' and 'mass' keys. The values
                corresponding to these keys indicate the indices in output_data where the absorbed energy and mass values are located.
    """
    # Get the keys and values from kwargs

    idx_1 = kwargs.get('absorbed_energy',None)
    idx_2 = kwargs.get('mass',None)

    if idx_1 is None:
        raise Exception("The absorbed energy value must be provided in the kwargs.")
    
    if idx_2 is None:
        raise Exception("The mass value must be provided in the kwargs.")
    
    return -1*output_data[idx_1]/output_data[idx_2]

def main():

    # Generate a logger
    triggers:list = [IOHLogger.trigger.ALWAYS]
    additional_properties = [IOHLogger.property.RAWYBEST]

    logger = IOHLogger.Analyzer( triggers=triggers,
                                 additional_properties=additional_properties,
                                 root="results/run2",
                                 folder_name="test_run",
                                 algorithm_name="Random_Search",
                                 store_positions=True
                                 )
    
    # Create the problem instance
    problem_instance = IOH_Single_Objective(model_number=2,
                                                    dimension=5,
                                                    runner_options=runnerOptions,
                                                    output_data_labels=["mass",'absorbed_energy'],
                                                    problem_name="Test_Problem",
                                                    functional_definition=mult1,
                                                    root_folder=logger.output_directory
                                   )
    
    # Attach the logger to the problem instance
    problem_instance.attach_logger(logger)

    for i in range(50):
        vector = np.random.uniform(-5,5,(5,)).tolist()  # Vector where the objective function is evaluated, it has as many components as the second input argument in get_problem below
        print(f"Evaluating vector: {vector}")
        obj_value = problem_instance(vector)
        print(obj_value)
    
    #sampler_instance.detach_observer()

    # Detach the logger from the problem instance
    problem_instance.reset()
    problem_instance.detach_logger()

    

if __name__ == '__main__':
    main()
