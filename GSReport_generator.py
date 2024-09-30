from src import sob
import numpy as np
import pandas as pd
import platform
#import glob
import os
import multiprocessing
import json


# ------------------------------------------------------------------
#------------------ CONSTANTS FOR DIRECTORIES ----------------------
# ------------------------------------------------------------------

CURRENT_DIRECTORY:str = os.getcwd()
INPUT_DIRECTORY:str = "INPUT"
OUTPUT_DIRECTORY:str = "OUTPUT"
JSON_FILE_PATH:str= "problem.json"



# Set the problem dimension here
PROBLEM_DIMENSION:int = 2

# Number of Samples
NUM_SAMPLES:int = 10

# Set the output data type for the problem
POSSIBLE_OUTPUT_DATA:tuple = ("mass","absorbed_energy","intrusion")
OUTPUT_DATA:str = POSSIBLE_OUTPUT_DATA[2]


def mkdir(new_dir_name:str)->None:
    joined_path = os.path.join(CURRENT_DIRECTORY,new_dir_name)
    if not os.path.exists(joined_path):
        # Create the directory
        os.mkdir(joined_path)
    else:
        import shutil

        # Delete everything 
        try:
            shutil.rmtree(path=joined_path,ignore_errors=True)
        except shutil.ExecError as e:
            # Print the error
            print(e.args)
        
        # Regenerate the directory
        os.mkdir(joined_path)

def generate_samples(num_samples:int,prob_def:dict)->dict:

    
    # import the scipy libraries
    from scipy.stats.qmc import LatinHypercube 
    from scipy.stats.qmc import Sobol
    from SALib.sample.morris import sample

    # LHS Sampler
    lhs_sampler = LatinHypercube(d=prob_def['num_vars'],scramble=True,
                                optimization="lloyd")
    
    # Sobol Sampler
    sobol_sampler = Sobol(d=prob_def['num_vars'],
                          scramble=True,
                          optimization="random-cd"
                          )
    
    # Morris samples (These samples are just between the range of -5 to 5)
    morris_samples = sample(prob_def,N=num_samples,num_levels=6)


    # Generate the LHS Samples
    lhs_samples = lhs_sampler.random(n=num_samples)

    # Get the Sobol Samples
    sobol_samples = sobol_sampler.random_base2(int(np.ceil(np.log2(num_samples))))

    # Rescale the samples to -5 to 5
    lhs_samples = -5 + (10*lhs_samples)
    sobol_samples = -5 + (10*sobol_samples)

    # Generate a directory with the samples
    return {"sobol":sobol_samples,
            "morris":morris_samples,
            "lhs":lhs_samples}
    



def set_up_problem(problem_type:int=1,
                   dimension:int=2,
                   output_data:str="intrusion")->(sob.problems.CrashTube|sob.problems.StarBox|sob.problems.ThreePointBending):
    linux_system = not platform.system()=="Windows"
    if linux_system:
        #batch_file_path = "/media/feifan/TSHIBA/12_GitHub/OpenRadioss/linux_scripts_mk3/openradioss_run_script_ps.sh"
        batch_file_path = "/home/ivanolar/Documentos/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui/runopenradioss.py"
        
    else:
        #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
        batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
    
    return sob.get_problem(problem_type,dimension,output_data,
                           batch_file_path=batch_file_path)


## TODO: This function shall be modified further to get the names of the variables
def generate_json(problem:(sob.problems.CrashTube|sob.problems.StarBox|sob.problems.ThreePointBending))->None:
    
    if isinstance(problem,sob.problems.CrashTube):
        
        prob_dict = {
        'num_vars': problem.dimension,
        'names': ['X'+str(x) for x in range(problem.dimension)],
        'bounds': [[-5.0, 5.0]] * problem.dimension
        }
    elif isinstance(problem,sob.problems.StarBox):
        prob_dict = {
        'num_vars': problem.dimension,
        'names': ['X'+str(x) for x in range(problem.dimension)],
        'bounds': [[-5.0, 5.0]] * problem.dimension
        }
    else:
        prob_dict = {
        'num_vars': problem.dimension,
        'names': ['X'+str(x) for x in range(problem.dimension)],
        'bounds': [[-5.0, 5.0]] * problem.dimension
        }


    # Generate the json file
    try:
        fp_obj = open(os.path.join(CURRENT_DIRECTORY,JSON_FILE_PATH),mode="w+",encoding="utf8")
    except Exception as e:
        raise("something wrong happened!",
              e.args)
        
    dumper = json.dump(obj=prob_dict,fp=fp_obj,indent=6)
    fp_obj.close()

    return prob_dict


"""
Once the optimization problem instance has been generate, 
the model is determined (mesh and fem data loaded) only when the variable array has been input.
"""

def main():
    
    # Create the directories
    mkdir(INPUT_DIRECTORY) # Input
    mkdir(OUTPUT_DIRECTORY) # Output

    # Set the problem
    a = set_up_problem(1,dimension=PROBLEM_DIMENSION,output_data=OUTPUT_DATA)

    # Write the JSON File defining the problem for the GSA Report
    prob_def:dict = generate_json(a)

    # Generate the samples
    sample_dict:dict = generate_samples(num_samples=NUM_SAMPLES,prob_def=prob_def)


    # Run the experiments
    y_sobol = []
    y_morris = []
    y_lhs = []

    # Loop for all the sampling types
    for key,samples in sample_dict.items():
        # Write the corresponding files
        np.savetxt(os.path.join(CURRENT_DIRECTORY,INPUT_DIRECTORY,f"x_{key}"),np.array(samples))

        # Loop for all the samples
        for _, iSampl in enumerate(samples):
            # Feed the object definition
            res = a(iSampl)

            if key == "sobol":
                y_sobol.append(res)
            elif key== "morris":
                y_morris.append(res)
            elif key =="lhs":
                y_lhs.append(res)
            else:
                raise Exception("Something wrong happened")
            
        # Write the file
        if key == "sobol":
            np.savetxt(os.path.join(CURRENT_DIRECTORY,OUTPUT_DIRECTORY,f"y_{key}"),np.array(y_sobol))
        elif key== "morris":
            np.savetxt(os.path.join(CURRENT_DIRECTORY,OUTPUT_DIRECTORY,f"y_{key}"),np.array(y_morris))
        elif key =="lhs":
            np.savetxt(os.path.join(CURRENT_DIRECTORY,OUTPUT_DIRECTORY,f"y_{key}"),np.array(y_lhs))
        else:
            raise Exception("Something wrong happened")


if __name__ == '__main__':
    main()
