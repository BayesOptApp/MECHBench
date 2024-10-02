from src import sob
import numpy as np
import pandas as pd
import platform
#import glob
import os
import multiprocessing as mp
import json


# ------------------------------------------------------------------
#------------------ CONSTANTS FOR DIRECTORIES ----------------------
# ------------------------------------------------------------------

CURRENT_DIRECTORY:str = os.getcwd()
INPUT_DIRECTORY:str = "INPUT"
OUTPUT_DIRECTORY:str = "OUTPUT"

JSON_FILE_NAME:str= "problem.json"

#GSA_REPORT_PATH:str= "/home/ivanolar/Documentos/GSA_Report_Git/GSAreport/src/GSAreport.py"
GSA_REPORT_PATH:str= "/home/olarterodriguezi/GSAreport/src/GSAreport.py"

#OPENRADIOSS_SETUP_FILE_LINUX:str = "/home/ivanolar/Documentos/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui/runopenradioss.py"
OPENRADIOSS_SETUP_FILE_LINUX:str = "/home/olarterodriguezi/OpenRadioss_linux/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui/runopenradioss.py"

SLURM_ENVIRONMENT = "SLURM_JOB_ID" in os.environ

# DEFAULT_PROBLEM_DIMENSION
PROBLEM_DIMENSION:int = 5

# Default Number of Samples
DEFAULT_NUM_SAMPLES:int = 512

# Set the output data type for the problem
POSSIBLE_OUTPUT_DATA:tuple = ("","mass","absorbed_energy","intrusion")
DEFAULT_OUTPUT_DATA:str = POSSIBLE_OUTPUT_DATA[3]


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


def generate_samples2(num_samples:int,definition:dict)->dict:

    # Import the processing library to use the command window
    import subprocess
    
    # Construct the command to run the batch file with arguments
    command:list = [
        "python3", # Additional setup to call the Python3 Setup
        GSA_REPORT_PATH, # Path to GSA Report
        "-p",
        os.path.join(CURRENT_DIRECTORY,JSON_FILE_PATH),
        "-d",
        os.path.join(CURRENT_DIRECTORY,INPUT_DIRECTORY),
        "--sample",
        "--samplesize",
        str(num_samples)
    ]

    try:
        subprocess.run(command,check=True,shell=False)
    except subprocess.CalledProcessError as e:
        print("Error occurred with the activation", e.args)
    except subprocess.SubprocessError as e:
        print("Subprocess error", e.args)
    except Exception as e:
        print("Weird Error occured", e.args)

    
    # Generate a directory with the samples
    try:
        return {"sobol":np.loadtxt(os.path.join(CURRENT_DIRECTORY,INPUT_DIRECTORY,"x_sobol.csv"),dtype=float),
            "morris":np.loadtxt(os.path.join(CURRENT_DIRECTORY,INPUT_DIRECTORY,"x_morris.csv"),dtype=float),
            "lhs":np.loadtxt(os.path.join(CURRENT_DIRECTORY,INPUT_DIRECTORY,"x_lhs.csv"),dtype=float)}
    except FileNotFoundError as e:
        raise FileNotFoundError("One of the files wasn't found", e.args)
    except Exception as e:
        print("A weird exception has happened")



def set_up_problem(problem_type:int=1,
                   dimension:int=2,
                   output_data:str="intrusion")->(sob.problems):
    linux_system = not platform.system()=="Windows"
    if linux_system:
        batch_file_path = OPENRADIOSS_SETUP_FILE_LINUX
        
    else:
        #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
        batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
    
    return sob.get_problem(problem_type,dimension,output_data,
                           batch_file_path=batch_file_path)


## TODO: This function shall be modified further to get the names of the variables
def generate_json(problem:(sob.problems))->None:
    
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
        
    json.dump(obj=prob_dict,fp=fp_obj,indent=6)
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
    sample_dict:dict = generate_samples2(num_samples=NUM_SAMPLES,definition=prob_def)


    # Run the experiments
    y_sobol = []
    y_morris = []
    y_lhs = []

    # Loop for all the sampling types
    for key,samples in sample_dict.items():
        # Write the corresponding files
        
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
    import argparse
    import textwrap

    parser = argparse.ArgumentParser(
        prog="GSAreport_samples_adapter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """\
            Generate a global sensitivity analysis report for a given data set or function.
            Common uses cases:
            --------------------------------
            Generate samples for evaluation by a real world function / simulator
                > python GSAreport.py -p problem.json -d data_dir --sample --samplesize 1000
            Analyse the samples with their output stored in the data folder
                > python GSAreport.py -p problem.json -d data_dir -o output_dir
            Analyse a real-world data set and use a Random Forest model to interpolate (data_dir contains x.csv and y.csv)
                > python GSAreport.py -p problem.json -d data_dir -o output_dir --samplesize 10000
            """
        ),
    )
    parser.add_argument(
        "--problem",
        "-p",
        default="problem.json",
        type=str,
        help="File path to the problem definition in json format.",
    )
    parser.add_argument(
        "--data",
        "-d",
        type=str,
        default="/data/",
        help="Directory where the intermediate data is stored. (defaults to /data)",
    )  # will be accesible under args.data
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="/output/",
        help="Directory where the output report is stored. (defaults to /output/)",
    )  # will be accesible under args.output
    parser.add_argument(
        "--name",
        type=str,
        default="SA",
        help="Name of the experiment, will be used in the report output.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="The number of important features to focus on, default is 10.",
    )
    parser.add_argument(
        "--sample",
        action="store_true",
        help="When you use this flag, only the samples are generated to be used in the analyse phase.",
    )
    parser.add_argument(
        "--samplesize",
        "-n",
        type=int,
        help="Number of samples to generate.",
        default=1000,
    )
    parser.add_argument(
        "--modelsize", type=int, help="Number of samples for the model.", default=1000
    )
    parser.add_argument(
        "--demo",
        help="Demo mode, uses a BBOB function as test function.",
        action="store_true",
    )
    args = parser.parse_args()

    output_dir = args.output
    data_dir = args.data
    top = args.top
