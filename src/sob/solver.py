import subprocess
import os
import psutil
import platform
from pathlib import Path
from typing import Union
from .utils.run_openradioss import RunOpenRadioss




def determine_batch_file_path_type(batch_file_path:Union[str,Path])->list:
    r"""
    This function is written in order to define the command line to run the case

    Args:
    --------------
    - batch_file_path (`Union[str,Path]`): A batch_file path to set the OpenRadioss run
    """

    # Convert the file path to a path
    # NOTE: This handle is meant for UNIX
    if isinstance(batch_file_path,str):
        batch_file_path = Path(batch_file_path)
    


def run_radioss(input_file_path:Union[str,Path], 
                batch_file_path:Union[str,Path], 
                write_vtk:bool= True, 
                write_csv:bool = True, 
                runStarter:bool = False,
                d3plot:bool=False):
    
    # TODO: This is an input line to get the number ofcores available to run
    
    #if os.cpu_count() > 1: nt = str(int(os.cpu_count()/2))
    #if os.cpu_count() > 1: nt = str(8)
    #else: nt="1"

    nt = "1"

    # Get the directory of the folder storing the deck files
    curdir = os.path.dirname(input_file_path)

    np = "4"

    #######
    # ------------- WARNING ----------------------
    # Simulations 
    #######
    sp = "dp"  # "sp" or any other value for non-sp 

    if write_vtk:
        # "yes" or "no" depending on whether you want to convert Anim files to vtk (for ParaView)
        vtk_option = "yes"  
    else:
        vtk_option = "no"
    
    if write_csv:
        csv_option = "yes"  # "yes" or "no" depending on whether you want to convert TH files to csv
    else:
        csv_option = "no"
    
    if runStarter:
        starter_option = "yes"  # "no" or "yes" depending on whether you are running the starter
    else:
        starter_option = "no"

    if d3plot:
        d3_plot_option = "yes"
    else:
        d3_plot_option = "no"

    # Construct the command to run the batch file with arguments
    # command = [
    #     "python3", # Additional setup to call the Python3 Setup
    #     batch_file_path,
    #     input_file_path, # %1
    #     nt, # %2
    #     np, # %3
    #     sp, # %4
    #     vtk_option, # %5
    #     csv_option, # %6
    #     starter_option, # %7,
    #     d3_plot_option # TODO: This is for the 3d plot option which is not considered
    # ]

    
    command = [
        input_file_path, # %1
        nt, # %2
        np, # %3
        sp, # %4
        vtk_option, # %5
        csv_option, # %6
        starter_option, # %7,
        d3_plot_option # TODO: This is for the 3d plot option which is not considered
    ]


    # Activate the OpenRadiossRunningEnvironment

    run_env:RunOpenRadioss = RunOpenRadioss(command,
                                            0,
                                            batch_file_path) 
    
    # Clean the environment
    run_env.delete_previous_results()

    # Get run command
    starter_command = run_env.get_starter_command()

    # Run the starter
    output_starter = subprocess.run(args=starter_command,env=run_env.environment(),
                                    cwd=run_env.running_directory,
                                    #shell=isShell, 
                                    shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
    
    if output_starter.stderr is not None:
        raise Exception("Starter not working: \n"+\
                        output_starter.stderr)

    # ======================================
    # run the engine
    # ======================================
    if starter_option == 'no':
        # Get the files which are part of the engine
        engine_list = run_env.get_engine_input_file_list()

        engine_stdout:list = []

        for iFile in engine_list:
            # Get the command
            iCommand = run_env.get_engine_command(iFile,True)
            engine_stdout.append(subprocess.run(args=iCommand,#env=run_env.environment(),
                                    cwd=run_env.running_directory,
                                    #shell=isShell, 
                                    shell = False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT))




    # ======================================
    # Convert to TH
    # ======================================

    th_list = run_env.get_th_list()

    for i_th_file in th_list:
        run_env.convert_th_to_csv(i_th_file)

    # ======================================
    # Convert to VTK
    # ======================================

    if vtk_option:

        vtk_list = run_env.get_animation_list()

        for i_anim_file in vtk_list:
            run_env.convert_anim_to_vtk(i_anim_file)
        

