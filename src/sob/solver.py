import subprocess
import os
import psutil
import platform
from pathlib import Path
from typing import Union
from .utils.run_openradioss import RunOpenRadioss


def run_radioss(input_file_path:Union[str,Path], 
                batch_file_path:Union[str,Path], 
                write_vtk:bool= True, 
                write_csv:bool = True, 
                runStarter:bool = False,
                d3plot:bool=False,
                nt_int:int = 1,
                np_int:int = 1):
    
    # TODO: This is an input line to get the number ofcores available to run
    

    nt = str(nt_int)

    # Get the directory of the folder storing the deck files
    curdir = os.path.dirname(input_file_path)

    np = str(np_int)

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
    if not runStarter:
        # Get the files which are part of the engine
        engine_list = run_env.get_engine_input_file_list(repair=False)

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

        if write_vtk:

            vtk_list = run_env.get_animation_list()

            for i_anim_file in vtk_list:
                run_env.convert_anim_to_vtk(i_anim_file)
        

