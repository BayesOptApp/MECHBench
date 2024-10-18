import subprocess
import os


def run_radioss(input_file_path, 
                batch_file_path, 
                isShell = False, 
                write_vtk= False, 
                write_csv = True, 
                runStarter = False,
                d3plot=False):
    
    # TODO: This is an input line to get the number ofcores available to run
    
    #if os.cpu_count() > 1: nt = str(int(os.cpu_count()/2))
    if os.cpu_count() > 1: nt = str(8)
    else: nt="1"

    # Get the directory of the folder storing the deck files
    curdir = os.path.dirname(input_file_path)

    np = "1"
    sp = "sp"  # "sp" or any other value for non-sp

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
    command = [
        "python3", # Additional setup to call the Python3 Setup
        batch_file_path,
        input_file_path, # %1
        nt, # %2
        np, # %3
        sp, # %4
        vtk_option, # %5
        csv_option, # %6
        starter_option, # %7,
        d3_plot_option # TODO: This is for the 3d plot option which is not considered
    ]

    # Create input output Wrappers
    out_out = open(os.path.join(curdir,"output.txt"),
                   mode="w+",encoding="utf-8")
    
    error_out = open(os.path.join(curdir,"error.txt"),
                   mode="w+",encoding="utf-8")
    try:
        # Run the batch file and wait for it to complete
        subprocess.run(command, check=True, shell=isShell,
                       stderr=error_out,
                       stdout=out_out,
                       encoding="utf-8")
        print("Batch file executed successfully.")

        # Check if there was an error
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError(f"Error running batch file: {e}")
    except Exception as e:
        raise Exception(f"An error occurred: {e}")
    else:
        # Close the I/O connectors
        out_out.close()
        error_out.close()