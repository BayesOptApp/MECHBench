import subprocess

def run_radioss(input_file_path, batch_file_path, isShell = False, write_vtk=False, write_csv = True, runStarter = False):
    
    nt = "1"
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

    # Construct the command to run the batch file with arguments
    command = [
        batch_file_path,
        input_file_path, # %1
        nt, # %2
        np, # %3
        sp, # %4
        vtk_option, # %5
        csv_option, # %6
        starter_option # %7
    ]

    try:
        # Run the batch file and wait for it to complete
        subprocess.run(command, check=True, shell=isShell)
        print("Batch file executed successfully.")

        # Check if there was an error
    except subprocess.CalledProcessError as e:
        print(f"Error running batch file: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")