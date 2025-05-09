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
    #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
    orss_main_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"


runnerOptions = {"open_radioss_main_path":orss_main_path,
                 "write_vtk":True,
                 "np":4,
                 "nt":1,
                 "nt":1,
                 "h_level":1
}

def main():
    a = sob.get_problem(3,30,['intrusion'],runnerOptions,sequential_id_numbering=False)
    intrusion1 = a(np.hstack(([-5]*15,[5]*15)),169)
    a = sob.get_problem(3,30,['intrusion'],runnerOptions,sequential_id_numbering=False)
    intrusion1 = a(np.hstack(([-5]*15,[5]*15)),169)
    print(intrusion1)
    


if __name__ == '__main__':
    main()
