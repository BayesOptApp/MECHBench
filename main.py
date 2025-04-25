from src import sob
import numpy as np
import pandas as pd

# Import the platform
import platform


'''
Added the following block to identify if the current system 
'''

'''
Once the optimization problem instance has been generate, the model is determined (mesh and fem data loaded) only when the variable array has been input.
'''

def main():
    linux_system = not (platform.system() == 'Windows')
    if linux_system:
        #batch_file_path = "/media/feifan/TSHIBA/12_GitHub/OpenRadioss/linux_scripts_mk3/openradioss_run_script_ps.sh"
        #batch_file_path = "/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/openradioss_gui/runopenradioss.py"
        batch_file_path = "/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/"
        
    else:
        #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
        batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"

    a = sob.get_problem(3,4,'intrusion',batch_file_path,sequential_id_numbering=True)
    intrusion1 = a([1,2,-2,1])
    print(intrusion1)
    


if __name__ == '__main__':
    main()
