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
        batch_file_path = "/home/ivanolar/Documentos/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui/runopenradioss.py"
        
    else:
        #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
        batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"

    a = sob.get_problem(1,5,["specific_energy",'intrusion'],batch_file_path,sequential_id_numbering=False)
    intrusion1 = a([-4.420697, -3.580725, -0.745822, -0.600424, 0.073193],4)
    print(intrusion1)
    


if __name__ == '__main__':
    main()
