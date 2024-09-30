from src import sob
import numpy as np
import pandas as pd


'''
Added the following block to identify if the current system 
'''

'''
Once the optimization problem instance has been generate, the model is determined (mesh and fem data loaded) only when the variable array has been input.
'''

def main():
    linux_system = True
    if linux_system:
        #batch_file_path = "/media/feifan/TSHIBA/12_GitHub/OpenRadioss/linux_scripts_mk3/openradioss_run_script_ps.sh"
        batch_file_path = "/home/ivanolar/Documentos/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui/runopenradioss.py"
        
    else:
        #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
        batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"

    a = sob.get_problem(1,3,'intrusion',batch_file_path)
    intrusion1 = a([-2,2,4.1])
    print(intrusion1)
    


if __name__ == '__main__':
    main()
