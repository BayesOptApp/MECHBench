from src import sob
import numpy as np
import pandas as pd

'''
Once the optimization problem instance has been generate, the model is determined (mesh and fem data loaded) only when the variable array has been input.
'''

def main():
    linux_system = False
    if linux_system:
        batch_file_path = "/media/feifan/TSHIBA/12_GitHub/OpenRadioss/linux_scripts_mk3/openradioss_run_script_ps_2.sh"
    else:
        batch_file_path = "C:\\Users\\iolar\\Downloads\\win_scripts_mk4\\win_scripts_mk4\\openradioss_run_script_ps.bat"

<<<<<<< HEAD
    a = sob.get_problem(2,5,'intrusion',batch_file_path)
    intrusion1 = a([1,2,3,3,1])
    print(intrusion1)
=======
    a = sob.get_problem(1,3,'intrusion',batch_file_path)
    intrusion1 = a([1,2,3])
    # print(intrusion1)
>>>>>>> 3c74c6dbaa301e307c81d45c3cf4dcdcb6d3d9ec
    


if __name__ == '__main__':
    main()
