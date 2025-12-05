from src.sob.physical_models.fem_settings.abstractFEMSettings import AbstractFEMSettings
from src.sob.physical_models.meshes import CrashTubeMesh
from src.sob.physical_models.fem_settings.starBoxModel import StarBoxModel
import numpy as np

class CrashTubeModel(StarBoxModel):
    def __init__(self, mesh: CrashTubeMesh, **kwargs) -> None:
        self.mesh = mesh
        
        mesh.units =  '  kg  mm  ms  kN  GPa  kN-mm'
        self.units = mesh.units
        mesh.write_mesh_file()
        # Define default parameters for load_impactor
        self.impactor_defaults = {
            'wall_n_id': 999999,
            'wall_loc': self.mesh.extrusion_length+1,
            'wall_mass': 300.0,
            'wall_vel': 8.33,
        }
        self.rigid_mass = self.impactor_defaults['wall_mass']
        
        # Define default parameters for load_database
        self.database_defaults = {
            'end_time': 50.0,
            'database_dtime': 0.05,
            'consider_d3plot': True, # cons_d3plot
            'd3plot_dtime': 0.5, # d3plot_dtime
            'cons_d3thdt' : False,
            'd3thdt_dtime': 10000, # d3thdt_dt
            'shell_warping': 1,
            'binary_ascii': 2,
            'write_cshell':True,
            'tdel': 0.0
            # Add other default parameters here
        }
        
        self.shell_warping = 1    # BWC, lsdyna default is 2. if there is warping set it to 1
        self.binary_ascii = 2  # 1: only ascii   2: only binary   3: both ascii and binary
        self.write_cshell = True   
        self.tdel = 0.0

        self.intout = 'STRESS'
        self.nodout = 'STRESS' 

        self.intp_db = 3
        self.sigflg = 1
        self.epsflg = 1
        self.rltflg = 1

        # Define default parameters for load_material
        self.material_defaults = {
            'mat_id': 999,
            'mat_density': 7.83E-6,
            'mat_young_mod': 200.0,
            'mat_poisson_r': 0.3,
            'mat_yield_initial': 0.366,
            'mat_tang_mod': 0.0,
            'mat_failure_pstrain': 1.0E+21,
            'mat_cowper_symond_c': 40.0,
            'mat_cowper_symond_p': 5.0,
            'mat_vp_rate_efffect': 1,
            'mat_load_curve_id': 1,
            'mat_effective_plastic_strain_stress': np.array([[0., 0.366], [2.5e-2, 0.4240],
                                                             [4.9e-2, 0.476], [7.2e-2, 0.507],
                                                             [9.5e-2, 0.529], [0.118, 0.546],
                                                             [0.140, 0.559], [0.182, 0.584]]),
            # Add other default parameters here
        }

        self._load_impactor(**kwargs)
        self._load_database(**kwargs)
        self._load_material(**kwargs)

        # ------- nodal forces
        self.write_nod_force_top = False
        self.write_nod_force_bottom = True