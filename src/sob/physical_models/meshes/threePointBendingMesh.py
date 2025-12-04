from src.sob.physical_models.meshes.abstractMeshSettings import AbstractMeshSettings
import numpy as np
import json
from src.sob.physical_models.meshes.routines.gmsh.three_point_bending_gmsh import ThreePointBending_GMSH
from typing import List
from copy import deepcopy


class ThreePointBendingMesh(AbstractMeshSettings):
    def __init__(self, variable_array, h_level:int=1, 
                 gmsh_verbosity:bool=False,
                 **kwargs) -> None:
         # Set the h_level
        self.h_level = h_level

        self.gmsh_verbosity = gmsh_verbosity
        
        # Optimization problem parameters
        self.variable_array:list = variable_array
        self.dimension = len(variable_array)

        self.units =  '  Mg mm s'
        
        # Default values for configuration parameters
        self.default_parameters = {
            'thickness_tube': 1.8,
            'thickness_layer': 0.7,
            'elsize': 10,
            'nelx': 80,
            'nely_div':2,
            'nelz_div':2,
            'extrusion_length': 800,
            'height': 80,
            'width': 120,
            'node_starting_id': 1001,
            'part_starting_id': 101,
            'id_min': 1000001,
            'elform': 24,
            'nip': 5,
            'shrf': 0.83333,
            'mat_id_tube': 1,
            'mat_id_layer': 2,
            'mat_id_v': 3,
            #'mid_point_node':44721
        }
        self._set_parameters(**kwargs)
        self._generate_database_and_grid()



    def generate_thickness_profiles(self)->List[dict]:
        r"""
        This function generates the thickness profiles to generate the
        varying thicknesses of the ribs

        Returns
        --------------
        A list with (5) dictionaries allocating the thickness profiles
        """


        def return_list_of_thicknesses(points:list)->np.ndarray:
            r"""
            This is a helper function in order to assign the different thicknesses 
            with varying interpolation points
            """

            # Set the special case
            if len(points) == 1:
                return np.asarray(points*8*2**(self.h_level-1))
            
            # Set the general case
            else:
                compare_array = np.linspace(1,8*2**(self.h_level-1),
                                            num=len(points),endpoint=True)
                
                interp_array = np.arange(1,8*2**(self.h_level-1)+1,
                                         1)
                
                interp_thicknesses = np.interp(interp_array,xp=compare_array,
                                               fp=np.asarray(points))
                

                return interp_thicknesses
        
        def generate_dictionaries(thickness_profile:np.ndarray)->dict:

            return {ii+1:float(prof) for ii,prof in enumerate(thickness_profile) }


        # The cases depend on dimensionality
        variable_array = deepcopy(self.variable_array)

        if len(variable_array)==1:
            variable_array = [variable_array[0]]*5
        elif len(variable_array)==2:
            variable_array = [1.7,1.7,1.7, variable_array[0],variable_array[1]]
        elif len(variable_array)==3:
            variable_array = [variable_array[0],1.7,1.7,variable_array[1],variable_array[2]]
        elif len(variable_array)==4:
            variable_array = [1.7,variable_array[0],variable_array[1],variable_array[2],variable_array[3]]


        # Initialize t
        # 
        # he storing bins
        bins = [[] for _ in range(5)]

        multiplier = len(variable_array)//5

        for i, val in enumerate(variable_array):
            # These changes are to match the definition
            # described in the MECHBench paper in ArXiv.

            if i+1 <= 5*multiplier:
                bins[i % 5].append(val)
            else:
                break
        
        # Assign the thicknesses to each bin for the missing values
        if self.dimension > 5:
            unassigned_variables = variable_array[5*multiplier:]
            if len(unassigned_variables)==1:
                bins[0].append(unassigned_variables[0])
            elif len(unassigned_variables)==2:
                bins[3].append(unassigned_variables[0])
                bins[4].append(unassigned_variables[1])
            elif len(unassigned_variables)==3:
                bins[0].append(unassigned_variables[0])
                bins[3].append(unassigned_variables[1])
                bins[4].append(unassigned_variables[2])
            elif len(unassigned_variables)==4:
                bins[1].append(unassigned_variables[0])
                bins[2].append(unassigned_variables[1])
                bins[3].append(unassigned_variables[2])
                bins[4].append(unassigned_variables[3])

        thickness_array = []
        for i,iBin in enumerate(bins):
            thickness_array.append(return_list_of_thicknesses(iBin))

        # Initialize the thickness dictionaries
        thickness_dicts = []

        for ii in range(len(thickness_array)):
            thickness_dicts.append(generate_dictionaries(thickness_array[ii]))
        
        return thickness_dicts


    def volume(self):
        # Calculate the distance between consecutive vertices
        distances = np.linalg.norm(np.diff(self.grid_pts, axis=0, append=[self.grid_pts[0]]), axis=1)
        # Sum the distances to get the perimeter
        perimeter = np.sum(distances)
        result = perimeter*self.extrusion_length*self.thickness
        return result
    
    @property
    def characteristic_length(self)->float:
        r"""Returns the characteristic length of the CrashTubeMesh,
        namely the extrusion length.
        """

        return self.extrusion_length
            
    
    def _set_parameters(self, **kwargs):
        """
        Set parameters for the Three Point Bending based on the provided keyword arguments.

        This method dynamically adds attributes to the instance of the class based on the key-value pairs provided
        as keyword arguments (kwargs). It also sets default values for parameters not provided in kwargs.

        Parameters:
            **kwargs (dict): Keyword arguments to set parameters for the Three Point Bending model.

        Notes:
            - If a parameter is provided in kwargs, it overrides the default value.
            - If 'extrusion_length' is provided in kwargs, it sets the extrusion length of the star box model.
            Otherwise, it sets the extrusion length to 30 times the element size (elsize).
            - If 'trigger_height' is provided in kwargs, it sets the trigger height of the star box model and
            calculates the number of trigger rows based on the provided trigger height and element size (elsize).
            Otherwise, it sets the trigger rows to 3.
            - If the star box model is 3D or 5D, it sets the thickness based on the last value in the variable array.

        """
        # Set parameters dynamically based on kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)
        # Set default parameters
        for key, value in self.default_parameters.items():
            setattr(self, key, kwargs.get(key, value))


    def _generate_database_and_grid(self):
        
        # size of the Tube with Layers
        height= self.height / 2
        width = self.width / 2

        # Initializing the grid points
        grid_list = []
        for aa in np.linspace(-height,height,5,endpoint=True):
            for bb in np.linspace(-width,width,7, endpoint=True):
                grid_list.append((bb,aa))

        self.grid_pts = np.asarray(grid_list)

    
    def write_py_mesh_input_2(self):
        # These would normally be computed or read from elsewhere
        units = "kg mm ms kN GPa kN-mm"
        extrusion_length = self.extrusion_length

        id_min = self.id_min
        #trigger_depth = self.trigger_depth
        #trigger_rows = self.trigger_rows
        #mid = self.mat_id

        grid = [{"gid": int(gid+1), "x": pts[0], "y": pts[1]} for gid, pts in enumerate(zip(self.grid_pts[:,0],  
                                                                       self.grid_pts[:,1]))]


        THICKNESS_MAPS = self.generate_thickness_profiles()


        # Build the full data structure
        data = {
            "units": units,
            "extrusion_length": extrusion_length,
            "h_level":self.h_level,
            "elform": self.elform,
            "nip": self.nip,
            #"shrf": self.shrf,
            "elsize": self.elsize,
            "id_min": id_min,
            #"mid": mid,
            "grid": grid,
            #"cell": cell,
            'thickness_maps':THICKNESS_MAPS,
            #'trigger_dict_1': dict_1,
            #'trigger_dict_2': dict_2,
            "gmsh_verbosity": self.gmsh_verbosity,
            'mat_id_tube': self.mat_id_tube,
            'mat_id_layer': self.mat_id_layer,
            'mat_id_v': self.mat_id_v,
            'nelx': self.nelx,
            'nely_div':self.nely_div,
            'nelz_div':self.nelz_div
        }

        # Output JSON file
        with open("py_mesh_input.json", "w") as f:
            json.dump(data, f, indent=4)



    def write_mesh_file(self):
        # ---- running py_mesh
        #self.write_py_mesh_input()
        #py_mesh_v2('py_mesh.input')

        # Use the GMSH pipeline
        self.write_py_mesh_input_2()
        cl = ThreePointBending_GMSH("three_point_bending_mesh")
        cl("py_mesh_input.json",True)