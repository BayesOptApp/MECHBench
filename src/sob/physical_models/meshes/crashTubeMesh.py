from src.sob.physical_models.meshes.abstractMeshSettings import AbstractMeshSettings
import numpy as np
import os
import json
from src.sob.physical_models.meshes.routines.gmsh.crashtube_gmsh import Crashtube_GMSH
from typing import Union

class CrashTubeMesh(AbstractMeshSettings):
    def __init__(self, variable_array, h_level:int=1,
                 gmsh_verbosity:bool=False, 
                 **kwargs) -> None:
        
        # Set the h_level
        self.h_level = h_level

        self.gmsh_verbosity = gmsh_verbosity
        
        # Optimization problem parameters
        self.variable_array = variable_array
        self.dimension = len(variable_array)

        self.units =  '  kg  mm  ms  kN  GPa  kN-mm'

        # size of the crash tube
        length = 120 / 2
        width = 80 / 2
        self.grid_pts = np.array([[-length, width], [length, width],
                                    [length, -width], [-length, -width]])
        (self.trigger_positions, 
        self.trigger_depths, 
        self.trigger_heights) = self._determine_trigger_depths_heights_and_position()

        self._set_parameters(**kwargs)
        self._generate_database_and_grid()
        
    def volume(self):
        # Calculate the distance between consecutive vertices
        distances = np.linalg.norm(np.diff(self.grid_pts, axis=0, append=[self.grid_pts[0]]), axis=1)
        # Sum the distances to get the perimeter
        perimeter = np.sum(distances)
        result = perimeter*self.extrusion_length*self.thickness
        return result
    
    @property
    def default_parameters(self)->dict:
        r"""
        Return the default parameters for some default configuration
        and mesh setup.
        """
        return {
            'thickness': 1.2,
            'n_elements_side': 12,
            'elsize': 4,
            'extrusion_length': 800,
            'node_starting_id': 1001,
            'part_starting_id': 101,
            'id_min': 1000001,
            'trigger_height': 3,  # 3 * elsize = 3 * 4 = 12
            'elform': 2,
            'nip': 3,
            'shrf': 0.83333,
            'mat_id': 999,
        }
    
    @property
    def characteristic_length(self)->float:
        r"""Returns the characteristic length of the CrashTubeMesh,
        namely the extrusion length.
        """

        return self.extrusion_length
    
    def _determine_trigger_depths_heights_and_position(self)->Union[np.ndarray,
                                                                    np.ndarray,
                                                                    np.ndarray]:
        '''
            Generalized for dimensions from 2 up to 30.

            variable_array is assumed to follow the repeating pattern:
            [depth_1, position_1, height_1, depth_2, position_2, height_2, ...]

            Returns
            ------------------------
            - trigger_positions (np.array): up to 5 vertical positions
            - trigger_depths (np.array): up to 5 depths
            - trigger_heights (np.array): up to 5 heights
        '''

        if self.dimension < 1 or self.dimension > 30:
            raise ValueError('Invalid dimension. The dimension must be in [1, 30]')

        if self.dimension <=15:
            trigger_positions = [0.0] * 5
            trigger_depths = [0.0] * 5
            trigger_heights = [self.default_parameters['trigger_height']] * 5
        else:
            trigger_positions = [0.0] * 10
            trigger_depths = [0.0] * 10
            trigger_heights = [self.default_parameters['trigger_height']] * 10

        num_vars = min(self.dimension, len(self.variable_array))

        for i in range(num_vars):
            group_index = i//3   # 0, 1, 2 for each trigger group
            role_index = np.remainder(i, 3)                # 0: depth, 1: position, 2: height

            if role_index == 0:
                trigger_depths[group_index] = self.variable_array[i]
            elif role_index == 1:
                trigger_positions[group_index] = self.variable_array[i]
            elif role_index == 2:
                trigger_heights[group_index] = self.variable_array[i]

            

        if self.dimension <=15:
            return (
                np.tile(np.array(trigger_positions),(2,)),
                np.tile(np.array(trigger_depths),(2,)),
                np.tile(np.array(trigger_heights),(2,))
            )
        else:
            return (
                np.array(trigger_positions),
                np.array(trigger_depths),
                np.array(trigger_heights)
            )
            
    
    def _set_parameters(self, **kwargs):
        """
        Set parameters for the star box model based on the provided keyword arguments.

        This method dynamically adds attributes to the instance of the class based on the key-value pairs provided
        as keyword arguments (kwargs). It also sets default values for parameters not provided in kwargs.

        Parameters:
            **kwargs (dict): Keyword arguments to set parameters for the star box model.

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

        # if 'extrusion_length' in kwargs:
        #     self.extrusion_length = float(kwargs['extrusion_length'])
        #     self.extrusion_length = int(self.extrusion_length / self.elsize) * self.elsize
        # else:
        #     self.extrusion_length = self.elsize * 30

        # if 'trigger_height' in kwargs:
        #     self.trigger_height = float(kwargs['trigger_height'])
        #     self.trigger_rows = int(self.trigger_height / self.elsize)
        # else:
        #     self.trigger_rows = 3
        
    def _generate_database_and_grid(self):
        # Array stored node ids
        node_hist = np.array([i for i in range(self.node_starting_id, self.node_starting_id + np.size(self.grid_pts, 0))])
        self.database = node_hist
        
        # Stack two arrays horizontally
        self.grid = np.hstack((node_hist.reshape(np.size(self.database), 1), self.grid_pts))
        
        # Generate cells
        self.cell = np.array([self.part_starting_id, self.node_starting_id + 1, self.node_starting_id, self.thickness])
        for i in range(1, np.size(self.grid_pts, 0) - 1):
            if np.mod(i, 2) == 1:
                cell_row = np.array([self.part_starting_id + i, self.node_starting_id + i, self.node_starting_id + i + 1, self.thickness])
            else:
                cell_row = np.array([self.part_starting_id + i, self.node_starting_id + i + 1, self.node_starting_id + i, self.thickness])
            self.cell = np.vstack((self.cell, cell_row))
        
        # Last row of cells
        cell_row = np.array([self.part_starting_id + i + 1, self.node_starting_id + i + 1, self.node_starting_id, self.thickness])
        self.cell = np.vstack((self.cell, cell_row))
        

    def write_py_mesh_input(self): 
        """ 
        writes py_mesh.input file 
        
        Inputs: 
        element_size
        element_form
        intergration_points
        element_shear_factor
        extrusion_length
        starting_id
        trigger_depth
        trigger_rows
        material_id
        database = node_hist

        Outputs
        
        """  
        adr = os.path.join(os.getcwd(),'py_mesh.input')
        inf = open(adr,'w')
        inf.writelines('#  units:' + self.units + '\n')
        inf.writelines('\n# ----- height of the structure (i.e. extrusion length)\n')
        inf.writelines('extrusion_length, %f\n' %self.extrusion_length)
        # ----- element related
        inf.writelines('\n# ----- element related\n')
        inf.writelines('elform, %d\n' %self.elform)
        inf.writelines('nip, %d\n' %self.nip)  
        inf.writelines('shrf, %f\n' %self.shrf)   
        inf.writelines('elsize, %f\n' %self.elsize) 
        inf.writelines('\n# ----- The id number for the first node and shell\n')
        inf.writelines('id_min, %d\n' %self.id_min)
        # ----- trigger 
        inf.writelines('\n# ----- Trigger\n')
        inf.writelines('trigger_rows, %d\n' %self.trigger_rows)
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #  added for trigger variables (Li)
        trigger_positions_str = 'trigger_positions, ' + ', '.join(map(str, self.trigger_positions))
        inf.writelines(f'{trigger_positions_str}\n')
        trigger_depths_str = 'trigger_depths, ' + ', '.join(map(str, self.trigger_depths))
        inf.writelines(f'{trigger_depths_str}\n')
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        inf.writelines('\n')  
        for i in range(len(self.cell[:,0])):
            inf.writelines('trigger, %d\n' %self.cell[i,0])
        # ----- material id
        inf.writelines('\n# ----- Material id\n')
        inf.writelines('mid, %d\n' %self.mat_id)
        # ----- writing database history node
        inf.writelines('\n# ----- Define a *DATABASE_HISTORY_NODE keyword\n')
        for i in range(len(self.database)):
            inf.writelines('database, %d\n' %self.database[i])     
        # ----- writing a grid
        inf.writelines('\n# ----- Define nodes that define the geometry\n')  
        max_grid_id = len(str(int(np.max(self.grid[:,0]))))  #grid_?
        str_grd_title = "{:>5}{:>"+str(max_grid_id+2)+"}{:>17}{:>16}\n"
        inf.writelines(str_grd_title.format('#   ,','gid,', 'x,', 'y'))
        str_grd = "{:>4}{:>"+str(max_grid_id+2)+"}{:>17}{:>16}\n"
        for i in range(np.size(self.grid,0)):
            inf.writelines(str_grd .format('grid,',str(int(self.grid[i,0]))+',',
                            str(int(self.grid[i,1]))+',', str(int(self.grid[i,2])))) #grid_?
        # ----- writing cells
        inf.writelines('\n# ----- Define parts') 
        max_cell_id = len(str(int(np.max(self.cell[:,0]))))
        str_cell_title = "{:>5}{:>"+str(max_cell_id+3)+"}{:>17}{:>16}{:>16}\n"  
        inf.writelines(str_cell_title.format('\n#   ,','cid,','g0,','g1,','t')) 
        str_cell = "{:>4}{:>"+str(max_cell_id+3)+"}{:>17}{:>16}{:>16}\n"
        for i in range(np.size(self.cell,0)):
            inf.writelines(str_cell.format('cell,',str(int(self.cell[i,0]))+',',
                            str(int(self.cell[i,1]))+',', str(int(self.cell[i,2]))+',',str(self.cell[i,3])))          
        inf.close()
    
    def build_trigger_directories(self):
        r"""
        This function builds the 'trigger directory' in order to set the properties
        needed to build-up the triggers in the mesh 
        """

        # Transform the list of trigger positions and trigger heights into by mapping
        # into "physical space" and not defined by the mesh.

        n_elements_z = int(self.extrusion_length / self.elsize)


        # Make a list of the positions in "mesh definition" of each trigger
        
        end_points_segments = np.linspace(0, n_elements_z, num=6, endpoint=True)

        # Get the neutral points of each segment
        neutral_points_segments = end_points_segments[:-1] + (end_points_segments[1:] - end_points_segments[:-1]) / 2
        neutral_points_segments = np.tile(np.flip(neutral_points_segments), (2, 1)).flatten()

        # Get the trigger positions in the mesh definition
        trigger_positions_complete = np.add(np.round(self.trigger_positions,
                                                     decimals=0), neutral_points_segments)

        # Get splitted 
        # NOTE: The flips is to ensure that the triggers are defined from top to bottom
        # 02/12/2025 - Adopted to match with the definition shown in the first
        #              version of the MECHBench paper in ArXiv
        trigger_positions_1 = np.flip((trigger_positions_complete[:5]*self.elsize)).astype(int)
        trigger_positions_2 = np.flip((trigger_positions_complete[5:]*self.elsize)).astype(int)

        trigger_heights_1 = np.flip((np.round(self.trigger_heights[:5],0)*self.elsize).astype(int))
        trigger_heights_2 = np.flip((np.round(self.trigger_heights[5:], 0)*self.elsize).astype(int))

        trigger_depths_1 = np.flip(self.trigger_depths[:5]).astype(float)
        trigger_depths_2 = np.flip(self.trigger_depths[5:]).astype(float)


        dict_1 = {
                i + 1: {
                    'Y': int(trigger_positions_1[i]),
                    'HEIGHT': int(trigger_heights_1[i]),
                    'DEPTH': float(trigger_depths_1[i])
                }
                for i in range(5)
            }
        
        dict_2 = {
                i + 1: {
                    'Y': int(trigger_positions_2[i]),
                    'HEIGHT': int(trigger_heights_2[i]),
                    'DEPTH': float(trigger_depths_2[i])
                }
                for i in range(5)
            }
        
        return dict_1, dict_2


    def write_py_mesh_input_2(self):
        # These would normally be computed or read from elsewhere
        units = "kg mm ms kN GPa kN-mm"
        extrusion_length = self.extrusion_length

        id_min = self.id_min
        #trigger_depth = self.trigger_depth
        #trigger_rows = self.trigger_rows
        mid = self.mat_id

        grid = [{"gid": int(gid), "x": x, "y": y} for gid, x, y in zip(self.grid[:,0], 
                                                                       self.grid[:,1], 
                                                                       self.grid[:,2])]

        thick_array = [self.thickness]

        cell = [{"cid": cid, "t": t} for cid, t in zip(range(self.part_starting_id, 
                                                             self.part_starting_id + len(thick_array)), 
                                                             thick_array)]
        
        (dict_1, 
         dict_2) = self.build_trigger_directories()

        # Build the full data structure
        data = {
            "units": units,
            "extrusion_length": extrusion_length,
            "h_level":self.h_level,
            "elform": self.elform,
            "nip": self.nip,
            "shrf": self.shrf,
            "elsize": self.elsize,
            "id_min": id_min,
            "mid": mid,
            "grid": grid,
            "cell": cell,
            'n_elements_side':self.n_elements_side,
            'thickness':self.thickness,
            'trigger_dict_1': dict_1,
            'trigger_dict_2': dict_2,
            "gmsh_verbosity": self.gmsh_verbosity,
        }

        # Output JSON file
        with open("py_mesh_input.json", "w") as f:
            json.dump(data, f, indent=4)


    def write_mesh_file(self):
        # ---- running py_mesh
        # self.write_py_mesh_input()
        # py_mesh_v2('py_mesh.input')

        # Use the GMSH pipeline
        self.write_py_mesh_input_2()
        cl = Crashtube_GMSH("crash_tube_mesh")
        cl("py_mesh_input.json",True)