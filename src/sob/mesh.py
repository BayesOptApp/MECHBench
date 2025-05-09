import numpy as np
import os
from abc import ABC, abstractmethod
from .lib.py_mesh import py_mesh
from .lib.py_mesh_v2 import py_mesh_v2
from .lib.starbox_gmsh import Starbox_GMSH
from .lib.crashtube_gmsh import Crashtube_GMSH
from typing import Union, List, Optional
from .lib.crashtube_gmsh import Crashtube_GMSH
from typing import Union, List, Optional
import json



class AbstractMeshSettings(ABC):
    r"""
    This is an abstract definition for the mesher settings
    """

    @abstractmethod
    def __init__(self, variable_array, h_level,**kwargs)->None:
        """
        Actual class initializer
        """
        pass

    @abstractmethod
    def volume(self)->float:
        r"""Returns the volume (in object units) of the structural object.       
        """
        pass

    @abstractmethod
    def write_mesh_file(self)->None:
        r"""
        Writes the mesh file for the simulator (LS-Dyna/OpenRadioss)
        """
        pass

    @property
    @abstractmethod
    def characteristic_length(self)->None:
        r"""
        Return the characteristic length of the structure. This is a property to set up
        computations of other Crashworthiness Objective metrics. Namely, for future additions,
        the characteristic length is the length of the structure in the direction of the impact.
        """
        pass

    @property
    def h_level(self)->int:
        r"""
        Mesh refinement level
        """
        return self.__h_level
    
    @h_level.setter
    def h_level(self,new_h_level:int)->None:
        r"""
        Mesh refinement setter

        Args
        ----------------------
        new_h_level (`int`): An integer defining the mesh refinement level
        """

        if isinstance(new_h_level,(int,float)) and int(new_h_level) >= 1:
            # Set the h_level in this case
            self.__h_level = int(new_h_level)
        else:
            raise ValueError("The mesh refinement level `h_level` must be a positive integer greater than 1.")

class StarBoxMesh(AbstractMeshSettings):
    def __init__(self, variable_array, h_level:int=1, crossing_wall=False, **kwargs) -> None:
        r"""
        Star-Box Mesh Initializer

        Args
        --------------------------
        - variable_array (`list`): An array with the variable physical mapping
        - h_level (`int`): Mesh Refinement level
        - crossing_wall (`bool`): A switch to allow the wall to cross the material

        """
        self.h_level = h_level
        # Optimization problem parameters
        self.variable_array = variable_array
        self.dimension = len(variable_array)

        self.units =  '  kg  mm  ms  kN  GPa  kN-mm'

        self.grid_pts = self._determine_grid_points()
        self._determine_trigger_depth(**kwargs)

        # Default values for configuration parameters
        self.default_parameters = {
            'thickness': 1.2,
            'elsize': 4,
            'extrusion_length': 120,
            'node_starting_id': 1001,
            'part_starting_id': 101,
            'id_min': 1000001,
            'trigger_height': 12,  # 3 * elsize = 3 * 4 = 12
            'elform': 2,
            'nip': 3,
            'shrf': 0.83333,
            'mat_id': 999,
            'node_starting_id': 1001,
            'n_elements_side':12
        }
        self._set_parameters(**kwargs)
        self._generate_database_and_grid(crossing_wall)

    def volume(self):
        # Calculate the distance between consecutive vertices
        distances = np.linalg.norm(np.diff(self.grid_pts, axis=0, append=[self.grid_pts[0]]), axis=1)
        # Sum the distances to get the perimeter
        perimeter = np.sum(distances)

        if self.dimension <=5:
            result = perimeter*self.extrusion_length*self.thickness
        elif self.dimension > 5 and self.dimension < 36:
            # Divide the extrusion length into 30 parts
            multiplier = self.extrusion_length/self.num_layers
            
            # Initialise the result to be 0.00
            result:float = 0.00

            for ii in range(30):
                result += multiplier*perimeter*self.thickness[ii]

        return result
    
    @property
    def characteristic_length(self)->float:
        r"""Returns the characteristic length of the StarBoxMesh,
        namely the extrusion length.
        """

        return self.extrusion_length

    def _determine_grid_points(self):
        """
        Determine the coordinates of the grid points to define the shape of the star box.

        This function calculates the coordinates of the vertices based on the input variable array
        and stores them into the attribute `self.grid_pts`, which is used as input for mesh generation.

        Raises:
            ValueError: If the dimension is invalid or not supported.

        Notes:
            - For 1D problems, the grid points are calculated based on the length of the star box.
            - For 2D and 3D problems, the grid points are calculated based on the length and width of the star box.
            - For 4D and 5D problems, a specialized method `_create_grid_pts_Hunkler` is called to determine
            the grid points.
        """
        if self.dimension == 1:
            length = self.variable_array[0] / 2
            return np.array([[-length, length], [length, length],
                                      [length, -length], [-length, -length]])
        elif self.dimension in [2, 3]:
            length = self.variable_array[0] / 2
            width = self.variable_array[1] / 2
            return np.array([[-length, width], [length, width],
                                      [length, -width], [-length, -width]])
        elif self.dimension >= 4:
            return self._create_grid_pts_Hunkler()
        else:
            raise NotImplementedError('Invalid dimension. The dimensions of the problem can only be in [1, 35]')
        
    def _create_grid_pts_Hunkler(self):
        """
        this function creates vertices for the hunkler star model. these vertices will
        be connected by the py_mesh to create the cross section of the crash box, which then 
        will be extruded later
        """
        a = self.variable_array[0]
        b = self.variable_array[1]
        u = self.variable_array[2]
        v = self.variable_array[3]
        return  np.array([[b/2.,a/2.],[0,a/2.-u],[-b/2.,a/2.],[-b/2.+v,0],
                        [-b/2.,-a/2.],[0,-a/2.+u],[b/2.,-a/2.],[b/2.-v,0]])
    
    def _determine_trigger_depth(self, **kwargs):
        """
        Determine the default trigger depth based on the size of the star box.

        This method calculates the default trigger depth based on the dimensions of the star box
        and any optional user-defined values provided via keyword arguments.

        Parameters:
            **kwargs (dict): Optional keyword arguments to override the default trigger depth.

        Notes:
            - The default trigger depth is calculated based on the dimensions of the star box.
            - For 1D problems, the trigger depth is a percentage (default 5%) of the length of the star box.
            - For 2D and 3D problems, the trigger depth is a percentage (default 5%) of the average of the
            length and width of the star box.
            - For 4D and 5D problems, the trigger depth is a percentage (default 5%) of the average of the
            diagonal lengths of the star box's cross-sections.

        Returns:
            float: The calculated default trigger depth.

        """
        if self.dimension == 1:
            length = self.variable_array[0] / 2
            self.trigger_depth = float(kwargs.get('trigger_depth', 0.05 * length))
        elif self.dimension in [2, 3]:
            length = self.variable_array[0] / 2
            width = self.variable_array[1] / 2
            self.trigger_depth = float(kwargs.get('trigger_depth', 0.05 * (length + width) / 2))
        elif self.dimension >= 4:
            a = self.variable_array[0]
            b = self.variable_array[1]
            u = self.variable_array[2]
            v = self.variable_array[3]
            self.trigger_depth = float(kwargs.get('trigger_depth', 0.05 * (np.sqrt((a / 2.) ** 2 + v ** 2) + np.sqrt(((b / 2.) ** 2 + u ** 2))) / 2)) 

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
        
        

        if 'extrusion_length' in kwargs:
            self.extrusion_length = float(kwargs['extrusion_length'])
            self.extrusion_length = int(self.extrusion_length / self.elsize) * self.elsize
        else:
            self.extrusion_length = self.elsize * 30
        
        # Set the number of layers of the mesh
        self.num_layers = 30*2**(self.h_level-1)

        if 'trigger_height' in kwargs:
            self.trigger_height = float(kwargs['trigger_height'])
            self.trigger_rows = int(self.trigger_height / self.elsize)
        else:
            self.trigger_rows = 3

        temp_thickness = self.default_parameters['thickness']
        self.thickness = []

        if self.dimension in [3, 5]:
            self.thickness.append(self.variable_array[-1])
        
        elif self.dimension in [1,2,4]:
            self.thickness.append(temp_thickness)

        elif self.dimension >5 and self.dimension <= 35:
            # Set some arrays which will be interpolated

            init_pos = self.elsize/2/2**(self.h_level-1)


            base_array = np.linspace(init_pos,self.extrusion_length-init_pos,
                                     num=self.dimension-4,
                                     endpoint=True)
            
            range_to_interpolate = self.variable_array[4:]


            
            # List of positions
            pos_array = np.linspace(init_pos,self.extrusion_length-init_pos,
                                    num=self.num_layers,endpoint=True)
            
            # TODO: Modify interpolation // This is a temporary operation
            # Perform the interpolation (set to be linear)

            
            list_of_thicknesses = np.interp(pos_array,
                                            base_array,
                                            range_to_interpolate).tolist()
            

            # Assign the thicknesses to the variable
            self.thickness:list = list_of_thicknesses


        
    def _generate_database_and_grid(self, crossing_wall=False):
        # Array stored node ids
        node_hist = np.array([i for i in range(self.node_starting_id, 
                                               self.node_starting_id + 
                                               np.size(self.grid_pts, 0))])
        self.database = node_hist
        
        # Stack two arrays horizontally
        self.grid = np.hstack((node_hist.reshape(np.size(node_hist), 1), self.grid_pts))
        
        # Generate cells

        ### NOTE: THIS IS A TEST TO CONTINUE DEBUGGING THE CODE
        
        self.cell = np.array([self.part_starting_id, 
                                self.node_starting_id + 1, 
                                self.node_starting_id, 
                                np.asarray(self.thickness[0]).ravel()[0]])
        

        ### PREVIOUS LINE TO BE REMOVED
        for i in range(1, np.size(self.grid_pts, 0) - 1):
            if np.mod(i, 2) == 1:
                cell_row = np.array([self.part_starting_id + i, self.node_starting_id + i, self.node_starting_id + i + 1, self.thickness[0]])
            else:
                cell_row = np.array([self.part_starting_id + i, self.node_starting_id + i + 1, self.node_starting_id + i, self.thickness[0]])
            self.cell = np.vstack((self.cell, cell_row))
        
        # Last row of cells
        cell_row = np.array([self.part_starting_id + i + 1, self.node_starting_id + i + 1, self.node_starting_id, self.thickness[0]])
        self.cell = np.vstack((self.cell, cell_row))

        # -- star box with crossing wall
        if crossing_wall:
            # 1. Redefine the grid data frame, add additional column of the center point coordinate  
            node_hist = np.array([j for j in range(self.node_starting_id,self.node_starting_id+np.size(self.grid_pts,0)+1)])
            # node_hist2 = node_hist.reshape(np.size(node_hist),1)
            self.database = node_hist
            new_grid_pts = np.vstack([self.grid_pts, [0.0, 0.0]])
            self.grid=np.hstack((node_hist.reshape(np.size(self.database),1), new_grid_pts))
            
            # 2. Add additional parts (line), stacking new rows to the cell data frame
            cell_id = self.part_starting_id+7
            for j in range(1,np.size(self.grid_pts,0)+1):
                if np.mod(j,2) == 1:
                    cell_id += 1
                    self.cell = np.vstack((self.cell,np.array([cell_id, node_hist[-1], self.node_starting_id+j, self.thickness[0]])))

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
        inf.writelines('trigger_depth, %f\n' %self.trigger_depth)
        inf.writelines('trigger_rows, %d\n' %self.trigger_rows)
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

    def write_py_mesh_input_2(self):
        # These would normally be computed or read from elsewhere
        units = "kg mm ms kN GPa kN-mm"
        extrusion_length = self.extrusion_length

        id_min = self.id_min
        trigger_depth = self.trigger_depth
        trigger_rows = self.trigger_rows
        mid = self.mat_id

        grid = [{"gid": int(gid), "x": x, "y": y} for gid, x, y in zip(self.grid[:,0], 
                                                                       self.grid[:,1], 
                                                                       self.grid[:,2])]

        if len(self.thickness)==1:
            thick_array = [self.thickness[0]]*self.num_layers
        else:
            thick_array = self.thickness.copy()

        cell = [{"cid": cid, "t": t} for cid, t in zip(range(self.part_starting_id, 
                                                             self.part_starting_id + len(thick_array)), 
                                                             thick_array)]

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
            "trigger_depth": trigger_depth,
            "trigger_rows": trigger_rows,
            "mid": mid,
            "grid": grid,
            "cell": cell,
            "num_layers":self.num_layers,
            'n_elements_side':self.n_elements_side
        }

        # Output JSON file
        with open("py_mesh_input.json", "w") as f:
            json.dump(data, f, indent=4)

    def write_mesh_file(self):
        # ---- running py_mesh
        # if self.dimension <=5:
        #     self.write_py_mesh_input()
        #     py_mesh('py_mesh.input')
        # else:
        self.write_py_mesh_input_2()
        cl = Starbox_GMSH("star_box_mesh")
        cl("py_mesh_input.json",True)
        
        
class CrashTubeMesh(AbstractMeshSettings):
    def __init__(self, variable_array, h_level:int=1, **kwargs) -> None:
        
        # Set the h_level
        self.h_level = h_level
        
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
            'node_starting_id': 1001,
            'part_starting_id': 101
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
            [position_1, depth_1, height_1, position_2, depth_2, height_2, ...]

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
            role_index = np.remainder(i, 3)                # 0: position, 1: depth, 2: height

            if role_index == 0:
                trigger_positions[group_index] = self.variable_array[i]
            elif role_index == 1:
                trigger_depths[group_index] = self.variable_array[i]
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

        n_elements_z = self.extrusion_length / self.elsize


        # Make a list of the positions in "mesh definition" of each trigger
        
        end_points_segments = np.linspace(0, n_elements_z, num=6, endpoint=True)

        neutral_points_segments = end_points_segments[:-1] + (end_points_segments[1:] - end_points_segments[:-1]) / 2
        neutral_points_segments = np.tile(neutral_points_segments, (2, 1)).flatten()

        # Get the trigger positions in the mesh definition
        trigger_positions_complete = np.add(np.round(self.trigger_positions,
                                                     decimals=0), neutral_points_segments)

        # Get splitted 
        trigger_positions_1 = (trigger_positions_complete[:5]*self.elsize).astype(int)
        trigger_positions_2 = (trigger_positions_complete[5:]*self.elsize).astype(int)

        trigger_heights_1 = (np.ceil(self.trigger_heights[:5])*self.elsize).astype(int)
        trigger_heights_2 = (np.ceil(self.trigger_heights[5:])*self.elsize).astype(int)

        trigger_depths_1 = (self.trigger_depths[:5])
        trigger_depths_2 = (self.trigger_depths[5:])



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
            'trigger_dict_2': dict_2
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




class ThreePointBendingMesh(AbstractMeshSettings):
    def __init__(self, variable_array, h_level:int=1, **kwargs) -> None:
         # Set the h_level
        self.h_level = h_level
        
        # Optimization problem parameters
        self.variable_array = variable_array
        self.dimension = len(variable_array)

        self.units =  '  kg  mm  ms  kN  GPa  kN-mm'
        
        # Default values for configuration parameters
        self.default_parameters = {
            'thickness': 1.2,
            'elsize': 4,
            'extrusion_length': 800,
            'node_starting_id': 1001,
            'part_starting_id': 101,
            'id_min': 1000001,
            'elform': 2,
            'nip': 5,
            'shrf': 0.83333,
            'mat_id': 999,
            'node_starting_id': 1001,
            'part_starting_id': 101,
            'mid_point_node':44721
        }
        self._set_parameters(**kwargs)
        self._generate_database_and_grid()

        # size of the Tube with Layers
        length = 120 / 2
        width = 80 / 2
        self.grid_pts = np.array([[-length, width], [length, width],
                                    [length, -width], [-length, -width]])
        
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



    def write_mesh_file(self):
        # ---- running py_mesh
        self.write_py_mesh_input()
        py_mesh_v2('py_mesh.input')