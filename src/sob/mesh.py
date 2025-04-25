import numpy as np
import os

from .lib.py_mesh import py_mesh
from .lib.py_mesh_v2 import py_mesh_v2

class StarBoxMesh():
    def __init__(self, variable_array, crossing_wall=False, **kwargs) -> None:
        # Optimization problem parameters
        self.variable_array = variable_array
        self.dimension = len(variable_array)

        self.units =  '  kg  mm  ms  kN  GPa  kN-mm'

        self.grid_pts = self._determine_grid_points()
        self._determine_trigger_depth(**kwargs)

        # Default values for configuration parameters
        self.default_parameters = {
            'thickness': [1.2],
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
            'part_starting_id': 101
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
            multiplier = self.extrusion_length/30
            
            # Initialise the result to be 0.00
            result:float = 0.00

            for ii in range(30):
                result += multiplier*perimeter*self.thickness[ii]






        return result

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

        if 'trigger_height' in kwargs:
            self.trigger_height = float(kwargs['trigger_height'])
            self.trigger_rows = int(self.trigger_height / self.elsize)
        else:
            self.trigger_rows = 3

        if self.dimension in [3, 5]:
            self.thickness[0] = self.variable_array[-1]
        
        elif self.dimension in [1,2,4]:
            self.thickness[0] = self.thickness

        elif self.dimension >=5 and self.dimension <= 35:
            # Set some arrays which will be interpolated
            base_array = np.linspace(0.0,
                                     self.extrusion_length,
                                     num=self.dimension-4,
                                     endpoint=True)
            
            range_to_interpolate = self.variable_array[4:]

            # Instantiate a list of thicknesses and positions
            #list_of_thicknesses = np.zeros((30,)).tolist()

            # List of positions
            pos_array = np.linspace(2,self.extrusion_length-2,
                                    num=30,endpoint=True)
            
            # TODO: Modify interpolation // This is a temporary operation
            # Perform the interpolation (set to be linear)

            
            list_of_thicknesses = np.interp(pos_array,
                                            base_array,
                                            range_to_interpolate)
            

            # Assign the thicknesses to the variable
            self.thickness:list = list_of_thicknesses


        
    def _generate_database_and_grid(self, crossing_wall=False):
        # Array stored node ids
        node_hist = np.array([i for i in range(self.node_starting_id, self.node_starting_id + np.size(self.grid_pts, 0))])
        self.database = node_hist
        
        # Stack two arrays horizontally
        self.grid = np.hstack((node_hist.reshape(np.size(self.database), 1), self.grid_pts))
        
        # Generate cells

        ### NOTE: THIS IS A TEST TO CONTINUE DEBUGGING THE CODE
        self.cell = np.array([self.part_starting_id, self.node_starting_id + 1, self.node_starting_id, self.thickness[0]])

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

    def write_mesh_file(self):
        # ---- running py_mesh
        self.write_py_mesh_input()
        py_mesh('py_mesh.input')
        
class CrashTubeMesh():
    def __init__(self, variable_array, **kwargs) -> None:
        # Optimization problem parameters
        self.variable_array = variable_array
        self.dimension = len(variable_array)

        self.units =  '  kg  mm  ms  kN  GPa  kN-mm'

        # size of the crash tube
        length = 120 / 2
        width = 80 / 2
        self.grid_pts = np.array([[-length, width], [length, width],
                                    [length, -width], [-length, -width]])
        self.trigger_positions, self.trigger_depths = self._determine_trigger_depths_and_position()
        # Default values for configuration parameters
        self.default_parameters = {
            'thickness': 1.2,
            'elsize': 4,
            'extrusion_length': 800,
            'node_starting_id': 1001,
            'part_starting_id': 101,
            'id_min': 1000001,
            'trigger_height': 12,  # 3 * elsize = 3 * 4 = 12
            'elform': 2,
            'nip': 3,
            'shrf': 0.83333,
            'mat_id': 999,
            'node_starting_id': 1001,
            'part_starting_id': 101
        }
        self._set_parameters(**kwargs)
        self._generate_database_and_grid()
        
    def volume(self):
        # Calculate the distance between consecutive vertices
        distances = np.linalg.norm(np.diff(self.grid_pts, axis=0, append=[self.grid_pts[0]]), axis=1)
        # Sum the distances to get the perimeter
        perimeter = np.sum(distances)
        result = perimeter*self.extrusion_length*self.thickness
        return result
    
    def _determine_trigger_depths_and_position(self):
        '''
        2 -> three positions and depths vary together with same value. 
        3 -> all three trigger positions fixed, the three trigger depth vary. 
        4 -> except for middle trigger position and depth, other 4 vary ().
        5 -> except for middle trigger depth, other 5 vary. 
        6 -> all three positions and depths vary.
        '''
        
        if self.dimension == 2:
            position = self.variable_array[0]
            depth = self.variable_array[1]
            trigger_positions = [position] * 3
            trigger_depths = [depth] * 3

        elif self.dimension == 3:
            trigger_positions = [0] *3
            trigger_depths = self.variable_array

        elif self.dimension == 4:
            trigger_positions = [self.variable_array[0],0,self.variable_array[2]]
            trigger_depths = [self.variable_array[1],0,self.variable_array[3]]

        elif self.dimension == 5:
            trigger_positions = self.variable_array[:3]
            trigger_depths = self.variable_array[3],0,self.variable_array[4]

        elif self.dimension == 6:
            trigger_positions = self.variable_array[:3]
            trigger_depths = self.variable_array[3:]

        else:
            raise ValueError('Invalid dimension. The dimensions of the problem can only be in [1, 5]')
        
        return np.array(trigger_positions),np.array(trigger_depths)
            
    
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

        if 'trigger_height' in kwargs:
            self.trigger_height = float(kwargs['trigger_height'])
            self.trigger_rows = int(self.trigger_height / self.elsize)
        else:
            self.trigger_rows = 3

        if self.dimension in [3, 5]:
            self.thickness = self.variable_array[-1]
        
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