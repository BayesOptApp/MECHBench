import os
import pandas as pd
import shutil
from .solver import run_radioss
from .mesh import *
from .fem import *

'''
Remark: following phases are important 

1.  
'''

class OptiProblem():
    '''
    The abstract of the structural optimization problem, with three subclass (three problem types)
    The instance is initialized by problem dimension, output data type, and the batch file path of the solver OpenRadioss. 
    The core idea is to generate a specific type of instance, input the variable array then output the evaluation of the desired data type like mass, intrusion, etc.
    The input variables should be in an universal search space, default (-5,5). For the evaluation those variables will be mapped to the real FEM problem space.
    '''
    def __init__(self, dimension:int, output_data:(str|list), batch_file_path:str) -> None:
        self.dimension = dimension
        self.output_data:(str|list) = output_data
        self.batch_file_path = batch_file_path

        # The attributes need to be overwritten in teh subclass
        self.problem_id = None
        self.variable_ranges = None # constraints of the problem
        self.input_file_name = None # input deck name
        self.output_file_name = None # output result name
        self.starter_out_file_name = None

        # The attributes will be loaded if the function _write_input_file has been called. 
        # Used for mass calculation.
        self.mesh = None
        self.model = None

        # The attributes will be loaded if the function run_simulation has been called.
        self.output_data_frame = None
        
        # Universal design variable search space
        self.search_space = (-5.0, 5.0)

        # TODO: This is a dummy variable for defining the simulation status:
        #   0 -> Simulation not started
        #   1 -> Just starter launched
        #   2 -> Full simulation

        self.sim_status:int = 0

    def _validate_variable_array(self, variable_array):
        """
        Validate the variable array against the search space.

        Parameters:
            variable_array (list): The array of variables to be validated.

        Raises:
            ValueError: If the size of the variable array does not match the problem dimension or if any variable is out of range.
        """
        if self.variable_ranges is None:
            raise ValueError("variable_ranges must be provided or defined in the subclass.")
        
        if len(variable_array) != self.dimension:
            raise ValueError('The size of variable array does not match the problem dimension')

        for i, (value, (lower, upper)) in enumerate(zip(variable_array, [self.search_space]*self.dimension)):
            if not (lower <= value <= upper):
                raise ValueError(f"Value at position {i} in variable_array is out of range: {value}. " f"Allowed range is [{lower}, {upper}].")
    
    def linear_maping_variable(self, search_space_variable, problem_space_range:tuple):
        """
        Map a variable from the search space to the problem space for use in FEM simulation.

        Parameters:
            search_space_variable (float): Variable in the search space (for optimization).
            problem_space_range (tuple): Range of variables in the problem space (for FEM simulation).

        Returns:
            float: Variable mapped to the problem space.
        """
        lower = problem_space_range[0]
        upper = problem_space_range[1]
        scale = (upper-lower)/(self.search_space[1]-self.search_space[0])
        problem_space_variable = lower + (search_space_variable-self.search_space[0])*scale
        return problem_space_variable

    def generate_input_deck(self, variable_array):
        # Change the simulation status
        self.sim_status = 0
        self._validate_variable_array(variable_array)

        fem_space_variable_array = [] # to get the variables in the FEM space
        for i, var in enumerate(variable_array):
            mapped_var = self.linear_maping_variable(var, self.variable_ranges[i])
            fem_space_variable_array.append(mapped_var)

        original_dir = os.getcwd()
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        print('######################################################\n')
        print(dir_name)
        working_dir = os.path.join(os.getcwd(), dir_name)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        os.chdir(working_dir)
        self._write_input_file(fem_space_variable_array)
        os.chdir(original_dir)

        self.problem_id += 1 # update the problem id for the input deck generation
        
    def _write_input_file(self, fem_space_variable_array):
        raise NotImplementedError("Subclasses must implement _write_input_file method")
    
    def run_simulation(self,runStarter=False):
        if self.input_file_name is None:
            raise ValueError("input_file_name must be provided or defined in the subclass.")
        # make problem id back to original, since it has been updated when generate_input_deck has been called
        self.problem_id -= 1
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        input_file_path = os.path.join(working_dir, self.input_file_name)
        run_radioss(input_file_path, self.batch_file_path, runStarter=runStarter)

    def load_output_data_frame(self):
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        # load simulation result dataframe and make it cleaner
        output_file_path = os.path.join(working_dir, self.output_file_name)
        self.output_data_frame = pd.read_csv(output_file_path)
        self.output_data_frame.columns = self.output_data_frame.columns.str.replace(' ', '')
        # update problem id again
        self.problem_id += 1

    def extract_mass_from_file(self):
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        # load simulation result dataframe and make it cleaner
        starter_out_file_path = os.path.join(working_dir, self.starter_out_file_name)
        # Open the file and read its contents
        with open(starter_out_file_path, 'r') as file:
            lines = file.readlines()
        # Define the marker for where the mass value starts
        mass_marker = "TOTAL MASS AND MASS CENTER"
        
        # Iterate through each line to find the marker
        for i, line in enumerate(lines):
            if mass_marker in line:
                # The mass value is two lines after the marker
                mass_line = lines[i+4].strip()
                # Split the mass line and extract the first value, which is the mass
                mass_value = mass_line.split()[0]
                return float(mass_value)

        return ValueError("Mass value output failed!")  # Return None if the mass value wasn't found 

    def __call__(self, variable_array):

        ## --------------------------------------------
        ## TODO: This functions are just to be recycled
        ## --------------------------------------------

        def mass_calculation(obj_:OptiProblem)->float:
            if obj_.sim_status < 1: 
                obj_.run_simulation(runStarter=True)
                # Change the status
                obj_.sim_status = 1

            val:float = obj_.extract_mass_from_file() - obj_.model.rigid_mass
            obj_.problem_id +=1    
            return val

        def absorbed_energy_calculation(obj_:OptiProblem)->float:
            return obj_.model.absorbed_energy()
        
        def instrusion_calculation(obj_:OptiProblem)->float:
            obj_.run_simulation()
            obj_.load_output_data_frame()
            matching_columns = [col for col in obj_.output_data_frame.columns if obj_.track_node_key in col]
            max_z:float = max(obj_.output_data_frame[matching_columns[2]], key=abs)

            # Change the sim status
            self.sim_status = 2
            return abs(max_z)

        # This line is everpresent...
        # First generate all the input deck for OpenRadioss
        self.generate_input_deck(variable_array)

        # This structure just checks if the kind of output data is a string
        # or lone variable
        if isinstance(self.output_data,str):
            if self.output_data == 'mass':
                result = mass_calculation(self,True)
            elif self.output_data == 'absorbed_energy':
                result = absorbed_energy_calculation(self)
            elif self.output_data == 'intrusion':
                result = instrusion_calculation(self)
            elif self.output_data == 'specific_energy':
                # Compute first the absorbed energy
                abs_ene:float = absorbed_energy_calculation(self)

                # Compute the mass
                mass:float = mass_calculation(self)

                result = abs_ene/mass       
            else:
                result = np.nan
        elif isinstance(self.output_data,list):

            resp_arr = list()
            looping_array = self.output_data.copy()
            idx:int = -1
            # Check is intrusion is in the list
            if 'intrusion' in self.output_data:
                # Get the position of "instrusion" in the array
                idx:int = self.output_data.index('intrusion')

                looping_array.pop(idx)
            
            # Now loop the calculations 
            for out_ in looping_array:
                if out_ == 'mass':
                    resp_arr.append( mass_calculation(self) )
                elif out_ == 'absorbed_energy':
                    resp_arr.append( absorbed_energy_calculation(self) )
                elif out_ == 'specific_energy':
                    # Compute first the absorbed energy
                    abs_ene:float = absorbed_energy_calculation(self)

                    # Compute the mass
                    mass:float = mass_calculation(self)

                    resp_arr.append( abs_ene/mass )      
                else:
                    resp_arr.append( np.nan )

            if idx !=-1:
                # Append
                # perform the intrusion calculation first
                instrusion_calc = instrusion_calculation(self)
                resp_arr.insert( idx,instrusion_calc)

            return resp_arr
            
            
        else:
            result = None
        
        return result        


class StarBox(OptiProblem):
    instance_counter = 1
    def __init__(self, dimension, output_data, batch_file_path) -> None:
        super().__init__(dimension, output_data, batch_file_path)
        # 1 -> square
        # 2 -> rectangular
        # 3 -> rectangular with varying thickness
        # 4 -> star shape
        # 5 -> star shape with varying thickness
        variable_ranges_map = {
            1: [(60, 120)],
            2: [(60, 120), (60, 120)],
            3: [(60, 120), (60, 120), (0.7, 3)],
            4: [(60, 120), (60, 120), (0, 30), (0, 30)],
            5: [(60, 120), (60, 120), (0, 30), (0, 30), (0.7, 3)]
        }
        self.variable_ranges = variable_ranges_map[self.dimension]
        self.input_file_name = 'combine.k'
        self.output_file_name = 'combineT01.csv'
        self.starter_out_file_name = 'combine_0000.out'
        self.track_node_key = 'DATABASE_HISTORY_NODE1001' # the key of the intrusion in the output csv file

        self.problem_id = StarBox.instance_counter
        StarBox.instance_counter+=1

    def _write_input_file(self, fem_space_variable_array):
        self.mesh = StarBoxMesh(fem_space_variable_array) 
        self.model = StarBoxModel(self.mesh)
        self.model.write_input_files()
    

class ThreePointBending(OptiProblem):
    '''
    This optimization problem focuses on adjusting thickness, unlike the Star Box and Crash Tube problems that involve remeshing. 
    The key process involves extracting the section of the input file related to thickness, modifying it, and merging it back with the original input to generate a new starter file.
    '''
    instance_counter = 1
    def __init__(self, dimension, output_data, batch_file_path) -> None:
        super().__init__(dimension, output_data, batch_file_path)
        # 1 -> all 5 shell thickness vary with same value. 
        # 2 -> only first and the last shell thickness vary, other fixed with middle value. 
        # 3 -> the first, middle, and last shell thickness vary. 
        # 4 -> expect for middle shell, the other 4 shell thickness vary.
        # 5 -> all three shell thickness vary.
        variable_ranges_map = {
            1: [(0.5, 3)],
            2: [(0.5, 3)]*2,
            3: [(0.5, 3)]*3,
            4: [(0.5, 3)]*4,
            5: [(0.5, 3)]*5
        }
        self.variable_ranges = variable_ranges_map[self.dimension]
        self.variable_range = (0.5, 3)
        self.input_file_name = 'ThreePointBending_0000.rad'
        self.output_file_name = 'ThreePointBendingT01.csv'
        self.starter_out_file_name = 'ThreePointBending_0000.out'
        self.track_node_key = 'intrusionTrack44721' # the key of the intrusion in the output csv file

        self.problem_id = ThreePointBending.instance_counter
        ThreePointBending.instance_counter+=1

    def generate_input_deck(self, variable_array):
        '''
        Special generate_input_deck function for Three point bending model, transform it into thickness mapping.
        '''
        self._validate_variable_array(variable_array)
        
        thickness_array = [] # to get the thickness in the FEM space
        if len(variable_array)==1:
            variable_array = [variable_array[0]]*5
        elif len(variable_array)==2:
            variable_array = [variable_array[0],0,0,0,variable_array[1]]
        elif len(variable_array)==3:
            variable_array = [variable_array[0],0,variable_array[1],0,variable_array[2]]
        elif len(variable_array)==4:
            variable_array = [variable_array[0],variable_array[1],0,variable_array[2],variable_array[3]]

        for i, var in enumerate(variable_array):
            mapped_var = self.linear_maping_variable(var, self.variable_range)
            thickness_array.append(mapped_var)

        original_dir = os.getcwd()
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        os.chdir(working_dir)
        self._write_input_file(thickness_array)
        os.chdir(original_dir)

        self.problem_id += 1 # update the problem id for the input deck generation

    def _copy_files_to_deck(self):
        # Get the path to the lib directory relative to the current file
        deck_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'

        lib_dir = os.path.join(os.path.dirname(__file__), 'lib')
        
        # Source file paths
        source_files = ['ThreePointBending_base.rad', 'ThreePointBending_0001.rad']  # Copy base starter file and engine files
        # Destination folder path
        dest_folder = os.getcwd()
        
        # Create the destination folder if it doesn't exist
        os.makedirs(dest_folder, exist_ok=True)
        
        # Copy each source file to the destination folder
        for file_name in source_files:
            source_path = os.path.join(lib_dir, file_name)
            dest_path = os.path.join(dest_folder, file_name)
            shutil.copyfile(source_path, dest_path)

    def _write_input_file(self, thickness_array):
        self._copy_files_to_deck()
        self.model = ThreePointBendingModel()
        self.model.write_input_file(thickness_array)


class CrashTube(OptiProblem):
    instance_counter = 1
    def __init__(self, dimension, output_data, batch_file_path) -> None:
        super().__init__(dimension, output_data, batch_file_path)
        # 2 -> three positions and depths vary together with same value. 
        # 3 -> all three trigger positions fixed, the three trigger depth vary. 
        # 4 -> except for middle trigger position and depth, other 4 vary ().
        # 5 -> except for middle trigger depth, other 5 vary. 
        # 6 -> all three positions and depths vary.
        variable_ranges_map = {
            2: [(-10, 10), (-4, 4)],
            3: [(-4, 4)]*3,
            4: [(-10, 10), (-4, 4), (-10, 10), (-4, 4)],
            5: [(-10, 10), (-10, 10), (-10, 10), (-4, 4), (-4, 4)],
            6: [(-10, 10), (-10, 10), (-10, 10), (-4, 4), (-4, 4), (-4, 4)]
        }
        self.variable_ranges = variable_ranges_map[self.dimension]
        self.input_file_name = 'combine.k'
        self.output_file_name = 'combineT01.csv'
        self.starter_out_file_name = 'combine_0000.out'
        self.track_node_key = 'DATABASE_HISTORY_NODE1001' # the key of the intrusion in the output csv file

        self.problem_id = CrashTube.instance_counter
        CrashTube.instance_counter+=1

    def _write_input_file(self, fem_space_variable_array):
        self.mesh = CrashTubeMesh(fem_space_variable_array) 
        self.model = CrashTubeModel(self.mesh)
        self.model.write_input_files()