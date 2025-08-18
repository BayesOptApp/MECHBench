import os
import pandas as pd
from abc import ABC, abstractmethod
from typing import List, Iterable, Union
import shutil
from .solver import run_radioss
from .utils.solver_setup import RunnerOptions
from .mesh import *
from .fem import *




### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### THESE ARE CONSTANTS TO DETERMINE THE OUTPUTS (OBJECTIVES) OF THE PROBLEM
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

__PRINCIPAL_OUTPUTS:tuple = ("mass", "absorbed_energy","intrusion","mean_impact_force","max_impact_force")
__COMPOSITE_OUTPUTS:tuple = ("specific_energy","usage_ratio","load_uniformity")

class OptiProblem(ABC):
    r'''
    The abstract of the structural optimization problem, with three subclass (three problem types)
    The instance is initialized by problem dimension, output data type, and the batch file path of the solver OpenRadioss. 
    The core idea is to generate a specific type of instance, input the variable array then output the evaluation of the desired data type like mass, intrusion, etc.
    The input variables should be in an universal search space, default (-5,5). For the evaluation those variables will be mapped to the real FEM problem space.
    '''
    def __init__(self, dimension:int, 
                 output_data:Union[Iterable,str], 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool=True) -> None:
        
        
        self.__dimension:int = dimension
        self.__output_data:Union[Iterable,str] = output_data
        self.__sequential_id_numbering:bool = sequential_id_numbering

        # Assign the Properties to the case
        self._runner_options = RunnerOptions(runner_options)

        self._runner_options.complement()
        self._runner_options.amend_integer_options(runner_options)


        # The attributes need to be overwritten in teh subclass
        self._problem_id:int = -1 
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

        if self.__sequential_id_numbering:
            self.problem_id += 1 # update the problem id for the input deck generation
    
    @abstractmethod
    def _write_input_file(self, fem_space_variable_array):
        pass
    
    def run_simulation(self,runStarter=False):
        if self.input_file_name is None:
            raise ValueError("input_file_name must be provided or defined in the subclass.")
        # make problem id back to original, since it has been updated when generate_input_deck has been called
        if self.__sequential_id_numbering:
            self.problem_id -= 1

        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        input_file_path = os.path.join(working_dir, self.input_file_name)

        if runStarter:
            # This is just one bypass in order to avoid setting MP settings
            # for just the mass computation
            run_radioss(input_file_path, 
                        self.batch_file_path, 
                        runStarter=runStarter,
                        write_vtk=False,
                        np_int=1,
                        nt_int=1)
        else:
            run_radioss(input_file_path, 
                        self.batch_file_path, 
                        runStarter=runStarter,
                        write_vtk=bool(self._runner_options('write_vtk')),
                        np_int=self._runner_options('np'),
                        nt_int=self._runner_options('nt'))

    def load_output_data_frame(self):
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        # load simulation result dataframe and make it cleaner
        output_file_path = os.path.join(working_dir, self.output_file_name)
        self.output_data_frame = pd.read_csv(output_file_path)
        self.output_data_frame.columns = self.output_data_frame.columns.str.replace(' ', '')
        
        # update problem id again
        if self.__sequential_id_numbering:
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


    def instrusion_calculation(self)->float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.load_output_data_frame()
            self.sim_status = 2

        col = [col for col in self.output_data_frame.columns if self.track_node_key in col][2]
 

        return abs(self.output_data_frame[col].abs().max()) - self.model.impactor_offset
    
    def mass_calculation(self)->float:
        if self.sim_status < 1: 
            self.run_simulation(runStarter=True)
            # Change the status
            self.sim_status = 1

        val:float = self.extract_mass_from_file() - self.model.rigid_mass

        if self.__sequential_id_numbering:
            self.problem_id +=1  

        return val
    
    def absorbed_energy_calculation(self)->float:
        return self.model.absorbed_energy()
    
    def _get_max_intrusion_index(self):
        col = [col for col in self.output_data_frame.columns if self.track_node_key in col][2]
        return self.output_data_frame[col].abs().idxmax()

    def _get_force_data(self):
        force_col = [col for col in self.output_data_frame.columns if self.impactor_force_key in col][2]
        max_idx = self._get_max_intrusion_index()
        df = self.output_data_frame.loc[:max_idx]
        
        if self.model.material_card_type == "OpenRadioss":
            return df[force_col].values
        else:
            time_vect = df["time"].values
            impulse_vec = df[force_col].values
            return np.gradient(impulse_vec, time_vect)

    def peak_force_calculation(self) -> float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.sim_status = 2
            self.load_output_data_frame()
        force_data = self._get_force_data()
        
        return np.abs(np.max(force_data).astype(float))

    def mean_force_calculation(self) -> float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.load_output_data_frame()
            self.sim_status = 2


        force_data = self._get_force_data()
        
        return np.abs(np.mean(force_data).astype(float))
    


    def __call__(self, variable_array: list, problem_id: int = -1) -> Union[float, List[float]]:
        
        if not self.__sequential_id_numbering:
            self.problem_id = problem_id
        
        # Set the sim status to 0
        self.sim_status = 0

        # Generate the input deck
        self.generate_input_deck(variable_array)

        def handle_single(key: str) -> float:
            return {
                'mass': self.mass_calculation,
                'absorbed_energy': self.absorbed_energy_calculation,
                'intrusion': self.instrusion_calculation,
                'mean_impact_force': self.mean_force_calculation,
                'max_impact_force': self.peak_force_calculation,
                'specific_energy': lambda: self.absorbed_energy_calculation() / self.mass_calculation(),
                'load_uniformity': lambda: abs(self.peak_force_calculation() / self.mean_force_calculation())
            }.get(key, lambda: np.nan)()

        if isinstance(self.output_data, str):
            return handle_single(self.output_data)

        elif isinstance(self.output_data, list):
            result = []
            intrusion_index = None
            output_keys = self.output_data.copy()

            if 'intrusion' in output_keys:
                intrusion_index = output_keys.index('intrusion')
                output_keys.pop(intrusion_index)

            for key in output_keys:
                result.append(handle_single(key))

            if intrusion_index is not None:
                result.insert(intrusion_index, self.instrusion_calculation())

            return result

        return None

    @property
    def problem_id(self)->int:
        return self._problem_id
    
    @property
    def dimension(self)->int:
        return self.__dimension
    
    @property
    def batch_file_path(self)->str:
        return self._runner_options('open_radioss_main_path')
    
    @property
    def output_data(self)->Union[List[str],str]:
        return self.__output_data
    
    @property
    def sequential_id_numbering(self):
        return self.__sequential_id_numbering
    

    @problem_id.setter
    def problem_id(self,new_problem_id:int)->None:
        
        if not isinstance(new_problem_id,int):
            raise TypeError("The type of the object is not correct")
        
        else:
            if new_problem_id <= 0:
                raise ValueError("The value must be greater than 0")
            # Assign the new problem ID
            self._problem_id = new_problem_id
    

    @output_data.setter
    def output_data(self, new_output_data:Union[Iterable,str])->None:
        #self.__output_data = new_output_data

        # CASE 1: Check the output data is a string
        if isinstance(new_output_data,str):
            if ((new_output_data.strip().lower() in __COMPOSITE_OUTPUTS) or 
                (new_output_data.strip().lower() in __PRINCIPAL_OUTPUTS)):
                # Assign
                self.__output_data = new_output_data
            else:
                raise ValueError("The output data is not defined as an output")
        elif isinstance(new_output_data,(list,List[str])):
            # Loop all over the elements and check the elements are part 
            # of the defined group of outputs

            idx_to_remove = []
            for idx, elem in enumerate(new_output_data):

                try:
                    if not ((elem.strip().lower() in __COMPOSITE_OUTPUTS) or 
                    (elem.strip().lower() in __PRINCIPAL_OUTPUTS)):
                        idx_to_remove.append(idx)
                except Exception as e:
                    print("The list of output data is not correctly set as a list of strings")
            

            for idx in idx_to_remove:
                try:
                    new_output_data.pop(idx)
                except IndexError as e:
                    print("The output data is empty")
                except Exception as e:
                    print("Something else happened...")
                
            
            # Now try to set the output data
            self.__output_data = new_output_data
                    
    @property
    def model(self)->AbstractModel:
        return self.__model
    
    @model.setter
    def model(self, new_model:AbstractModel)->None:
        if not issubclass(type(new_model),AbstractModel) and (new_model is not None):
            raise TypeError("The type of the object is not correct")
        else:
            self.__model = new_model

    
    @output_data.deleter
    def output_data(self)->None:
        del self.__output_data




class StarBox(OptiProblem):
    instance_counter = 1

    def __init__(self, 
                 dimension, 
                 output_data, 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool,**kwargs) -> None:
        
        super().__init__(dimension, 
                         output_data, 
                         runner_options,
                         sequential_id_numbering)
        # 1 -> square
        # 2 -> rectangular
        # 3 -> rectangular with varying thickness
        # 4 -> star shape
        # 5 -> star shape with varying thickness
        # 6 - 35 -> star shape with different thickness profiles.
        
        self.variable_ranges = self._generate_variable_ranges_map(self.dimension)
        self.input_file_name = 'combine.k'
        self.output_file_name = 'combineT01.csv'
        self.starter_out_file_name = 'combine_0000.out'

        # the key of the intrusion in the output csv file
        #self.track_node_key = 'DATABASE_HISTORY_NODE1000'
        self.track_node_key = 'DATABASE_HISTORY_NODE99999' 
        self.impactor_force_key = "TH-RWALL1"

        if self.sequential_id_numbering:
            self.problem_id = StarBox.instance_counter
            StarBox.instance_counter+=1
    
    @staticmethod
    def _generate_variable_ranges_map(dimension:int)->List[tuple]:
        r"""
        Generates the ranges map given a dimensionality; This is a 
        static method, which is provided with the class to call methods
        outsude the framework.

        Args
        ----------------------
        - dimension: `int`: An integer denoting the dimensionality set by the user

        Returns
        ----------------------
        - `List[tuple]`: A list of tuples indicating the physical ranges of the problem.
        """

        variable_ranges_map = {
        1: [(60, 120)],
        2: [(60, 120), (60, 120)],
        3: [(60, 120), (60, 120), (0.7, 3)],
        4: [(60, 120), (60, 120), (0, 30), (0, 30)],
        5: [(60, 120), (60, 120), (0, 30), (0, 30), (0.7, 3)] 
                        }
        
        if dimension > 0 and dimension<=5:
         
            return variable_ranges_map[dimension]
        
        else:
            base_ranges:List[tuple] = variable_ranges_map[5]

            for _ in range(dimension-5):
                base_ranges.append((0.7,3))
            
            return base_ranges


    def _write_input_file(self, fem_space_variable_array):
        self.mesh = StarBoxMesh(fem_space_variable_array,
                                h_level=self._runner_options('h_level'),
                                gmsh_verbosity=self._runner_options('gmsh_verbosity')
        )
        self.model = StarBoxModel(self.mesh)
        self.model.write_input_files()
    

class ThreePointBending(OptiProblem):
    '''
    This optimization problem focuses on adjusting thickness, of 5 
    sheets of material, to optimize the performance of a three-point
    bending test. The design variables are the thickness of the sheets,
    which can be varied within a specified range. The objective is to
    minimize the load uniformity and maximize the absorbed energy during
    the bending test. The problem is defined by a set of constraints and
    a search space for the design variables. The optimization process
    involves generating input decks for the simulation, running the
    simulations, and analyzing the results to find the optimal thickness
    configuration for the sheets.

    The design variables are the thickness of the sheets, which can be
    varied within a specified range. 
    '''
    instance_counter = 1
    def __init__(self, 
                 dimension, 
                 output_data, 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool) -> None:
        super().__init__(dimension, output_data, runner_options,sequential_id_numbering)
        # 1 -> all 5 shell thickness vary with same value. 
        # 2 -> only first and the last shell thickness vary, other fixed with middle value. 
        # 3 -> the first, middle, and last shell thickness vary. 
        # 4 -> expect for middle shell, the other 4 shell thickness vary.
        # 5 -> all three shell thickness vary.
        variable_ranges_map = { ii + 1 : [(-5, 5)]*ii for ii in range(40)}
        
        self.variable_ranges = variable_ranges_map[self.dimension]
        self.variable_range = (0.5, 3)
        self.input_file_name = 'ThreePointBending_0000.rad'
        self.output_file_name = 'ThreePointBendingT01.csv'
        self.starter_out_file_name = 'ThreePointBending_0000.out'
        self.track_node_key = 'intrusionTrack99999' # the key of the intrusion in the output csv file
        self.impactor_force_key = "TH_RWALL1"

        if self.sequential_id_numbering:
            self.problem_id = ThreePointBending.instance_counter
            ThreePointBending.instance_counter+=1

    def generate_input_deck(self, variable_array):
        '''
        Special generate_input_deck function for Three point bending model, transform it into thickness mapping.
        '''
        self._validate_variable_array(variable_array)
        
        thickness_array = [] # to get the thickness in the FEM space

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

        if self.sequential_id_numbering:
            self.problem_id = ThreePointBending.instance_counter
            ThreePointBending.instance_counter+=1

    def _copy_files_to_deck(self):
        # Get the path to the lib directory relative to the current file
        deck_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'

        lib_dir = os.path.join(os.path.dirname(__file__), 'lib')
        
        # Source file paths
        source_files = ['ThreePointBending_0001.rad']  # Copy base starter file and engine files
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
        self.mesh = ThreePointBendingMesh(thickness_array,
                                          h_level=self._runner_options('h_level'),
                                          gmsh_verbosity=self._runner_options('gmsh_verbosity')
                                          )
        
        self.model = ThreePointBendingModel(self.mesh)
        self.model.write_input_file(thickness_array)
    
    def mass_calculation(self):
        r"""
        This is a refactoring of the mass calculation function, which is
        used to calculate the mass of the model. The function is called
        when the simulation is run, and it returns the mass of the model.

        Since the mass is computed in metric tons, then this function
        converts the mass to kg. 
        """
        return super().mass_calculation()*1000.0


class CrashTube(OptiProblem):
    instance_counter = 1
    def __init__(self, dimension, output_data, 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool) -> None:
        super().__init__(dimension, output_data, runner_options ,sequential_id_numbering)

        ### NOTE: THIS IS THE PREVIOUS DEFINITION
        ### NOTE: JUST WRITTEN FOR COMPARISON 
        # 2 -> three positions and depths vary together with same value. 
        # 3 -> all three trigger positions fixed, the three trigger depth vary. 
        # 4 -> except for middle trigger position and depth, other 4 vary ().
        # 5 -> except for middle trigger depth, other 5 vary. 
        # 6 -> all three positions and depths vary.
        # variable_ranges_map = {
        #     2: [(-10, 10), (-4, 4)],
        #     3: [(-4, 4)]*3,
        #     4: [(-10, 10), (-4, 4), (-10, 10), (-4, 4)],
        #     5: [(-10, 10), (-10, 10), (-10, 10), (-4, 4), (-4, 4)],
        #     6: [(-10, 10), (-10, 10), (-10, 10), (-4, 4), (-4, 4), (-4, 4)]
        # }

        ### # NOTE: THIS IS THE NEW DEFINITION
        r"""
        - Dimensionality : 1 - 3 -> 1 Trigger on the 4 sides with 1-3 degrees of freedom (Priority: Vertical Position; Depth; height)
        - Dimensionality: 4 - 6 Two triggers
        - Dimensionality: 7 - 9 Three triggers
        - Dimensionality: 10 - 12 Four triggers
        - Dimensionality: 13 - 15 Five triggers
        - Dimensionality: 16 - 18 Three independent variables for neighbouring faces of first segment
        - Dimensionality: 19 - 21 Three independent variables for neighbouring faces of first & second segment
        - Dimensionality: 21 - 24 Three independent variables for neighbouring faces of first, second and third segment
        - Dimensionality: 25 - 27 Three independent variables for neighbouring faces of first, second, third and fourth
                                 segment
        - Dimensionality: 28 - 30 Three independent variables for neighbouring faces of all segments
        ...
        """

        ### # NOTE: END OF NEW DEFINITION
        
        self.variable_ranges = self._generate_variable_ranges_map(dimension)
        self.input_file_name = 'combine.k'
        self.output_file_name = 'combineT01.csv'
        self.starter_out_file_name = 'combine_0000.out'
        self.track_node_key = 'DATABASE_HISTORY_NODE99999' # the key of the intrusion in the output csv file
        self.impactor_force_key = "TH-RWALL1IMPACTOR"

        if self.sequential_id_numbering:
            self.problem_id = CrashTube.instance_counter
            CrashTube.instance_counter+=1

    def _write_input_file(self, fem_space_variable_array):
        self.mesh = CrashTubeMesh(fem_space_variable_array,
                                  h_level=self._runner_options('h_level'),
                                  gmsh_verbosity=self._runner_options('gmsh_verbosity')
                                  ) 
        self.model = CrashTubeModel(self.mesh)
        self.model.write_input_files()
    

    @staticmethod
    def _generate_variable_ranges_map(dimension:int)->List[tuple]:
        r"""
        Generates the ranges map given a dimensionality; This is a 
        static method, which is provided with the class to call methods
        outsude the framework.

        Args
        ----------------------
        - dimension: `int`: An integer denoting the dimensionality set by the user

        Returns
        ----------------------
        - `List[tuple]`: A list of tuples indicating the physical ranges of the problem.
        """

        patterns = [(-10, 10), (-4, 4), (0, 4)]

        variable_ranges_map = {}

        for dim in range(1, dimension+1):
            variable_ranges_map[dim] = [patterns[i % 3] for i in range(dim)]
        
        return variable_ranges_map[dimension]
       
import os
import pandas as pd
from abc import ABC, abstractmethod
from typing import List, Iterable, Union
import shutil
from .solver import run_radioss
from .utils.solver_setup import RunnerOptions
from .mesh import *
from .fem import *




### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### THESE ARE CONSTANTS TO DETERMINE THE OUTPUTS (OBJECTIVES) OF THE PROBLEM
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_PRINCIPAL_OUTPUTS:tuple = ("mass", "absorbed_energy","intrusion","mean_impact_force","max_impact_force")
_COMPOSITE_OUTPUTS:tuple = ("specific_energy_absorbed","usage_ratio","load_uniformity")
_PROBLEM_SPECIFIC_OUTPUTS:tuple = ("penalized_sea","penalized_mass")

class OptiProblem(ABC):
    r'''
    The abstract of the structural optimization problem, with three subclass (three problem types)
    The instance is initialized by problem dimension, output data type, and the batch file path of the solver OpenRadioss. 
    The core idea is to generate a specific type of instance, input the variable array then output the evaluation of the desired data type like mass, intrusion, etc.
    The input variables should be in an universal search space, default (-5,5). For the evaluation those variables will be mapped to the real FEM problem space.
    '''
    def __init__(self, dimension:int, 
                 output_data:Union[Iterable,str], 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool=True) -> None:
        
        
        self.dimension:int = dimension
        self.output_data:Union[Iterable,str] = output_data
        self.__sequential_id_numbering:bool = sequential_id_numbering

        # Assign the Properties to the case
        self._runner_options = RunnerOptions(runner_options)

        self._runner_options.complement()
        self._runner_options.amend_integer_options(runner_options)


        # The attributes need to be overwritten in teh subclass
        self._problem_id:int = -1 
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

        if self.__sequential_id_numbering:
            self.problem_id += 1 # update the problem id for the input deck generation
    
    @abstractmethod
    def _write_input_file(self, fem_space_variable_array):
        pass
    
    def run_simulation(self,runStarter=False):
        if self.input_file_name is None:
            raise ValueError("input_file_name must be provided or defined in the subclass.")
        # make problem id back to original, since it has been updated when generate_input_deck has been called
        if self.__sequential_id_numbering:
            self.problem_id -= 1

        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        input_file_path = os.path.join(working_dir, self.input_file_name)

        if runStarter:
            # This is just one bypass in order to avoid setting MP settings
            # for just the mass computation
            run_radioss(input_file_path, 
                        self.batch_file_path, 
                        runStarter=runStarter,
                        write_vtk=False,
                        np_int=1,
                        nt_int=1)
        else:
            run_radioss(input_file_path, 
                        self.batch_file_path, 
                        runStarter=runStarter,
                        write_vtk=bool(self._runner_options('write_vtk')),
                        np_int=self._runner_options('np'),
                        nt_int=self._runner_options('nt'))

    def load_output_data_frame(self):
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        # load simulation result dataframe and make it cleaner
        output_file_path = os.path.join(working_dir, self.output_file_name)
        self.output_data_frame = pd.read_csv(output_file_path)
        self.output_data_frame.columns = self.output_data_frame.columns.str.replace(' ', '')
        
        # update problem id again
        if self.__sequential_id_numbering:
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


    def instrusion_calculation(self)->float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.load_output_data_frame()
            self.sim_status = 2

        col = [col for col in self.output_data_frame.columns if self.track_node_key in col][2]
 

        return abs(self.output_data_frame[col].abs().max()) - self.model.impactor_offset
    
    def mass_calculation(self)->float:
        if self.sim_status < 1: 
            self.run_simulation(runStarter=True)
            # Change the status
            self.sim_status = 1

        val:float = self.extract_mass_from_file() - self.model.rigid_mass

        if self.__sequential_id_numbering:
            self.problem_id +=1  

        return val
    
    def absorbed_energy_calculation(self)->float:
        return self.model.absorbed_energy()
    
    def _get_max_intrusion_index(self):
        col = [col for col in self.output_data_frame.columns if self.track_node_key in col][2]
        return self.output_data_frame[col].abs().idxmax()

    def _get_force_data(self):
        force_col = [col for col in self.output_data_frame.columns if self.impactor_force_key in col][2]
        max_idx = self._get_max_intrusion_index()
        df = self.output_data_frame.loc[:max_idx]
        
        if self.model.material_card_type == "OpenRadioss":
            return df[force_col].values
        else:
            time_vect = df["time"].values
            impulse_vec = df[force_col].values
            return np.gradient(impulse_vec, time_vect)

    def peak_force_calculation(self) -> float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.sim_status = 2
            self.load_output_data_frame()
        force_data = self._get_force_data()
        
        return np.abs(np.max(force_data).astype(float))

    def mean_force_calculation(self) -> float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.load_output_data_frame()
            self.sim_status = 2


        force_data = self._get_force_data()
        
        return np.abs(np.mean(force_data).astype(float))
    


    def __call__(self, variable_array: list, problem_id: int = -1) -> Union[float, List[float]]:
        
        if not self.__sequential_id_numbering:
            self.problem_id = problem_id
        
        # Set the sim status to 0
        self.sim_status = 0

        # Generate the input deck
        self.generate_input_deck(variable_array)

        def handle_single(key: str) -> float:
            return {

                'mass': self.mass_calculation,

                'absorbed_energy': self.absorbed_energy_calculation,

                'intrusion': self.instrusion_calculation,

                'mean_impact_force': self.mean_force_calculation,

                'max_impact_force': self.peak_force_calculation,

                'specific_energy_absorbed': lambda: self.absorbed_energy_calculation() / self.mass_calculation(),

                'load_uniformity': lambda: abs(self.peak_force_calculation() / self.mean_force_calculation()),

                'penalized_sea': lambda: -(

                    self.absorbed_energy_calculation() / self.mass_calculation()

                    if self.instrusion_calculation() <= 60

                    else -100 * (self.instrusion_calculation() - 60)),

                'penalized_mass': lambda: (self.mass_calculation()

                                           if self.instrusion_calculation() <= 50

                                           else 4.25952 + 10 * (self.instrusion_calculation() / 50 - 1))



            }.get(key, lambda: np.nan)()

        if isinstance(self.output_data, str):
            return handle_single(self.output_data)

        elif isinstance(self.output_data, list):
            result = []
            intrusion_index = None
            output_keys = self.output_data.copy()

            if 'intrusion' in output_keys:
                intrusion_index = output_keys.index('intrusion')
                output_keys.pop(intrusion_index)

            for key in output_keys:
                result.append(handle_single(key))

            if intrusion_index is not None:
                result.insert(intrusion_index, self.instrusion_calculation())

            return result

        return None

    @property
    def problem_id(self)->int:
        return self._problem_id
    
    @property
    def dimension(self)->int:
        return self.__dimension
    
    @dimension.setter
    def dimension(self, new_dimension:int)->None:
        if not isinstance(new_dimension,int):
            raise TypeError("The type of the object is not correct")
        
        else:
            if new_dimension <= 0:
                raise ValueError("The value must be greater than 0")
            # Assign the new dimension
            self.__dimension = new_dimension
    
    @property
    def batch_file_path(self)->str:
        return self._runner_options('open_radioss_main_path')
    
    @property
    def output_data(self)->Union[List[str],str]:
        return self._output_data
    
    @property
    def sequential_id_numbering(self):
        return self.__sequential_id_numbering
    

    @problem_id.setter
    def problem_id(self,new_problem_id:int)->None:
        
        if not isinstance(new_problem_id,int):
            raise TypeError("The type of the object is not correct")
        
        else:
            if new_problem_id <= 0:
                raise ValueError("The value must be greater than 0")
            # Assign the new problem ID
            self._problem_id = new_problem_id
    

    @output_data.setter
    def output_data(self, new_output_data:Union[Iterable,str])->None:
        #self.__output_data = new_output_data

        # CASE 1: Check the output data is a string
        if isinstance(new_output_data,str):
            if ((new_output_data.strip().lower() in _COMPOSITE_OUTPUTS) or 
                (new_output_data.strip().lower() in _PRINCIPAL_OUTPUTS) or
                (new_output_data.strip().lower() in _PROBLEM_SPECIFIC_OUTPUTS)):
                # Assign
                self._output_data = new_output_data
            else:
                raise ValueError("The output data is not defined as an output")
        elif isinstance(new_output_data,(list,List[str])):
            # Loop all over the elements and check the elements are part 
            # of the defined group of outputs

            idx_to_remove = []
            for idx, elem in enumerate(new_output_data):

                try:
                    if not ((elem.strip().lower() in _COMPOSITE_OUTPUTS) or 
                    (elem.strip().lower() in _PRINCIPAL_OUTPUTS) or 
                    (elem.strip().lower() in _PROBLEM_SPECIFIC_OUTPUTS)):
                        idx_to_remove.append(idx)
                except Exception as e:
                    print("The list of output data is not correctly set as a list of strings")
            

            for idx in idx_to_remove:
                try:
                    new_output_data.pop(idx)
                except IndexError as e:
                    print("The output data is empty")
                except Exception as e:
                    print("Something else happened...")
                
            
            # Now try to set the output data
            self._output_data = new_output_data
                    
    @property
    def model(self)->AbstractModel:
        return self.__model
    
    @model.setter
    def model(self, new_model:AbstractModel)->None:
        if not issubclass(type(new_model),AbstractModel) and (new_model is not None):
            raise TypeError("The type of the object is not correct")
        else:
            self.__model = new_model

    
    @output_data.deleter
    def output_data(self)->None:
        del self._output_data




class StarBox(OptiProblem):
    instance_counter = 1

    def __init__(self, 
                 dimension, 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool,
                 output_data:Optional[Union[Iterable,str]],
                 **kwargs) -> None:
        r"""
        The StarBox optimization problem is a structural optimization problem
        that involves a star-shaped box with varying thicknesses and the dimensions of the star shape.
        """

        # Get the output data if it is not provided
        if output_data is None:
            output_data = "penalized_sea"
        elif isinstance(output_data, str):
            if output_data in self.forbidden_output_data:
                raise ValueError(f"The output data {output_data} is not allowed for the StarBox problem.")
        elif isinstance(output_data, tuple) or isinstance(output_data, list):
            for elem in output_data:
                if elem in self.forbidden_output_data:
                    raise ValueError(f"The output data {elem} is not allowed for the StarBox problem.")


        
        super().__init__(dimension, 
                         output_data, 
                         runner_options,
                         sequential_id_numbering)
        # 1 -> square
        # 2 -> rectangular
        # 3 -> rectangular with varying thickness
        # 4 -> star shape
        # 5 -> star shape with varying thickness
        # 6 - 34 -> star shape with different thickness profiles.
        
        self.variable_ranges = self._generate_variable_ranges_map(self.dimension)
        self.input_file_name = 'combine.k'
        self.output_file_name = 'combineT01.csv'
        self.starter_out_file_name = 'combine_0000.out'

        # the key of the intrusion in the output csv file
        self.track_node_key = 'DATABASE_HISTORY_NODE99999' 
        self.impactor_force_key = "TH-RWALL1"

        if self.sequential_id_numbering:
            self.problem_id = StarBox.instance_counter
            StarBox.instance_counter+=1
    
    @staticmethod
    def _generate_variable_ranges_map(dimension:int)->List[tuple]:
        r"""
        Generates the ranges map given a dimensionality; This is a 
        static method, which is provided with the class to call methods
        outsude the framework.

        Args
        ----------------------
        - dimension: `int`: An integer denoting the dimensionality set by the user

        Returns
        ----------------------
        - `List[tuple]`: A list of tuples indicating the physical ranges of the problem.
        """

        variable_ranges_map = {
        1: [(60, 120)],
        2: [(60, 120), (60, 120)],
        3: [(60, 120), (60, 120), (0.7, 3)],
        4: [(60, 120), (60, 120), (0, 30), (0, 30)],
        5: [(60, 120), (60, 120), (0, 30), (0, 30), (0.7, 3)] 
                        }
        
        if dimension > 0 and dimension<=5:
         
            return variable_ranges_map[dimension]
        
        else:
            base_ranges:List[tuple] = variable_ranges_map[5]

            for _ in range(dimension-5):
                base_ranges.append((0.7,3))
            
            return base_ranges


    def _write_input_file(self, fem_space_variable_array)->None:
        '''
        Write the input files for the StarBox problem to be processed by OpenRadioss.
        Args:
            fem_space_variable_array (list): The array of variables in the actual physical space.
        Returns:
            None
        '''
        self.mesh = StarBoxMesh(fem_space_variable_array,
                                h_level=self._runner_options('h_level'),
                                gmsh_verbosity=self._runner_options('gmsh_verbosity')
        )
        self.model = StarBoxModel(self.mesh)
        self.model.write_input_files()
    
    @property
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the StarBox problem.
        """
        return ['penalized_mass']
    
    
    
    

class ThreePointBending(OptiProblem):
    '''
    This optimization problem focuses on adjusting thickness, of 5 
    sheets of material, to optimize the performance of a three-point
    bending test. The design variables are the thickness of the sheets,
    which can be varied within a specified range. The objective is to
    minimize the load uniformity and maximize the absorbed energy during
    the bending test. The problem is defined by a set of constraints and
    a search space for the design variables. The optimization process
    involves generating input decks for the simulation, running the
    simulations, and analyzing the results to find the optimal thickness
    configuration for the sheets.

    The design variables are the thickness of the sheets, which can be
    varied within a specified range. 
    '''
    instance_counter = 1
    def __init__(self, 
                 dimension, 
                 output_data, 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool) -> None:
        # Get the output data if it is not provided
        if output_data is None:
            output_data = "penalized_mass"
        elif isinstance(output_data, str):
            if output_data in self.forbidden_output_data:
                raise ValueError(f"The output data {output_data} is not allowed for the StarBox problem.")
        elif isinstance(output_data, tuple) or isinstance(output_data, list):
            for elem in output_data:
                if elem in self.forbidden_output_data:
                    raise ValueError(f"The output data {elem} is not allowed for the StarBox problem.")

        
        super().__init__(dimension, output_data, runner_options,sequential_id_numbering)
        # 1 -> all 5 shell thickness vary with same value. 
        # 2 -> only first and the last shell thickness vary, other fixed with middle value. 
        # 3 -> the first, middle, and last shell thickness vary. 
        # 4 -> expect for middle shell, the other 4 shell thickness vary.
        # 5 -> all three shell thickness vary.
        variable_ranges_map = { ii + 1 : [(-5, 5)]*ii for ii in range(40)}
        
        self.variable_ranges = variable_ranges_map[self.dimension]
        self.variable_range = (0.5, 3)
        self.input_file_name = 'ThreePointBending_0000.rad'
        self.output_file_name = 'ThreePointBendingT01.csv'
        self.starter_out_file_name = 'ThreePointBending_0000.out'
        self.track_node_key = 'intrusionTrack99999' # the key of the intrusion in the output csv file
        self.impactor_force_key = "TH_RWALL1"

        if self.sequential_id_numbering:
            self.problem_id = ThreePointBending.instance_counter
            ThreePointBending.instance_counter+=1

    def generate_input_deck(self, variable_array):
        '''
        Special generate_input_deck function for Three point bending model, transform it into thickness mapping.
        '''
        self._validate_variable_array(variable_array)
        
        thickness_array = [] # to get the thickness in the FEM space

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

        if self.sequential_id_numbering:
            self.problem_id = ThreePointBending.instance_counter
            ThreePointBending.instance_counter+=1

    def _copy_files_to_deck(self):
        # Get the path to the lib directory relative to the current file
        deck_name = f'{self.__class__.__name__.lower()}_deck{self.problem_id}'

        lib_dir = os.path.join(os.path.dirname(__file__), 'lib')
        
        # Source file paths
        source_files = ['ThreePointBending_0001.rad']  # Copy base starter file and engine files
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
        self.mesh = ThreePointBendingMesh(thickness_array,
                                          h_level=self._runner_options('h_level'),
                                          gmsh_verbosity=self._runner_options('gmsh_verbosity')
                                          )
        
        self.model = ThreePointBendingModel(self.mesh)
        self.model.write_input_file(thickness_array)
    
    def mass_calculation(self):
        r"""
        This is a refactoring of the mass calculation function, which is
        used to calculate the mass of the model. The function is called
        when the simulation is run, and it returns the mass of the model.

        Since the mass is computed in metric tons, then this function
        converts the mass to kg. 
        """
        return super().mass_calculation()*1000.0
    
    @property
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the StarBox problem.
        """
        return ['penalized_sea']



class CrashTube(OptiProblem):
    instance_counter = 1
    def __init__(self, dimension, output_data, 
                 runner_options:RunnerOptions,
                 sequential_id_numbering:bool) -> None:
        
        # Get the output data if it is not provided
        if output_data is None:
            output_data = "load_uniformity"
        elif isinstance(output_data, str):
            if output_data in self.forbidden_output_data:
                raise ValueError(f"The output data {output_data} is not allowed for the StarBox problem.")
        elif isinstance(output_data, tuple) or isinstance(output_data, list):
            for elem in output_data:
                if elem in self.forbidden_output_data:
                    raise ValueError(f"The output data {elem} is not allowed for the StarBox problem.")
 
        
        super().__init__(dimension=dimension, 
                         output_data=output_data, 
                         runner_options=runner_options ,
                         sequential_id_numbering=sequential_id_numbering)

        ### NOTE: THIS IS THE PREVIOUS DEFINITION
        ### NOTE: JUST WRITTEN FOR COMPARISON 
        # 2 -> three positions and depths vary together with same value. 
        # 3 -> all three trigger positions fixed, the three trigger depth vary. 
        # 4 -> except for middle trigger position and depth, other 4 vary ().
        # 5 -> except for middle trigger depth, other 5 vary. 
        # 6 -> all three positions and depths vary.
        # variable_ranges_map = {
        #     2: [(-10, 10), (-4, 4)],
        #     3: [(-4, 4)]*3,
        #     4: [(-10, 10), (-4, 4), (-10, 10), (-4, 4)],
        #     5: [(-10, 10), (-10, 10), (-10, 10), (-4, 4), (-4, 4)],
        #     6: [(-10, 10), (-10, 10), (-10, 10), (-4, 4), (-4, 4), (-4, 4)]
        # }

        ### # NOTE: THIS IS THE NEW DEFINITION
        r"""
        - Dimensionality : 1 - 3 -> 1 Trigger on the 4 sides with 1-3 degrees of freedom (Priority: Vertical Position; Depth; height)
        - Dimensionality: 4 - 6 Two triggers
        - Dimensionality: 7 - 9 Three triggers
        - Dimensionality: 10 - 12 Four triggers
        - Dimensionality: 13 - 15 Five triggers
        - Dimensionality: 16 - 18 Three independent variables for neighbouring faces of first segment
        - Dimensionality: 19 - 21 Three independent variables for neighbouring faces of first & second segment
        - Dimensionality: 21 - 24 Three independent variables for neighbouring faces of first, second and third segment
        - Dimensionality: 25 - 27 Three independent variables for neighbouring faces of first, second, third and fourth
                                 segment
        - Dimensionality: 28 - 30 Three independent variables for neighbouring faces of all segments
        ...
        """

        ### # NOTE: END OF NEW DEFINITION
        
        self.variable_ranges = self._generate_variable_ranges_map(dimension)
        self.input_file_name = 'combine.k'
        self.output_file_name = 'combineT01.csv'
        self.starter_out_file_name = 'combine_0000.out'
        self.track_node_key = 'DATABASE_HISTORY_NODE99999'
        self.impactor_force_key = "TH-RWALL1IMPACTOR"

        if self.sequential_id_numbering:
            self.problem_id = CrashTube.instance_counter
            CrashTube.instance_counter+=1

    def _write_input_file(self, fem_space_variable_array):
        self.mesh = CrashTubeMesh(fem_space_variable_array,
                                  h_level=self._runner_options('h_level'),
                                  gmsh_verbosity=self._runner_options('gmsh_verbosity')
                                  ) 
        self.model = CrashTubeModel(self.mesh)
        self.model.write_input_files()
    

    @staticmethod
    def _generate_variable_ranges_map(dimension:int)->List[tuple]:
        r"""
        Generates the ranges map given a dimensionality; This is a 
        static method, which is provided with the class to call methods
        outsude the framework.

        Args
        ----------------------
        - dimension: `int`: An integer denoting the dimensionality set by the user

        Returns
        ----------------------
        - `List[tuple]`: A list of tuples indicating the physical ranges of the problem.
        """

        patterns = [(-10, 10), (-4, 4), (0, 4)]

        variable_ranges_map = {}

        for dim in range(1, dimension+1):
            variable_ranges_map[dim] = [patterns[i % 3] for i in range(dim)]
        
        return variable_ranges_map[dimension]
    
    @property
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the Crash Tube problem.
        """
        return ['penalized_sea', "penalized_mass", "absorbed_energy","specific_energy_absorbed"]
     
