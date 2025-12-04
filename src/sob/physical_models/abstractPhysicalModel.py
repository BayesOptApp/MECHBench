from abc import ABC, abstractmethod
import os
from typing import Iterable, List, Union
import numpy as np

from src.sob.physical_models.fem_settings.abstractFEMSettings import AbstractFEMSettings
from src.sob.physical_models.utils.solver_setup import RunnerOptions
from src.sob.physical_models.solvers.openRadioss_runner import run_OpenRadioss
import pandas as pd

_PRINCIPAL_OUTPUTS = [
    'mass',
    'absorbed_energy',
    'intrusion',
    'mean_impact_force',     
    'max_impact_force'
]

_COMPOSITE_OUTPUTS = [
    'specific_energy_absorbed',
    'load_uniformity',
    'penalized_sea',
    'penalized_mass'
]

_PROBLEM_SPECIFIC_OUTPUTS:tuple = ("penalized_sea","penalized_mass")


class AbstractPhysicalModel(ABC):
    r'''
    The abstract of the structural optimization models, with three subclass (three problem types)
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
        self._deck_id:int = -1 
        self.variable_ranges = None # constraints of the problem
        self.input_file_name = None # input deck name
        self.output_file_name = None # output result name
        self.starter_out_file_name = None

        # The attributes will be loaded if the function _write_input_file has been called. 
        # Used for mass calculation.
        self.mesh = None
        self._fem_model = None

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
    
    def linear_mapping_variable(self, search_space_variable, problem_space_range:tuple):
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
            mapped_var = self.linear_mapping_variable(var, self.variable_ranges[i])
            fem_space_variable_array.append(mapped_var)

        original_dir = os.getcwd()
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.deck_id}'
        print('######################################################\n')
        print(dir_name)
        working_dir = os.path.join(os.getcwd(), dir_name)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        os.chdir(working_dir)
        self._write_input_file(fem_space_variable_array)
        os.chdir(original_dir)

        if self.__sequential_id_numbering:
            self.deck_id += 1 # update the problem id for the input deck generation
    
    @abstractmethod
    def _write_input_file(self, fem_space_variable_array):
        pass
    
    def run_simulation(self,runStarter=False):
        if self.input_file_name is None:
            raise ValueError("input_file_name must be provided or defined in the subclass.")
        # make problem id back to original, since it has been updated when generate_input_deck has been called
        if self.__sequential_id_numbering:
            self.deck_id -= 1

        dir_name = f'{self.__class__.__name__.lower()}_deck{self.deck_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        input_file_path = os.path.join(working_dir, self.input_file_name)

        if runStarter:
            # This is just one bypass in order to avoid setting MP settings
            # for just the mass computation
            run_OpenRadioss(input_file_path, 
                        self.batch_file_path, 
                        runStarter=runStarter,
                        write_vtk=False,
                        np_int=1,
                        nt_int=1)
        else:
            run_OpenRadioss(input_file_path, 
                        self.batch_file_path, 
                        runStarter=runStarter,
                        write_vtk=bool(self._runner_options('write_vtk')),
                        np_int=self._runner_options('np'),
                        nt_int=self._runner_options('nt'))

    def load_output_data_frame(self):
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.deck_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        # load simulation result dataframe and make it cleaner
        output_file_path = os.path.join(working_dir, self.output_file_name)
        self.output_data_frame = pd.read_csv(output_file_path)
        self.output_data_frame.columns = self.output_data_frame.columns.str.replace(' ', '')
        
        # update problem id again
        if self.__sequential_id_numbering:
            self.deck_id += 1

    def extract_mass_from_file(self):
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.deck_id}'
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
 

        return abs(self.output_data_frame[col].abs().max()) - self.fem_model.impactor_offset
    
    def mass_calculation(self)->float:
        if self.sim_status < 1: 
            self.run_simulation(runStarter=True)
            # Change the status
            self.sim_status = 1

        val:float = self.extract_mass_from_file() - self.fem_model.rigid_mass

        if self.__sequential_id_numbering:
            self.deck_id +=1  

        return val
    
    def absorbed_energy_calculation(self)->float:
        return self.fem_model.absorbed_energy()
    
    def _get_max_intrusion_index(self):
        col = [col for col in self.output_data_frame.columns if self.track_node_key in col][2]
        return self.output_data_frame[col].abs().idxmax()

    def _get_force_data(self):
        force_col = [col for col in self.output_data_frame.columns if self.impactor_force_key in col][2]
        max_idx = self._get_max_intrusion_index()
        df = self.output_data_frame.loc[:max_idx]

        time_vect = df["time"].values
        impulse_vec = df[force_col].values

        # Check if the impulse curve is monotonic
        if np.all(np.diff(impulse_vec) >= 0) or np.all(np.diff(impulse_vec) <= 0):
            # Monotonic curve which represents an impulse
            return np.diff(impulse_vec) / np.diff(time_vect)
        else:
            # If the impulse curve is not monotonic, then it corresponds to a force curve and not
            # an impulse curve and we return it as is.
            return impulse_vec

    def peak_force_calculation(self) -> float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.sim_status = 2
            self.load_output_data_frame()
        force_data = self._get_force_data()
        
        return np.abs(np.max(force_data)).astype(float).ravel()

    def mean_force_calculation(self) -> float:
        if self.sim_status < 2 or self.output_data_frame is None: 
            self.run_simulation()
            self.load_output_data_frame()
            self.sim_status = 2


        force_data = self._get_force_data()
        
        return np.abs(np.mean(force_data)).astype(float).ravel()
    


    def __call__(self, variable_array: list, deck_id: int = -1) -> Union[float, List[float]]:
        
        if not self.__sequential_id_numbering:
            self.deck_id = deck_id
        
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
    def deck_id(self)->int:
        return self._deck_id
    
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
    

    @deck_id.setter
    def deck_id(self,new_deck_id:int)->None:
        
        if not isinstance(new_deck_id,int):
            raise TypeError("The type of the object is not correct")
        
        else:
            if new_deck_id <= 0:
                raise ValueError("The value must be greater than 0")
            # Assign the new problem ID
            self._deck_id = new_deck_id
    

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
        elif isinstance(new_output_data,(list,tuple)):
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
    def fem_model(self)->AbstractFEMSettings:
        return self._fem_model
    
    @fem_model.setter
    def fem_model(self, new_model:AbstractFEMSettings)->None:
        if not issubclass(type(new_model),AbstractFEMSettings) and (new_model is not None):
            raise TypeError("The type of the object is not correct")
        else:
            self._fem_model = new_model

    
    @output_data.deleter
    def output_data(self)->None:
        del self._output_data
    
    @property
    @abstractmethod
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the problem.
        """
        pass