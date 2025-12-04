from typing import List, Optional, Union, Iterable
from src.sob.physical_models.abstractPhysicalModel import AbstractPhysicalModel
from src.sob.physical_models.meshes import ThreePointBendingMesh
from src.sob.physical_models.utils.solver_setup import RunnerOptions
from src.sob.physical_models.fem_settings import ThreePointBendingModel
import os, shutil


class ThreePointBending(AbstractPhysicalModel):
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
            self.deck_id = ThreePointBending.instance_counter
            ThreePointBending.instance_counter+=1

    def generate_input_deck(self, variable_array):
        '''
        Special generate_input_deck function for Three point bending model, transform it into thickness mapping.
        '''
        self._validate_variable_array(variable_array)
        
        thickness_array = [] # to get the thickness in the FEM space

        for i, var in enumerate(variable_array):
            mapped_var = self.linear_mapping_variable(var, self.variable_range)
            thickness_array.append(mapped_var)

        original_dir = os.getcwd()
        dir_name = f'{self.__class__.__name__.lower()}_deck{self.deck_id}'
        working_dir = os.path.join(os.getcwd(), dir_name)
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        os.chdir(working_dir)
        self._write_input_file(thickness_array)
        os.chdir(original_dir)

        if self.sequential_id_numbering:
            self.deck_id = ThreePointBending.instance_counter
            ThreePointBending.instance_counter+=1

    def _copy_files_to_deck(self):
        # Get the path to the lib directory relative to the current file
        deck_name = f'{self.__class__.__name__.lower()}_deck{self.deck_id}'

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
        
        self.fem_model = ThreePointBendingModel(self.mesh)
        self.fem_model.write_input_file(thickness_array)
    
    def mass_calculation(self):
        r"""
        This is a refactoring of the mass calculation function, which is
        used to calculate the mass of the model. The function is called
        when the simulation is run, and it returns the mass of the model.

        Since the mass is computed in metric tons, then this function
        converts the mass to kg. 
        """
        return super().mass_calculation()*1000.0
    
    def peak_force_calculation(self) -> float:
        r"""
        This is a refactoring of the peak force calculation function, which is
        used to calculate the peak force of the model. The function is called
        when the simulation is run, and it returns the peak force of the model.

        Since the force is computed in N, then this function converts the
        force to kN.
        """
        return super().peak_force_calculation()/1000.0
    
    def mean_force_calculation(self) -> float:
        r"""
        This is a refactoring of the mean force calculation function, which is
        used to calculate the mean force of the model. The function is called
        when the simulation is run, and it returns the mean force of the model.

        Since the force is computed in N, then this function converts the
        force to kN.
        """
        return super().mean_force_calculation()/1000.0
    
    @property
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the Three Point Bending problem.
        """
        return ['penalized_sea']