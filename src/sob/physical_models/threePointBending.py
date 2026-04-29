from typing import List, Optional, Union, Iterable
from src.sob.physical_models.abstractPhysicalModel import AbstractPhysicalModel, Optional, Path, Union
from src.sob.physical_models.meshes import ThreePointBendingMesh
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
                 runner_options:dict,
                 root_folder:Optional[Union[str,Path]]=None) -> None:
        
        # Get the output data if it is not provided
        if output_data is None:
            output_data = "penalized_mass"
        elif isinstance(output_data, str):
            if output_data in self.forbidden_output_data:
                raise ValueError(f"The output data {output_data} is not allowed for the Three Point Bending problem.")
        elif isinstance(output_data, tuple) or isinstance(output_data, list):
            for elem in output_data:
                if elem in self.forbidden_output_data:
                    raise ValueError(f"The output data {elem} is not allowed for the Three Point Bending problem.")

        
        super().__init__(
                         dimension, 
                         output_data, 
                         runner_options, 
                         root_folder
                         )
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


    def generate_input_deck(self, 
                            variable_array, 
                            deck_id: int):
        '''
        Special generate_input_deck function for Three point bending model, transform it into thickness mapping.
        '''
        self._validate_variable_array(variable_array)
        
        thickness_array = [] # to get the thickness in the FEM space

        for i, var in enumerate(variable_array):
            mapped_var = self.linear_mapping_variable(var, self.variable_range)
            thickness_array.append(mapped_var)

        original_dir = self.root_folder.absolute()
        dir_name = f'{self.__class__.__name__.lower()}_deck{deck_id}'
        working_dir = original_dir.joinpath( dir_name)
        if not working_dir.exists():
            working_dir.mkdir(parents=True, exist_ok=True)

        os.chdir(working_dir.absolute().as_posix())
        self._write_input_file(thickness_array, deck_id=deck_id)
        os.chdir(original_dir.absolute().as_posix())


    def _copy_files_to_deck(self,
                            deck_id:int):
        
        # Get the path to the lib directory relative to the current file
        deck_name = f'{self.__class__.__name__.lower()}_deck{deck_id}'

        lib_dir = Path(os.path.join(os.path.dirname(__file__), 'lib'))
        
        # Source file paths
        source_files = ['ThreePointBending_0001.rad']  # Copy base starter file and engine files
        # Destination folder path
        dest_folder = Path(os.getcwd())
        
        # Create the destination folder if it doesn't exist
        if not dest_folder.exists():
            dest_folder.mkdir(parents=True, exist_ok=True)
        
        # Copy each source file to the destination folder
        for file_name in source_files:
            source_path = lib_dir.joinpath(file_name)
            dest_path = dest_folder.joinpath(file_name)
            shutil.copyfile(source_path.as_posix(), dest_path.as_posix())

    def _write_input_file(self, 
                          thickness_array,
                          deck_id: int):
        self._copy_files_to_deck(deck_id=deck_id)
        self.mesh = ThreePointBendingMesh(thickness_array,
                                          h_level=self._runner_options.h_level,
                                          gmsh_verbosity=self._runner_options.gmsh_verbosity
                                          )
        
        self.fem_model = ThreePointBendingModel(self.mesh)
        self.fem_model.write_input_file(thickness_array)
    
    def mass_calculation(self, deck_id:int) -> float:
        r"""
        This is a refactoring of the mass calculation function, which is
        used to calculate the mass of the model. The function is called
        when the simulation is run, and it returns the mass of the model.

        Since the mass is computed in metric tons, then this function
        converts the mass to kg. 
        """
        return super().mass_calculation(deck_id=deck_id)*1000.0
    
    def peak_force_calculation(self, deck_id:int) -> float:
        r"""
        This is a refactoring of the peak force calculation function, which isroot_folder:Optional[Union[str,Path]]
        used to calculate the peak force of the model. The function is called
        when the simulation is run, and it returns the peak force of the model.

        Since the force is computed in N, then this function converts the
        force to kN.
        """
        return super().peak_force_calculation(deck_id=deck_id)/1000.0
    
    def mean_force_calculation(self, deck_id:int) -> float:
        r"""
        This is a refactoring of the mean force calculation function, which is
        used to calculate the mean force of the model. The function is called
        when the simulation is run, and it returns the mean force of the model.

        Since the force is computed in N, then this function converts the
        force to kN.
        """
        return super().mean_force_calculation(deck_id=deck_id)/1000.0
    
    @property
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the Three Point Bending problem.
        """
        return ['penalized_sea']