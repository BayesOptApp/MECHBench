from typing import List, Optional, Union, Iterable
from src.sob.physical_models.abstractPhysicalModel import AbstractPhysicalModel
from src.sob.physical_models.meshes import StarBoxMesh
from src.sob.physical_models.fem_settings import StarBoxModel


class StarBox(AbstractPhysicalModel):
    instance_counter = 1

    def __init__(self, 
                 dimension, 
                 runner_options:dict,
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
            self.deck_id = StarBox.instance_counter
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
                                h_level=self._runner_options.h_level,
                                gmsh_verbosity=self._runner_options.gmsh_verbosity
        )
        self.fem_model = StarBoxModel(self.mesh)
        self.fem_model.write_input_files()
    
    @property
    def forbidden_output_data(self)->List[str]:
        """
        Returns a list of forbidden output data for the StarBox problem.
        """
        return ['penalized_mass']
    
    
    