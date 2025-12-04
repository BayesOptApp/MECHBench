from typing import List, Optional, Union, Iterable
from src.sob.physical_models.abstractPhysicalModel import AbstractPhysicalModel, Optional, Path, Union
from src.sob.physical_models.meshes import CrashTubeMesh
from src.sob.physical_models.fem_settings import CrashTubeModel

class CrashTube(AbstractPhysicalModel):
    r"""
    The Crash Tube optimization problem involves optimizing the positions, heights,
    and depths of triggers on a crash tube structure to improve its performance during impact tests.
    The design variables include the vertical positions, heights, and depths of the triggers,
    which can be adjusted within specified ranges
    """

    instance_counter = 1
    def __init__(self, dimension, output_data, 
                 runner_options:dict,
                 sequential_id_numbering:bool,
                 root_folder:Optional[Union[str,Path]]=None) -> None:
        
        # Get the output data if it is not provided
        if output_data is None:
            output_data = "load_uniformity"
        elif isinstance(output_data, str):
            if output_data in self.forbidden_output_data:
                raise ValueError(f"The output data {output_data} is not allowed for the Crash Tube problem.")
        elif isinstance(output_data, tuple) or isinstance(output_data, list):
            for elem in output_data:
                if elem in self.forbidden_output_data:
                    raise ValueError(f"The output data {elem} is not allowed for the Crash Tube problem.")
 
        
        super().__init__(dimension=dimension, 
                         output_data=output_data, 
                         runner_options=runner_options ,
                         sequential_id_numbering=sequential_id_numbering,
                         root_folder=root_folder)

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
            self.deck_id = CrashTube.instance_counter
            CrashTube.instance_counter+=1

    def _write_input_file(self, fem_space_variable_array):
        self.mesh = CrashTubeMesh(fem_space_variable_array,
                                  h_level=self._runner_options.h_level,
                                  gmsh_verbosity=self._runner_options.gmsh_verbosity
                                  ) 
        self.fem_model = CrashTubeModel(self.mesh)
        self.fem_model.write_input_files()
    

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

        patterns = [(-4, 4), (-10, 10), (0, 4)]

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