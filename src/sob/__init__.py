from src.sob.physical_models import AbstractPhysicalModel,StarBox, ThreePointBending, CrashTube
from typing import Optional, Iterable, Union

def get_problem(model_type:int, 
                dimension:int, 
                runner_options:dict,
                output_data:Optional[Union[Iterable,str]]=None,
                sequential_id_numbering:bool=True,**kwargs)->AbstractPhysicalModel:
    '''
    Generates a problem instance based on the specified model type and configuration.

    Parameters:
    ----------
    model_type : int
        Specifies the type of model to be used. The options are:
        - 1: Crash star box model - Design variables determine the shape and thickness of the star box (up to 5 variables).
        - 2: Three point bending model - Design variables determine the thickness of the reinforcement inner shell (up to 5 variables).
        - 3: Crash tube model - Design variables determine the positions and depths of the three triggers (up to 6 variables).

    dimension : int
        Defines the specific configuration of the model. The allowable ranges and configurations are as follows:
        - For model_type = 1 (Crash star box model):
            - 1: Square
            - 2: Rectangular
            - 3: Rectangular with varying thickness
            - 4: Star shape
            - 5: Star shape with varying thickness
            - 6-34: Star shape with different thickness profiles.
        
        - For model_type = 2 (Three point bending model):
            - 1: All 5 shell thicknesses vary with the same value.
            - 2: Only the first and last shell thicknesses vary; others are fixed at a middle value.
            - 3: The first, middle, and last shell thicknesses vary.
            - 4: All shell thicknesses except the middle vary.
            - 5: All shell thicknesses vary.
            - 6-40: Varying configurations of shell thicknesses with thickness profiles for each rib.

        - For model_type = 3 (Crash tube model):
            - 2: Three positions and depths vary together with the same value.
            - 3: All three trigger positions are fixed; the three trigger depths vary.
            - 4: All vary except for the middle trigger position and depth.
            - 5: All vary except for the middle trigger depth.
            - 6: All three positions and depths vary.

    output_data : str
        Specifies the type of output data required. The options are:
        - 'mass'
        - 'absorbed_energy'
        - 'intrusion' (Requires running FEM simulation while 'mass' and 'absorbed_energy' do not require one.)
        - 'specific-energy-absorbed' (Computes the SEA by operating the mass and absorbed energy)
        - 'mean_impact_force' (Requires running FEM simulation)
        - 'peak-impact-force' (Requires running FEM simulation)
        - 'usage_ratio' (Requires running FEM simulation)
        - 'load_uniformity' (Requires running FEM simulation)
        

    runner_options : dict
        A dictionary containing options for the simulation runner, such as paths and computational settings.
        At least the key "open_radioss_main_path" must be provided.
        
    sequential_id_numbering : bool
        If True, assigns sequential IDs to each problem instance for unique identification.
    '''
    if model_type==1:
        problem_instance = StarBox(dimension=dimension, 
                                   output_data=output_data, 
                                   runner_options=runner_options,
                                   sequential_id_numbering=sequential_id_numbering,
                                   **kwargs)
        return problem_instance
    elif model_type==2:
        problem_instance = ThreePointBending(dimension=dimension, 
                                   output_data=output_data, 
                                   runner_options=runner_options,
                                   sequential_id_numbering=sequential_id_numbering,
                                   **kwargs)
        return problem_instance
    elif model_type==3:
        problem_instance = CrashTube(dimension=dimension, 
                                   output_data=output_data, 
                                   runner_options=runner_options,
                                   sequential_id_numbering=sequential_id_numbering,
                                   **kwargs)
        return problem_instance
    else:
        raise ValueError('Invalid model type')


