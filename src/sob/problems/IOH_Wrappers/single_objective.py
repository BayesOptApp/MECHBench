from __future__ import annotations

import numpy as np
from ioh.iohcpp.problem import RealSingleObjective
from ioh.iohcpp import RealBounds, RealSolution
from src.sob.physical_models import get_model, AbstractPhysicalModel
from pathlib import Path
from typing import Callable, Union, List, Optional


class IOH_Single_Objective_Wrapper(RealSingleObjective):
    def __init__(self, 
                 model_number: int, 
                 dimension: int,
                 output_data_labels: List[str],
                 functional_definition:Callable,
                 problem_name: str,
                 runner_options: dict,
                 root_folder: Optional[Union[str,Path]]=None,
                 problem_id: int=10001,   
                 instance_id: int=1,   
                 ):
        r"""
        This class wraps the RealSingleObjective class from the IOH library to create a 
        single-objective optimization problem based on a physical model defined.

        Args:
            model_number (int): The identifier for the physical model to be used. (1 to 3)
            dimension (int): The dimensionality of the input space for the optimization problem.
            output_data_labels (list[str]): List of output data labels to be considered from the physical model.
            functional_definition (callable): A function that defines how to compute the objective value from the output data.
            problem_name (str): The name of the optimization problem.
            runner_options (dict): Options for the physical model runner. Defaults to None.
            root_folder (Union[strPath], optional): The root folder for storing results. Defaults to None.
            problem_id (int, optional): The IOH problem identifier. Defaults to 1.
            instance_id (int, optional): The instance identifier for the problem. Defaults to 1 TODO> Instances to be defined properly afterwards.

        """

        assert model_number in [1,2,3], "Model number must be 1, 2, or 3."
        assert functional_definition is not None, "A functional definition must be provided."
        assert callable(functional_definition), "The functional definition must be callable."

        self._functional_definition = functional_definition
        
        # Set the physical model based on the provided model number and output data labels
        self._physical_model: AbstractPhysicalModel = get_model(model_type=model_number,
                                                                dimension=dimension,
                                                                output_data=output_data_labels,
                                                                runner_options=runner_options,
                                                                root_folder=root_folder
                                                                )
        
        # Set the bounds
        bounds = RealBounds(size=dimension,lb=-5.0, ub=5.0) # Fixed bounds for now, can be modified later

        # Define an arbitrary optimum for the problem (can be adjusted later)
        optimum = RealSolution(x = [0.0 for _ in range(dimension)], y = -np.inf)

        super().__init__(name=problem_name,
                         n_variables=dimension,
                         instance=instance_id,
                         is_minimization=True,
                         bounds=bounds, # This is fixed for now, can be modified later
                         constraints=[], # Empty constraints for now
                         optimum= optimum)
        
        # Set the problem id
        super().set_id(problem_id)

        # Set the instance id
        super().set_instance(instance_id)

    
    def evaluate(self, 
                 x: List[float]) -> float:
        r"""
        This wraps the original `evaluate` method from the RealSingleObjective class.
        It evaluates the objective function at a given point `x` by utilizing the associated physical model
        and the provided functional definition.
        Args:
            x (np.ndarray): The input vector where the objective function is evaluated.
        Returns:
            float: The computed objective value at the input vector `x`.
        """
        
        # Get the output data from the physical model
        output_values:Union[List[float],float] = self._physical_model(x,
                                                                    deck_id=self.state.evaluations+1)
        
        kwargs = {self._physical_model.output_data[ii]: ii for ii in range(len(self._physical_model.output_data))} if isinstance(output_values, (list,tuple)) else {self._physical_model.output_data: output_values}
        # Compute the objective value using the functional definition
        objective_value = self._functional_definition(x,
                                                      *output_values, 
                                                      **kwargs)

        return objective_value
    


    @property
    def physical_model(self) -> AbstractPhysicalModel:
        r"""
        Returns the physical model associated with this optimization problem.

        Returns:
            AbstractPhysicalModel: The physical model instance.
        """
        return self._physical_model
    
    @property
    def dimension(self) -> int:
        r"""
        Returns the dimensionality of the input space for the optimization problem.

        Returns:
            int: The number of dimensions.
        """
        return self.meta_data.n_variables
    
    @property
    def root_folder(self) -> Union[str,Path]:
        r"""
        Returns the root folder for storing the input deck files.

        Returns:
            Union[str,Path]: The root folder path.
        """
        return self.physical_model.root_folder
    
    @root_folder.setter
    def root_folder(self, value: Optional[Union[str,Path]]):
        r"""
        Sets the root folder for storing the input deck files.

        Args:
            value (Union[str,Path]): The new root folder path.
        """
        self.physical_model.root_folder = value

    
    @property
    def output_data_labels(self) -> Union[List[str]]:
        r"""
        Returns the list of output data labels considered by the physical model.

        Returns:
            Union[List[str]]: The list of output data labels.
        """
        return self.physical_model.output_data
    
    @output_data_labels.setter
    def output_data_labels(self, labels: Union[List[str],str]):
        r"""
        Sets the list of output data labels considered by the physical model.

        Args:
            labels (List[str]): The new list of output data labels.
        """
        self.physical_model.output_data = labels
    
    @property
    def functional_definition(self) -> Callable:
        r"""
        Returns the functional definition used to compute the objective value.

        Returns:
            Callable: The functional definition function.
        """
        return self._functional_definition
    

    
    # NOTE: This section of the module is for not implemented features yet.
    # They will be implemented in future versions.

    def set_id(self, problem_id: int):
        r"""
        Sets the problem identifier and instance for the optimization problem.

        Args:
            problem_id (int): The new problem identifier.
        """

        raise NotImplementedError("Setting problem ID is not implemented yet.")
    
    def set_instance(self, instance_id: int):
        r"""
        Sets the instance identifier for the optimization problem.

        Args:
            instance_id (int): The new instance identifier.
        """

        raise NotImplementedError("Setting instance ID is not implemented yet.")
    
    def set_name(self, name: str):
        r"""
        Sets the name of the optimization problem.

        Args:
            name (str): The new problem name.
        """

        raise NotImplementedError("Setting problem name is not implemented yet.")


