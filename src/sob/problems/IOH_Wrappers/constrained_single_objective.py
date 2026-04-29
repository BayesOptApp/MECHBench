import numpy as np
from typing import List, Callable, Optional, Any, Union
from pathlib import Path
from ioh.iohcpp.problem import RealSingleObjective
from ioh.iohcpp import RealBounds, RealSolution, RealConstraint
from ioh.iohcpp import ConstraintEnforcement

#from src.sob.problems.IOH_Wrappers.real_constraint import IOH_Real_Constraint_Wrapper
from src.sob.physical_models import get_model, AbstractPhysicalModel



class IOH_Constrained_Single_Objective_Wrapper(RealSingleObjective):
    """
    A wrapper for IOH single-objective problems that adds constraint handling.
    """
    
    def __init__(
            self,
            model_number: int, 
            dimension: int,
            output_data_labels: List[str],
            functional_definition_objective:Callable,
            functional_definition_constraints: Union[List[Callable],Callable],
            problem_name: str,
            runner_options: dict,
            root_folder: Optional[Union[str,Path]]=None, 
            problem_id: int=10000,  
            instance_id: int=1,
            )->None:
        r"""
        Initialize the constrained single-objective problem.
        
        Args:
            model_number (int): The identifier for the physical model to be used. (1 to 3)
            dimension (int): The dimensionality of the input space for the optimization problem.
            output_data_labels (list[str]): List of output data labels to be considered from the physical model.
            functional_definition_objective (callable): Function to compute the objective value.
            functional_definition_constraints (list[callable] or callable): Functions to compute constraint values.
            problem_name (str): The name of the optimization problem.
            runner_options (dict): Options for the physical model runner.
            root_folder (Union[strPath], optional): The root folder for storing results. Defaults to None.
            instance_id (int, optional): The instance identifier for the problem. Defaults to 1.
        """

        # Assertions to ensure correct input types
        assert functional_definition_objective is not None, "A functional definition for the objective must be provided."
        assert callable(functional_definition_objective), "The functional definition for the objective must be callable."

        self._functional_definition_objective = functional_definition_objective

        # Normalize constraints to a list
        if callable(functional_definition_constraints):
            constraint_list = [functional_definition_constraints]
        else:
            constraint_list = functional_definition_constraints

        # Create the decorator once
        constraint_decorator = make_constraint_decorator(self)

        


        # Set the bounds
        bounds = RealBounds(size=dimension,lb=-5.0, ub=5.0) # Fixed bounds for now, can be modified later
        
        optimum = RealSolution(x = [0.0 for _ in range(dimension)], y = -np.inf) # Arbitrary optimum, can be modified later

        # Set the physical model based on the provided model number and output data labels for the objective
        self._physical_model: AbstractPhysicalModel = get_model(model_type=model_number,
                                                                          dimension=dimension,
                                                                          output_data=output_data_labels,
                                                                          runner_options=runner_options,
                                                                          root_folder=root_folder
                                                                        )

        # Initialize the base single-objective problem
        super().__init__(name = problem_name,
                         n_variables = dimension,
                         instance = instance_id,
                         is_minimization=True,
                         bounds=bounds, # This is fixed for now, can be modified later
                         #constraints=functional_definition_constraints,
                         constraints=[],
                         optimum= optimum)
        

        # Intialize the output values as an empty list
        self._output_values: Any = []

        # Set the problem ID
        super().set_id(problem_id)

        # Set the instance ID
        super().set_instance(instance_id)
        
        # Convert python constraints → IOH constraints
        for ii, func in enumerate(constraint_list):

            wrapped = constraint_decorator(func)

            ioh_constraint = RealConstraint(
                wrapped,
                name=f"constraint_{ii+1}",
                weight=1.0,
                exponent=1.0,
                enforced=ConstraintEnforcement.HIDDEN # Not enforced for now. Hidden constraints do not affect the optimization process.
            )

            # Register constraint with IOH
            self.add_constraint(ioh_constraint)


        
        
    
    def evaluate(self, 
                 x: List[float]) -> float:
        r"""
        Evaluate the objective function at the given point. This method runs the physical model
        and computes the objective value using the provided functional definition.

        The physical model is executed with the input vector `x`, and the output data is
        passed to the functional definition to compute the objective value.

        NOTE: Constraints are handled separately by the IOH framework.
        This method only computes the objective value and follows the standard IOH evaluation protocol.

        Args:
            x (List[float]): The input vector where the objective function is evaluated. (Unused here, but kept for compatibility)
            

        Returns:
            float: The computed objective value.
        """
        
        # Get the output data from the physical model
        output_values = self.output_values


        # Compute the objective value using the functional definition
        objective_value = self._functional_definition_objective(x,
                                                                *output_values,
                                                                **self.generate_kwargs())

        return objective_value

    def __call__(self, x:list) -> float:
        # Run the physical model to get output data
        self._output_values = self._physical_model(x,deck_id=self.state.evaluations+1)

        # Generate kwargs for the functional definition
        kwargs = self.generate_kwargs()

        # Call the super __call__
        val = super().__call__(x)
        return val
    
    def generate_kwargs(self) -> dict:
        r"""
        Generate keyword arguments for the functional definitions based on the output data labels.

        Returns:
            dict: A dictionary mapping output data labels to their corresponding indices.
        """
        return {self.output_data_labels[ii]: ii for ii in range(len(self.output_data_labels))}


    @property
    def output_values(self) -> Any:
        r"""
        Returns the output values from the last evaluation of the physical model.

        Returns:
            Any: The output values from the physical model.
        """
        return self._output_values


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
        return self.physical_model.output_data if isinstance(self.physical_model.output_data, list) else [self.physical_model.output_data]
    
    @output_data_labels.setter
    def output_data_labels(self, labels: Union[List[str],str]):
        r"""
        Sets the list of output data labels considered by the physical model.

        Args:
            labels (List[str]): The new list of output data labels.
        """
        self.physical_model.output_data = labels
    
    @property
    def functional_definition_objective(self) -> Callable:
        r"""
        Returns the functional definition used to compute the objective value.

        Returns:
            Callable: The functional definition function.
        """
        return self._functional_definition_objective
    
    
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

# The following is a function wrapper for each of the constraints 
def make_constraint_decorator(problem: IOH_Constrained_Single_Objective_Wrapper):
    """
    Creates a decorator-like wrapper that transforms a Python constraint function
    into an IOH-compatible constraint callback.
    """
    def constraint_wrapper(func: Callable):
        
        def wrapped_constraint(x: list) -> float:
            # Use latest model output
            output_values = problem.output_values
            kwargs = problem.generate_kwargs()
            return func(x, *output_values, **kwargs)
        
        return wrapped_constraint
    
    return constraint_wrapper