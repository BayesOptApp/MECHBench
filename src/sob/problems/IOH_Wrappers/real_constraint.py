import numpy as np
from typing import List, Callable, Optional, Any, Union
from ioh.iohcpp import RealConstraint
#ioh.iohcpp.ConstraintEnforcement



class IOH_Real_Constraint_Wrapper(RealConstraint):
    """
    A wrapper for IOH real-valued constraints.
    """
    
    def __init__(
            self,
            constraint_function: Callable,
            weight: float = 1.0,
            exponent: float = 1.0,
            name: str = "Real_Constraint",
            )->None:
        r"""
        Initialize the real-valued constraint.

        Args:
            constraint_function (callable): Function to compute the constraint value.
            weight (float, optional): Weight of the constraint. Defaults to 1.0.
            exponent (float, optional): Exponent for the constraint violation. Defaults to 1.
            name (str, optional): Name of the constraint. Defaults to "Real_Constraint".

        """

        assert constraint_function is not None, "A constraint function must be provided."
        assert callable(constraint_function), "The constraint function must be callable."

        super().__init__(constraint_function, # the function to compute the constraint value (check the name after debugging)
                         weight=weight,
                         exponent=exponent,
                         name=name,
                         enforced=0, # NOTE: Set to 0 for now, to set it to be not enforced. Can be modified later.
                            )
        
        #NOTE: Attempt to rewrite the `compute_violation` method if needed.

    def compute_violation(self, 
                          x: List[float], 
                          *args,
                          **kwargs) -> float:
        r"""
        Compute the constraint violation at the given point.

        Args:
            x (List[float]): The input vector where the constraint is evaluated.
            *args: Additional positional arguments for the constraint function.
            **kwargs: Additional keyword arguments for the constraint function.

        Returns:
            float: The computed constraint violation.
        """

        # Compute the constraint value using the provided function
        constraint_value = self.fn(x, *args, **kwargs)

        # Compute the violation based on the constraint value
        violation = max(0.0, constraint_value)

        return violation
