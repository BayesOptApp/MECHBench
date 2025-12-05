from src.sob.physical_models import get_model
from src.sob.physical_models.abstractPhysicalModel import AbstractPhysicalModel
import numpy as np
from src.sob.observer import Observer
from typing import Optional, Union, Dict, Iterable, List
from pathlib import Path


class Sampler:
    """
    A class that wraps a physical model and pairs it with an observer (logger).
    """
    
    def __init__(self,
                 model_number:int,
                 dimension:int,
                 runner_options:Dict,
                 output_data:Optional[Union[Iterable,str]]=["intrusion"], # Default output data
                 root_folder:Optional[Union[str,Path]]=None,
                 ):
        """
        Initialize the Sampler with a physical model and optional observer.
        
        Args:
            model_number: The model type number to instantiate the physical model.
            observer: The observer/logger instance (optional)
        """
        # Get the physical model
        self.model = get_model(model_type=model_number,
                                    dimension=dimension,
                                    runner_options=runner_options,
                                    output_data=output_data,
                                    sequential_id_numbering=True, # Set to False for Sampler usage
                                    root_folder=root_folder,
                                    )
        
        # Initialize observer to None
        self._observer:Optional[Observer] = None
    
    def __call__(self, 
                 vector:Union[List[float],np.ndarray], 
                 *args, 
                 **kwargs)->Union[float, List[float]]:
        """
        Allow the Sampler instance to be called like a function.
        
        Args:
            vector: Input vector for the model
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        
        Returns:
            The result from the model execution with the given input vector
        """

        # Transform vector into a list if it's a numpy array
        if isinstance(vector, np.ndarray):
            vector = vector.tolist()
        
        # Evaluate the model with the provided vector
        result = self.model(vector, )


        # If observer is set, log the result
        if self.observer:
            # Make a list with the input variables (marked from x0, x1, ...)
            input_fields = {f"x{i}": val for i, val in enumerate(vector)}
            # Make a list with the output variables (marked from y0, y1, ...)
            if isinstance(result, (list,tuple)):
                output_fields = {f"{self.model.output_data[i]}": val for i, val in enumerate(result)}
            else:
                output_fields = {f"{self.model.output_data}": result}

            # Combine input and output fields
            log_fields = {**input_fields, **output_fields}
            self.observer.log(**log_fields)
        
        return result
    
    
    def set_observer(self, observer):
        """
        Set or update the observer/logger.
        
        Args:
            observer: The new observer/logger instance
        """
        self.observer = observer
    
    @property
    def model(self)->AbstractPhysicalModel:
        """
        Get the wrapped physical model.
        
        Returns:
            The physical model instance
        """
        return self._model
    
    @model.setter
    def model(self, value:AbstractPhysicalModel)->None:
        """
        Set the wrapped physical model.
        
        Args:
            value: The physical model instance to set
        """
        
        if not isinstance(value, AbstractPhysicalModel):
            raise TypeError("model must be an instance of AbstractPhysicalModel")
        
        self._model:AbstractPhysicalModel = value
    
    @property
    def output_data(self)->Optional[Union[Iterable,str]]:
        """
        Get the output data configuration of the model.
        
        Returns:
            The output data configuration
        """
        return self._model.output_data
    
    @property
    def root_folder(self)->Path:
        """
        Get the root folder of the model.
        
        Returns:
            The root folder path
        """
        return self._model.root_folder
    
    # This is now to handle the observer

    @property
    def observer(self)->Union[Observer,None]:
        """
        Get the observer/logger instance.
        
        Returns:
            The observer/logger instance
        """
        return self._observer
    
    @observer.setter
    def observer(self, value:Observer)->None:
        """
        Set or update the observer/logger instance.
        
        Args:
            value: The observer/logger instance to set
        """
        if value is not None and not isinstance(value, Observer):
            raise TypeError("observer must be an instance of Observer or None")
        
        self._observer:Observer= value

        # Change the model directory if observer is set
        if self._observer is not None:
            self._model.root_folder = self._observer.main_folder_path

    
    @observer.deleter
    def observer(self)->None:
        """
        Remove the observer/logger instance.
        """

        if self._observer is not None:
            # Finalize the observer before removing
            self._observer.finalize(model=self._model.input_file_name)
        self._observer = None

        

    def attach_observer(self, observer:Observer)->None:
        """
        Attach an observer/logger to the sampler.
        
        Args:
            observer: The observer/logger instance to attach
        """
        self.observer = observer
    
    def detach_observer(self)->None:
        """
        Detach the current observer/logger from the sampler.
        """
        del self.observer
