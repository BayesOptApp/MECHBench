r"""
This is a module to generalize a template procedure to 
generate meshes with GMSH
"""

import numpy as np
from typing import Tuple, List, Union, Optional, Collection
import json
from pathlib import Path
import os

from abc import ABC, abstractmethod
import sys

try: 
    import gmsh
except ModuleNotFoundError:
    print("GMSH is not installed in the current Python environment!")


class Template_GMSH_Mesh_Constructor(ABC):
    r"""
    This is a template class to define the GMSH mesh class constructor
    for all the remaining cases
    """

    @abstractmethod
    def __init__(self, model_name:str="Default")->None:
        r"""
        Main class initializer;
        This initializer just sets a model name
        a model_name given by parameter
        """

        self.model_name = model_name
    
    def _end_gmsh(self)->None:
        try:
            gmsh.finalize()
        except:
            print("GMSH job cannot be closed!")
        
    @abstractmethod
    def __call__(self, *args, **kwds):
        pass
    
    @staticmethod
    def load_json_file(path_to_file: Union[str, Path]) -> dict:
        r"""
        This is just a member function to load a JSON file
        into memory to define a mesh
        """
        if isinstance(path_to_file, str):
            path_to_file_path = Path(path_to_file)
        elif isinstance(path_to_file, Path):
            path_to_file_path = path_to_file
        else:
            raise ValueError("The path given as a parameter is not a string or Path object")

        with open(path_to_file_path, 'r') as json_file:
            dictt = json.load(json_file)

        return dictt
    
    @property
    def model_name(self)->str:
        return self._model_name
    
    @property
    @abstractmethod
    def required_parameters(self)->Union[tuple,list]:
        r"""
        This is a property the derived classes must have which depends 
        on defining a set of required parameters to construct the mesh
        """
        pass

    def check_parameters(self,params_dict:dict)->bool:
        r"""
        This function just checks a parameter set from a dictionary and ensures 
        all the required parameters of the class are contained within the dictionary.
        """
        count = 0
        total_params = len(self.required_parameters)
        missing_params = []

        for key in params_dict.keys():
            if key.lower() in self.required_parameters:
                count+=1
            else:
                missing_params.append(key)

        
        if count == total_params:
            return True
        else:
            print("Missing parameters",missing_params)
            return False

        
    
    @model_name.setter
    def model_name(self,new_model_name:object)->None:
        try:
            model_name_str:str = str(new_model_name)
        except Exception as e:
            print(e.args)
        
        self._model_name = model_name_str
        


        

