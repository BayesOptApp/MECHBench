import os
from abc import ABC, abstractmethod
from src.sob.physical_models.meshes import AbstractMeshSettings



class AbstractFEMSettings(ABC):
    r"""
    This is a class used to define the methods and properties of
    the Setup Cards for each of the methods
    """

    @abstractmethod
    def __init__(self, mesh, **kwargs)->None:
        pass
    
    # @abstractmethod
    # def write_mesh_file(self):
    #     r"""
    #     A method to activate writing the Mesh part of the Simulation Card
    #     for either OpenRadioss or LS-Dyna defined cards
    #     """
    #     pass
    
    # @abstractmethod
    # def _load_impactor(self):
    #     r"""
    #     Loads the impactor properties
    #     """
    #     pass
    
    @abstractmethod
    def absorbed_energy(self):
        # initial kinetic energy
        pass

    
    @property
    def mesh(self)->AbstractMeshSettings:
        r"""
        Returns the mesh object linked to the model
        """
        return self._mesh
    
    @mesh.setter
    def mesh(self, new_mesh:AbstractMeshSettings)->None:
        r"""
        Sets a new mesh object by parameter
        """
        if not isinstance(new_mesh,AbstractMeshSettings):
            raise TypeError("The mesh must be an instance of AbstractMeshSettings or its subclasses.")
        self._mesh = new_mesh
    
    @property
    @abstractmethod
    def material_card_type(self)->str:
        r"""
        Determines the material card type (OpenRadioss, LS-Dyna)
        """
        pass

    @property
    @abstractmethod
    def impactor_offset(self)->float:
        r"""
        Defines the offset of the impactor at initial time from 
        impacting the structure
        """
        pass
