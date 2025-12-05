from abc import ABC, abstractmethod

class AbstractMeshSettings(ABC):
    r"""
    This is an abstract definition for the mesher settings
    """

    @abstractmethod
    def __init__(self, variable_array, h_level,**kwargs)->None:
        """
        Actual class initializer
        """
        pass

    @abstractmethod
    def volume(self)->float:
        r"""Returns the volume (in object units) of the structural object.       
        """
        pass

    @abstractmethod
    def write_mesh_file(self)->None:
        r"""
        Writes the mesh file for the simulator (LS-Dyna/OpenRadioss)
        """
        pass

    @property
    @abstractmethod
    def characteristic_length(self)->None:
        r"""
        Return the characteristic length of the structure. This is a property to set up
        computations of other Crashworthiness Objective metrics. Namely, for future additions,
        the characteristic length is the length of the structure in the direction of the impact.
        """
        pass

    @property
    def h_level(self)->int:
        r"""
        Mesh refinement level
        """
        return self.__h_level
    
    @h_level.setter
    def h_level(self,new_h_level:int)->None:
        r"""
        Mesh refinement setter

        Args
        ----------------------
        new_h_level (`int`): An integer defining the mesh refinement level
        """

        if isinstance(new_h_level,(int,float)) and int(new_h_level) >= 1:
            # Set the h_level in this case
            self.__h_level = int(new_h_level)
        else:
            raise ValueError("The mesh refinement level `h_level` must be a positive integer greater than 1.")
    
    @property
    def gmsh_verbosity(self)->bool:
        r"""
        GMSH verbosity switch
        """
        return self.__gmsh_verbosity
    
    @gmsh_verbosity.setter
    def gmsh_verbosity(self, new_gmsh_verbosity:bool)->None:
        r"""
        GMSH verbosity setter

        Args
        ----------------------
        new_gmsh_verbosity (`bool`): A boolean defining the verbosity of the gmsh mesher
        """

        if isinstance(new_gmsh_verbosity,(bool,int)):
            # Set the gmsh_verbosity in this case
            self.__gmsh_verbosity = bool(new_gmsh_verbosity)
        else:
            raise ValueError("The gmsh verbosity `gmsh_verbosity` must be a boolean.")