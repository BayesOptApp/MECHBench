### THIS IS A MODULE TO SET-UP THE CASE FROM THE START

# Module Properties
__author__ = "Ivan Olarte Rodriguez"


# Module Imports
import warnings
from collections import namedtuple, OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Union



from dataclasses import dataclass, fields
from typing import Any, Dict


@dataclass
class RunnerOptions:
    """
    Clean and safe configuration container for running OpenRadioss cases.
    Replaces the old dynamic RunnerOptions(dictionary) class.

    - No fuzzy matching
    - No eval()
    - No attribute merging
    - No recursion
    - Safe defaults
    """

    # ---- Core configuration fields ----
    open_radioss_main_path: Path = Path("/home/ivanolar/Documents/OpenRadioss2/" + \
                                        "OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui")
    

    write_vtk: int = 0
    h_level: int = 1
    nt: int = 1
    np: int = 1
    gmsh_verbosity: int = 0
    save_mesh_vtk: int = 0

    # -----------------------------------
    # Post-initialization validation
    # -----------------------------------
    def __post_init__(self):

        if isinstance(self.open_radioss_main_path, str):
            self.open_radioss_main_path = Path(self.open_radioss_main_path)
            
        # integer fields that must be >= 1
        for name in ("h_level", "nt", "np"):
            val = getattr(self, name)
            if not isinstance(val, int) or val < 1:
                raise ValueError(f"{name} must be a positive integer (>= 1).")

        # binary fields: 0 or 1
        for name in ("write_vtk", "save_mesh_vtk"):
            val = getattr(self, name)
            if val not in (0, 1):
                raise ValueError(f"{name} must be either 0 or 1.")

        # gmsh verbosity allowed values
        if self.gmsh_verbosity not in (0, 1):
            raise ValueError("gmsh_verbosity must be 0 or 1.")
        
        # Check that the open_radioss_main_path exists
        if not self.open_radioss_main_path.exists():
            raise ValueError(f"open_radioss_main_path '{self.open_radioss_main_path}' does not exist.")

    # -----------------------------------
    # Construction helpers
    # -----------------------------------
    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RunnerOptions":
        """
        Safe constructor to initialize from a dictionary.
        Unknown keys are ignored.
        """
        valid = {f.name for f in fields(cls)}
        filtered = {k: v for k, v in d.items() if k in valid}
        return cls(**filtered)

    def as_dict(self) -> Dict[str, Any]:
        """Return a standard dictionary representation of this configuration."""
        return {f.name: getattr(self, f.name) for f in fields(self)}

    # -----------------------------------
    # Utility
    # -----------------------------------
    def normalize(self):
        """
        Forces integer fields to be integers.
        Useful if values come from JSON or command-line parsing.
        """
        self.h_level = int(self.h_level)
        self.nt = int(self.nt)
        self.np = int(self.np)
        self.write_vtk = int(self.write_vtk)
        self.save_mesh_vtk = int(self.save_mesh_vtk)
        self.gmsh_verbosity = int(self.gmsh_verbosity)
