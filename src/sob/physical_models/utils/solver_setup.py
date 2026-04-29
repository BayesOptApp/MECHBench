### THIS IS A MODULE TO SET-UP THE CASE FROM THE START

# Module Properties
__author__ = "Ivan Olarte Rodriguez"


# Module Imports
from sys import platform
import os
import warnings
from collections import namedtuple, OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Union



from dataclasses import dataclass, fields
from typing import Any, Dict

import src.sob.physical_models.utils.platform_det as platform_det

@dataclass
class RunnerOptions:
    open_radioss_main_path: Optional[Path] = None

    write_vtk: int = 0
    h_level: int = 1
    nt: int = 1
    np: int = 1
    gmsh_verbosity: int = 0
    save_mesh_vtk: int = 0

    def __post_init__(self):
        self._validate()
        self._setup_openradioss_path()

    # ---------------------------
    # Validation
    # ---------------------------
    def _validate(self):
        for name in ("h_level", "nt", "np"):
            val = getattr(self, name)
            if not isinstance(val, int) or val < 1:
                raise ValueError(f"{name} must be >= 1")

        for name in ("write_vtk", "save_mesh_vtk"):
            if getattr(self, name) not in (0, 1):
                raise ValueError(f"{name} must be 0 or 1")

        if self.gmsh_verbosity not in (0, 1):
            raise ValueError("gmsh_verbosity must be 0 or 1")

    # ---------------------------
    # Path setup
    # ---------------------------
    def _setup_openradioss_path(self):
        system = platform_det.platform_detection()

        if self.open_radioss_main_path is None:
            platform_det.raise_if_not_allowed_platform()

            folder_map = {
                "Linux": "OpenRadioss_linux64",
                "Windows": "OpenRadioss_win64",
            }

            base_path = Path.cwd() / folder_map[system] / "OpenRadioss"

            if not base_path.exists():
                print(f"OpenRadioss not found at {base_path}. Downloading...")
                platform_det.download_zip_openradioss()

            self.open_radioss_main_path = base_path.resolve()

        else:
            self.open_radioss_main_path = Path(
                self.open_radioss_main_path
            ).resolve()

            if not self.open_radioss_main_path.exists():
                raise ValueError(
                    f"Path does not exist: {self.open_radioss_main_path}"
                )

    # ---------------------------
    # Helpers
    # ---------------------------
    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RunnerOptions":
        valid = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in valid})

    def as_dict(self) -> Dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    def normalize(self):
        for name in (
            "h_level", "nt", "np",
            "write_vtk", "save_mesh_vtk",
            "gmsh_verbosity"
        ):
            setattr(self, name, int(getattr(self, name)))