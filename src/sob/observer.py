from typing import Optional, Dict, Any, Union
from pathlib import Path
import csv
import json
import time
import os


class Observer:
    def __init__(self,
                 folder_name:str, 
                 root:Optional[Union[str,Path]]=None,

    ):
        """
        Parameters
        ----------
        csv_path : str
            Path to CSV file to write data rows.
        config_path : str
            Path to JSON file to write config metadata at the end.
        config : dict or None
            Metadata/configuration to store in final JSON.
        """
        self._root = Path(root) if root is not None else Path.cwd()
        self._folder_name = folder_name


        self.start_time = time.perf_counter()

        self._csv_initialized = False
        self._fieldnames = None
        self._rows_written = 0

    def log(self, **fields):
        """
        Log a set of named values to the CSV file.

        Example:
            logger.log(step=1, x=0.2, y1=0.04, y2=0.008)
        """
        # Add a timestamp if desired
        fields = {"time": time.time() - self.start_time, **fields}

        if not self._csv_initialized:
            self._initialize_csv(fields.keys())

        # Write row
        with open(self.csv_path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self._fieldnames)
            writer.writerow(fields)

        self._rows_written += 1

    def _initialize_csv(self, keys):
        """Create CSV with header row."""
        self._fieldnames = list(keys)

        # Initialize CSV file
        self.csv_path = self._main_folder_path.joinpath("data/log.csv")

        with open(self.csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self._fieldnames)
            writer.writeheader()

        self._csv_initialized = True

    # -----------------------------
    # Finalize and write JSON config
    # -----------------------------

    def finalize(self, **kwargs):
        """
        Write configuration metadata to the JSON config file.
        """

        model = kwargs.get("model", None)
        if model is None:
            raise NotImplementedError("Model serialization not implemented yet.")
        
        self.config_path = self._main_folder_path.joinpath("config.json")

        full_config = {
            "model": model,
            "rows_written": self._rows_written,
            "duration_seconds": time.time() - self.start_time,
            "csv_path": os.path.abspath(self.csv_path)
        }

        with open(self.config_path, "w") as f:
            json.dump(full_config, f, indent=2)

    # Optional convenience
    def summary(self):
        print(f"CSV rows written: {self._rows_written}")
        print(f"Config JSON path: {self.config_path}")

    @property
    def root(self) -> Path:
        return self._root
    
    @root.setter
    def root(self, value:Union[str,Path]):
        self._root = Path(value)
        self._root.mkdir(parents=True, exist_ok=True)

        self._main_folder_path = self._root.joinpath(self.folder_name)
        
        # If folder exists, find next available number
        counter = 1
        while self._main_folder_path.exists():
            self.folder_name = f"{self.folder_name}_{counter}"
            self._main_folder_path = self._root.joinpath(f"{self.folder_name}_{counter}")
            counter += 1
        
        # Create the folder
        self.folder_name = self._main_folder_path.name
    
    @property
    def folder_name(self) -> str:
        return self._folder_name
    
    @folder_name.setter
    def folder_name(self, value:str):
        self._folder_name = value

        self._main_folder_path = self._root.joinpath(self._folder_name)

        # If folder exists, find next available number
        counter = 1
        while self._main_folder_path.exists():
            self._folder_name = f"{value}_{counter}"
            self._main_folder_path = self._root.joinpath(f"{self._folder_name}_{counter}")
            counter += 1
        
        # Create the folder
        self._main_folder_path.mkdir(parents=True, exist_ok=True)
        self._folder_name = self._main_folder_path.name
    
    @property
    def rows_written(self) -> int:
        return self._rows_written
    
    @property
    def csv_initialized(self) -> bool:
        return self._csv_initialized
    


