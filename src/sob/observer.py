from typing import Optional, Union
from pathlib import Path
import csv
import json
import time
import os


class Observer:
    """
    Observer that creates a unique output folder, logs CSV rows, 
    and writes a JSON config at the end.
    """
    def __init__(self,
                 folder_name: str,
                 root: Optional[Union[str, Path]] = None):
        
        self._root = Path(root) if root is not None else Path.cwd().absolute()
        if not self._root.exists():
            self._root.mkdir(parents=True, exist_ok=True)

        self._folder_name = None
        self._main_folder_path = None

        # This will create the directory and ensure uniqueness
        self.folder_name = folder_name

        self.start_time = time.perf_counter()
        self._csv_initialized = False
        self._fieldnames = None
        self._rows_written = 0

        self.csv_path = None
        self.config_path = None

    # ============================================================
    # Logging
    # ============================================================

    def log(self, **fields):
        r"""Log a row of key/value data to CSV."""


        fields = {"time": time.perf_counter() - self.start_time, **fields}

        if not self._csv_initialized:
            self._initialize_csv(fields.keys())

        with open(self.csv_path, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self._fieldnames)
            writer.writerow(fields)

        self._rows_written += 1

    def _initialize_csv(self, keys):
        """Create CSV with header row."""
        self._fieldnames = list(keys)
        data_folder = self._main_folder_path.joinpath("data").absolute()
        if not data_folder.exists():
            data_folder.mkdir(exist_ok=True, parents=True)

        self.csv_path = data_folder.joinpath("log.csv").absolute()

        with open(self.csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self._fieldnames)
            writer.writeheader()

        self._csv_initialized = True

    # ============================================================
    # Finalization
    # ============================================================

    def finalize(self, **kwargs):
        """Write JSON summary with model description."""
        model = kwargs.get("model", None)
        if model is None:
            raise NotImplementedError("Model serialization not implemented yet.")

        self.config_path = self._main_folder_path.joinpath("config.json").absolute()

        full_config = {
            "model": model,
            "rows_written": self._rows_written,
            "duration_seconds": time.perf_counter() - self.start_time,
            "csv_path": str(self.csv_path.resolve()) if self.csv_path else None
        }

        with open(self.config_path, "w") as f:
            json.dump(full_config, f, indent=2)

    def summary(self):
        print(f"CSV rows written: {self._rows_written}")
        print(f"Folder: {self._main_folder_path}")
        if self.config_path:
            print(f"Config JSON path: {self.config_path}")

    # ============================================================
    # Properties
    # ============================================================

    @property
    def root(self) -> Path:
        return self._root

    @property
    def folder_name(self) -> str:
        return self._folder_name

    @folder_name.setter
    def folder_name(self, value: str):
        """
        Sets folder name and ensures uniqueness inside `root`.
        Example:
            output → output_1 → output_2 → ...
        """
        base = value
        candidate = self._root / base
        counter = 1

        while candidate.exists():
            candidate = self._root / f"{base}_{counter}"
            counter += 1

        candidate.mkdir(parents=True, exist_ok=True)

        self._main_folder_path = candidate.absolute()
        self._folder_name = candidate.name

    # Read-only helper properties
    @property
    def rows_written(self) -> int:
        return self._rows_written

    @property
    def csv_initialized(self) -> bool:
        return self._csv_initialized

    @property
    def main_folder_path(self) -> Path:
        return self._main_folder_path

    


