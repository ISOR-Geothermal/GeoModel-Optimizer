import layermesh.mesh as lm
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pywaiwera
import json
import shutil
from typing import Tuple, Optional, List, Callable, Any
import os, sys
from pathlib import Path
from copy import deepcopy
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import time
from pandas import Timestamp


class WaiweraRun:
    def __init__(self, params, json_file, temperature, pressure, timestamp, meta: dict | None = None):
        self.params = deepcopy(params)
        self.json_file = json_file
        self.temperature = temperature.copy()
        self.pressure = pressure.copy()
        self.timestamp = timestamp
        
        # keep some meta information about the specifics of this run
        self.meta = {}
        if meta is not None:
            self.meta.update(meta)

        mesh_path = Path(params["mesh"]["filename"])
        if not mesh_path.exists():
            raise FileNotFoundError(f"Could not find the mesh file at {mesh_path}. \
                                    Ensure that the file exists at the location specified in the passed parameters")
        self.mesh = self._load_mesh(mesh_path)
        self.h5_path = Path(params["output"]["filename"])
        if not self.h5_path.exists():
            raise Warning(f"could not find corresponding h5 file {self.h5_path}. This could just mean the run has not \
                          happened yet, or it could mean the path specified is wrong")

    @classmethod
    def run_from_file(cls, json_file: str, num_processes=2, meta: Optional[dict] = None):
        timestamp = Timestamp.now()

        env = pywaiwera.docker.DockerEnv()
        env.run_waiwera(json_file, noupdate=True, num_processes=num_processes)
        
        with open(json_file, "r") as f:
            params = json.load(f)

        h5_name = params["output"]["filename"]
        t = cls.get_temperature(h5_name, -1)
        p = cls.get_pressure(h5_name, -1)

        run = WaiweraRun(params, json_file, t, p, timestamp, meta=meta)
        return run

    @classmethod
    def run_from_dict(cls, params: dict, json_file: str, num_processes=2, meta: Optional[dict] = None):
        timestamp = Timestamp.now()
        with open(json_file, "w") as f:
            json.dump(params, f, indent=4)
            
        with nostdout():
            env = pywaiwera.docker.DockerEnv()
            env.run_waiwera(json_file, noupdate=True, num_processes=num_processes)
        
        h5_name = params["output"]["filename"]
        t = cls.get_temperature(h5_name, -1)
        p = cls.get_pressure(h5_name, -1)

        run = WaiweraRun(params, json_file, t, p, timestamp, meta=meta)
        return run

    @classmethod
    def from_file(cls, h5_filename: str, json_filename: str, meta: Optional[dict] = None):
        with open(json_filename) as f:
            params = json.load(f)
        t = cls.get_temperature(h5_filename, timestep=-1)
        p = cls.get_pressure(h5_filename, timestep=-1)
        timestamp = Timestamp.now()

        run = WaiweraRun(params, json_filename, t, p, timestamp, meta=meta)
        return run

    def plots(self, value: str = "temperature", depth=0):
        if value in ["temperature", "t"]:
            data = self.temperature
        elif value in ["pressure", "p"]:
            data = self.pressure
        else:
            raise ValueError
            
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 4))
        self.mesh.slice_plot("x", axes=ax1, value=data)
        self.mesh.layer_plot(depth, value=data, axes=ax2)
    
    def execute_function(self, function: Callable) -> Any:
        """
        Execute a callable `function` that takes an open h5py.File as input
        """
        with h5py.File(str(self.h5_path)) as data:
            return function(data)
    
    def _load_mesh(self, mesh_path: Path) -> 'layermesh.Mesh':
        if mesh_path.suffix == ".exo":
            mesh_filename = str(mesh_path).replace(".exo", ".h5")
        else:
            mesh_filename = str(mesh_path)
        mesh = lm.mesh(str(mesh_filename)) # TODO
        return mesh

    def __repr__(self):
        if not self.meta:
            r = f"WaiweraRun [{self.timestamp}]"
        else:
            r = f"WaiweraRun [{self.meta['run_name']}]"

        # r = f"WaiweraRun
        return r
    
    @staticmethod
    def get_values(value: str, filename: str, timestep: int) -> np.array:
        with h5py.File(filename) as data:
            cell_idx = [int(i[0]) for i in data["cell_index"]]
            temperature = data["cell_fields"][value][timestep]
            sorted_temp = np.zeros(temperature.shape)
            for i, idx in enumerate(cell_idx):
                sorted_temp[i] = temperature[idx]
        return sorted_temp

    @staticmethod
    def get_temperature(filename: str, timestep: int):
        return WaiweraRun.get_values("fluid_temperature", filename, timestep)
    
    @staticmethod
    def get_pressure(filename: str, timestep: int):
        return WaiweraRun.get_values("fluid_pressure", filename, timestep)

    @contextlib.contextmanager
    def open(self):
        data = h5py.File(self.h5_path)
        try:
            yield data
        finally:
            data.close()


class NullOutput:
    def write(self, x):
        pass

    def flush(self):
        pass


@contextlib.contextmanager
def nostdout():
    keep = sys.stdout
    sys.stdout = NullOutput()
    yield
    sys.stdout = keep
