import layermesh.mesh as lm
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pywaiwera
import json
import shutil
from typing import Tuple, Optional, List, Callable
import os, sys
from pathlib import Path
from copy import deepcopy
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
import time
from pandas import Timestamp


def get_values(value: str, filename: str, timestep: int) -> np.array:
    with h5py.File(filename) as data:
        cell_idx = [int(i[0]) for i in data["cell_index"]]
        temperature = data["cell_fields"][value][timestep]
        sorted_temp = np.zeros(temperature.shape)
        for i, idx in enumerate(cell_idx):
            sorted_temp[i] = temperature[idx]
    return sorted_temp


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
