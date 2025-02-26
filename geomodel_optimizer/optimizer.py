from .utils import nostdout
from .waiwerarun import WaiweraRun
import layermesh.mesh as lm
from pandas import DataFrame
import json
import shutil
from typing import Tuple, Optional, List, Callable
from pathlib import Path
from copy import deepcopy
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools


class LocationOptimizer:
    def __init__(self, 
                 base_input: dict | str, 
                 workdir: str, 
                 mesh_path_h5: str,
                 mesh_path_exo: Optional[str] = None, 
                 initial_state_path: Optional[str] = None):
        self._make_dirs(workdir)
        self.mesh = lm.mesh(mesh_path_h5)
        if isinstance(base_input, str):
            with open(base_input, "r") as f:
                base_input = json.load(f)

        initial_base_input = deepcopy(base_input)

        # this call modifies the input to use files in the local directory
        base_input = self._prepare_inputs(initial_base_input, mesh_path_h5, mesh_path_exo, initial_state_path)
        self.base_input = base_input

        self.meta: Optional[DataFrame] = None # metadata about runs

        self.source_terms: List[dict] = [] # info about sources added to base params
        self.sink_terms: List[dict] = [] # info about sinks added to base params
        self.combined_terms: List[dict] = [] # product of sink and source terms, generated when run files are output

        self.run_params: List[dict] = [] # run parameters in dict form
        self.run_file_paths: List[str] = [] # path to json files directly run by waiwera
        self.runs: List[WaiweraRun] = [] # completed runs

        self.YEAR_IN_SECONDS = 365 * 24 * 3600

    def add_run(self,  # TODO add support for component terms
                source_idx: Optional[int] = None, 
                sink_idx: Optional[int] = None, 
                source_rate: Optional[float] = None,
                sink_rate: Optional[float] = None, 
                source_enthalpy: float = 84e3,
                **unused
                ) -> None:
        if (source_idx is None) and (sink_idx is None):
            raise ValueError("TODO both cant be none") # TODO raise warning instead?
        
        params = deepcopy(self.base_input)
        if sink_idx is not None: # sink_idx can be 0
            assert sink_rate is not None, f"must provide a sink rate along with sink_idx {sink_idx}"
            self.add_sink(params=params, cell=sink_idx, rate=sink_rate)
        
        if source_idx is not None:
            assert source_rate is not None, f"must provide a source rate along with source_idx {source_idx}"
            self.add_source(params=params, cell=source_idx, rate=source_rate, enthalpy=source_enthalpy)
        
        self.run_params.append(params)

    def add_sink(self, 
                 params: dict, 
                 rate: float, 
                 cell: Optional[int] = None, 
                 location: Optional[Tuple[int]] = None,
                 component: str = "water",
                 ) -> None:
        """
        Add a sink at the provided cell. Modifies `params` in place.
        """
        if (cell is None) and (location is None):
            raise ValueError("Must supply either cell or location for sink")
        if cell is None:
            cell = self.mesh.find(location).index
        try:
            sources = params["source"]
        except KeyError:
            sources = []
        new_source = {
            "cell": cell,
            "rate": -abs(rate), # ensure negative
            "component": component
        }
        sources.append(new_source)
        params["sources"] = sources

    def moving_sink(self,
                rate: float,
                cells: Optional[List[int]] = None,
                coordinates_file: Optional[str] = None,
                component: str = "water"
            ) -> None:
        """
        TODO
        Adds sink in every cell provided in `cells` list. 
        # Every iteration adds a parameter dict to the optimizer object's run_params list. 
        Every iteration adds an entry to the sink_terms list. When the run list is finalized the product of 
        the source_terms and sink_terms is used to generate the full list of runs
        """
        if (cells is None) and (coordinates_file is None):
            raise ValueError("Missing cell information. Provide either cells: List[int] or coordinates_file: str argument")
        if coordinates_file:
            cells = self._cells_from_file(coordinates_file)

        for i, sink_cell in enumerate(cells):
            sink_meta = {
                "sink_idx": sink_cell,
                "sink_rate": rate,
                "sink_coords": self.mesh.cell[sink_cell].centre,
                "sink_name": f"sink_{sink_cell}",
                "sink_component": component
            }
            self.sink_terms.append(sink_meta)

    def permanent_sink(self, 
                       rate: float, 
                       cell: Optional[int] = None, 
                       location: Optional[Tuple[int]] = None,
                       component: str = "water"
                       ) -> bool:
        self.add_sink(cell=cell, location=location, rate=rate, params=self.base_input, component=component)
        return True

    def permanent_source(self, 
                         rate: float, 
                         cell: Optional[int] = None, 
                         location: Optional[Tuple[int]] = None,
                         enthalpy: float = 84.9e3, 
                         component: str = "water"
                         ) -> bool:
        self.add_source(cell=cell, location=location, rate=rate, params=self.base_input, enthalpy=enthalpy, component=component)
        return True # TODO why this return

    def add_source(self, 
                   params: dict, 
                   rate: float, 
                   cell: Optional[int] = None, 
                   location: Optional[Tuple[int]] = None,
                   enthalpy=84.9e3, 
                   component="water"
                   ) -> None:
        if (cell is None) and (location is None):
            raise ValueError("Must supply either cell or location for sink")
        if cell is None:
            cell = self.mesh.find(location).index
        try:
            sources = params["source"]
        except KeyError:
            sources = []
        new_source = {
            "cell": cell,
            "component": component,
            "enthalpy": enthalpy,
            "rate": abs(rate)
        }
        sources.append(new_source)
        params["sources"] = sources        

    def moving_source(self, 
                rate: float,
                enthalpy: float = 83.9e3,
                cells: Optional[List[int]] = None,
                component: str = "water",
                coordinates_file: Optional[str] = None,
            ) -> None:
        """
        Add a simulation with a source in each of the cells defined in `cells`. If a file 
        is specified with `coordinates_file`, the cells will be read from there instead. The 
        file should contain three comma separated numbers per line, representing the X, Y and Z
        coordinates of the source in the mesh in meters.
        """
        if (cells is None) and (coordinates_file is None):
            raise ValueError("Missing cell information. Provide either cells: List[int] or coordinates_file: str argument")

        if coordinates_file:
            cells = self._cells_from_file(coordinates_file)

        for i, source_cell in enumerate(cells):
            source_meta = {
                "source_idx": source_cell,
                "source_rate": rate,
                "source_coords": self.mesh.cell[source_cell].centre,
                "source_enthalpy": enthalpy,
                "source_name": f"source_{source_cell}",
                "source_component": component
            }
            self.source_terms.append(source_meta)

    def _cells_from_file(self, file_path: str) -> List[int]:
        """
        Read coords from file and turn them into cell indexes. Each line in the input file should 
        contain one 3d coordinate. Ignore lines commented out with "#"

        TODO warn user about duplicate cells? Warn about coordinates mismatch? (cell.centre vs input)
        """
        coordinates: List[float] = []

        with open(file_path) as f:
            for line in f.readlines():
                if "#" in line:
                    continue
                line_coords = tuple([float(i) for i in line.split(",")])
                coordinates.append(line_coords)

        cell_indexes = self._cells_from_coords(coordinates)
        return cell_indexes

    def _cells_from_coords(self, coords: List[Tuple[float]]) -> List[int]:
        cells = []
        for c in coords:
            if len(c) != 3:
                raise ValueError("TODO need 3d coord")
            cell = self.mesh.find(c)
            cells.append(cell.index)
        return cells
        
    def _make_dirs(self, workdir: str):
        wd = Path(workdir)
        if wd.is_absolute():
            raise ValueError("TODO path cannot be absolute")
        wd.mkdir(exist_ok=True, parents=True)
        json_dir = wd / "json"
        output_dir = wd / "outputs"
        input_dir = wd / "inputs"
        log_dir = wd / "logs"

        json_dir.mkdir(exist_ok=True)
        output_dir.mkdir(exist_ok=True)
        input_dir.mkdir(exist_ok=True)
        log_dir.mkdir(exist_ok=True)

        self.workdir = wd
        self.json_dir = json_dir
        self.output_dir = output_dir
        self.input_dir = input_dir
        self.log_dir = log_dir

    def _prepare_inputs(self, parameters: dict, mesh_h5_path: str, mesh_exo_path: Optional[str] = None, initial_state_path: Optional[str] = None) -> dict:
        """
        Copy all relevant inputs defined in the parameters file to the working directory
        """
        if mesh_exo_path is None:
            mesh_exo_path = Path(parameters["mesh"]["filename"])
        else:
            mesh_exo_path = Path(mesh_exo_path)

        if not mesh_exo_path.exists():
            raise FileNotFoundError("TODO could not find mesh h5 file")
        shutil.copy(mesh_exo_path, self.input_dir)
        shutil.copy(mesh_h5_path, self.input_dir)

        relative_mesh_file = self.input_dir / mesh_exo_path.name
        relative_mesh_file = relative_mesh_file.as_posix() # need to use linux paths for docker image
        parameters["mesh"]["filename"] = relative_mesh_file

        if initial_state_path is None:
            initial_state_path = Path(parameters["initial"]["filename"])
        else:
            initial_state_path = Path(initial_state_path)
        
        if not initial_state_path.exists():
            raise FileNotFoundError("TODO could not find state file")
        shutil.copy(initial_state_path, self.input_dir)
        
        relative_state_file = self.input_dir / initial_state_path.name
        relative_state_file = relative_state_file.as_posix() # need to use linux paths for docker image
        parameters["initial"]["filename"] = relative_state_file # TODO check if original initial params need to be removed

        base_parameters_output = self.input_dir / "base_params.json"
        with base_parameters_output.open("w") as f:
            json.dump(parameters, f, indent=4)
        return parameters
        
    def _get_combined_terms(self) -> List[dict]:
        """
        If the optimization has both moving sinks and sources, the total number of runs is 
        the product of all source locations and all sink locations
        """
        # if only either source or sink terms exist nothing needs to be done
        if len(self.source_terms) == 0:
            return self.sink_terms
        if len(self.sink_terms) == 0:
            return self.source_terms
        
        # else we need to create a product of all terms
        combined_terms = []
        for source, sink in itertools.product(self.source_terms, self.sink_terms):
            combined = source.copy()
            combined.update(sink)
            combined_terms.append(combined)

        return combined_terms

    def _get_run_name(self, meta: dict) -> str:
        sink_name = meta.get("sink_name", "")
        source_name = meta.get("source_name", "")

        name = f"{sink_name}_{source_name}"
        if name.startswith("_"):
            name = name[1:]
        if name.endswith("_"):
            name = name[:-1]
        return name

    def output_run_files(self):
        """
        Set the run parameters outputs and output the parameters to a json file to be 
        used by Waiwera. Also creates the self.meta DataFrame that holds information
        about run parameters
        """
        self.combined_terms = self._get_combined_terms()
        for terms in self.combined_terms:
            self.add_run(**terms)

        run_file_paths = []
        metadata = [] # gather meta information to be made into a dataframe
        for i, (params, meta) in enumerate(zip(self.run_params, deepcopy(self.combined_terms))):
            name = self._get_run_name(meta)
            self._set_output(params, name)
            json_path = self.json_dir / (name + ".json")
            with json_path.open("w") as f:
                json.dump(params, f, indent=4)
            run_file_paths.append(str(json_path).replace("\\", "/"))

            meta.update({
                "json_file": str(json_path),
                "run_index": i,
                "run_name": name,
                "loss": None
            })
            metadata.append(meta)

        self.run_file_paths = run_file_paths
        self.meta = self._generate_info_df(metadata)

    def _generate_info_df(self, meta: List[dict]) -> DataFrame:
        df = DataFrame(meta)
        if "source_coords" in df.columns:
            df["source_x"] = df.source_coords.apply(lambda x: x[0])
            df["source_y"] = df.source_coords.apply(lambda x: x[1])
            df["source_z"] = df.source_coords.apply(lambda x: x[2])
            del df["source_coords"]
        
        if "sink_coords" in df.columns:
            df["sink_x"] = df.sink_coords.apply(lambda x: x[0])
            df["sink_y"] = df.sink_coords.apply(lambda x: x[1])
            df["sink_z"] = df.sink_coords.apply(lambda x: x[2])
            del df["sink_coords"]

        return df

    def _set_output(self, params: dict, basename: str):
        h5_name = self.output_dir / (basename + ".h5")
        log_name = self.log_dir / (basename + ".yaml")
        params["output"]["filename"] = h5_name.as_posix()
        params["logfile"]["filename"] = log_name.as_posix()

    def _sequential_run(self, nproc=4) -> bool:
        total = len(self.run_file_paths)
        for run_index, path in tqdm(enumerate(self.run_file_paths), total=total):
            run = self._single_run(path, nproc=nproc, run_index=run_index)
            self.runs.append(run)
        return True

    def _single_run(self, path: str, nproc: int = 4, run_index: Optional[int] = None) -> WaiweraRun:
        with nostdout():
            if run_index is not None:
                meta = self.meta.iloc[run_index].to_dict()
            else:
                meta = {}
            run = WaiweraRun.run_from_file(path, num_processes=nproc, meta=meta)
        return run

    def _parallel_run(self, pool_size=4, nproc=1):
        # TODO
        with tqdm(total=len(self.run_params)) as progress:
            with ProcessPoolExecutor(max_workers=pool_size) as ppe:
                futures = [ppe.submit(self._single_run, path, nproc) for path in self.run_file_paths]
                for future in as_completed(futures):
                    progress.update(1)

    def list_run_output_files(self) -> List[Path]: # TODO drop?
        return list(self.output_dir.glob("*.h5"))

    def compute_loss(self, loss_function: Callable):
        if self.meta is None:
            raise ValueError("TODO no meta")
        calculated_loss = []
        for run in self.runs:
            # loss = run.execute_function(loss_function)
            loss = loss_function(run)
            calculated_loss.append(loss)
        self.meta["loss"] = calculated_loss
    
    def __len__(self):
        return len(self.run_params)
    
    def execute(self, nproc=4):
        """
        Perform a run with each set of parameters in `self.run_params`. The runs can then be accessed 
        from the self.runs list
        """
        self.output_run_files()
        self._sequential_run(nproc=nproc)

    def adjust_timesteps(self, step_size: int, n_steps: int) -> None:
        """
        Adjust the timesteps for the production simulation. Disables adaptable timestepping.
        step_size : int time step in seconds
        n_steps : int number of time steps
        """
        time = {
            'step': {
                'size': step_size,
                'adapt': {'on': False},
                'maximum': {'number': n_steps}
            },
            'stop': step_size * n_steps
        }
        self.base_input["time"] = time
