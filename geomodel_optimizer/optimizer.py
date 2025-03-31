from .utils import nostdout
from .waiwerarun import WaiweraRun
import layermesh.mesh as lm
from pandas import DataFrame
import json
import shutil
from typing import Tuple, Optional, List, Callable, Any
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
        self._meta: list[dict[str, Any]] | None = None

        self.source_terms: List[dict] = [] # info about sources added to base params
        self.sink_terms: List[dict] = [] # info about sinks added to base params
        self.standalone_runs: List[dict] = [] # runs added with add_run
        self.combined_terms: List[dict] = [] # product of sink and source terms, generated when run files are output

        self.run_params: List[dict] = [] # run parameters in dict form
        self.run_file_paths: List[str] = [] # path to json files directly run by waiwera
        self.runs: List[WaiweraRun] = [] # completed runs

        self._run_counter: int = 1 # for keeping track of runs added with `add_run`

        self.YEAR_IN_SECONDS = 365 * 24 * 3600

    def add_run(self, 
                source_idx: Optional[int] = None, 
                sink_idx: Optional[int] = None, 
                source_rate: Optional[float] = None,
                sink_rate: Optional[float] = None, 
                source_enthalpy: float = 84e3,
                source_component: str = "water",
                sink_component: str = "water",
                ) -> None:
        """
        Adds a run configuration to the simulation by defining a source, a sink, or both.

        Parameters:
        -----------
        source_idx : Optional[int], default=None
            The index of the source cell. Must be provided along with `source_rate` if a source is specified.
        sink_idx : Optional[int], default=None
            The index of the sink cell. Must be provided along with `sink_rate` if a sink is specified.
        source_rate : Optional[float], default=None
            The flow rate at the source. Required if `source_idx` is specified. Units are kg/s if component is "water", 
            or J/s if component is "energy"
        sink_rate : Optional[float], default=None
            The flow rate at the sink. Required if `sink_idx` is specified. Units are kg/s if component is "water", 
            or J/s if component is "energy"
        source_enthalpy : float, default=84e3
            The enthalpy of the fluid entering through the source (J/kg).
        source_component : str, default="water"
            Component type at the source, "water" or "energy".
        sink_component : str, default="water"
            Component type at the sink, "water" or "energy".

        Raises:
        -------
        ValueError:
            If neither `source_idx` nor `sink_idx` is provided.
        AssertionError:
            If `source_idx` is provided without `source_rate`, or `sink_idx` is provided without `sink_rate`.

        Notes:
        ------
        - The function deep copies `self.base_input` to create a new parameter set.
        - If a sink is specified, it is added using `self.add_sink()`.
        - If a source is specified, it is added using `self.add_source()`.
        - The modified parameters are appended to `self.run_params`.
        """
        if (source_idx is None) and (sink_idx is None):
            raise ValueError("Must have at least a source or a sink")
        
        params = deepcopy(self.base_input)
        if sink_idx is not None: # sink_idx can be 0
            assert sink_rate is not None, f"must provide a sink rate along with sink_idx {sink_idx}"
            self.add_sink(params=params, cell=sink_idx, rate=sink_rate, component=sink_component)
        
        if source_idx is not None:
            assert source_rate is not None, f"must provide a source rate along with source_idx {source_idx}"
            self.add_source(params=params, cell=source_idx, rate=source_rate, enthalpy=source_enthalpy, component=source_component)

        terms = vars()
        terms.pop("self")
        terms.pop("params")

        name = f"run_{self._run_counter}"
        terms["run_name"] = name
        self._run_counter += 1

        if source_idx:
            source_coords = self.mesh.cell[source_idx].centre
            terms["source_coords"] = source_coords
        if sink_idx:
            sink_coords = self.mesh.cell[sink_idx].centre
            terms["sink_coords"] = sink_coords
        
        self.standalone_runs.append(terms)

    def _commit_run(self, 
                source_idx: Optional[int] = None, 
                sink_idx: Optional[int] = None, 
                source_rate: Optional[float] = None,
                sink_rate: Optional[float] = None, 
                source_enthalpy: float = 84e3,
                source_component: str = "water",
                sink_component: str = "water",
                **unused
                ) -> None:
        if (source_idx is None) and (sink_idx is None):
            raise ValueError("Must have at least a source or a sink")
        
        params = deepcopy(self.base_input)
        if sink_idx is not None: # sink_idx can be 0
            assert sink_rate is not None, f"must provide a sink rate along with sink_idx {sink_idx}"
            self.add_sink(params=params, cell=sink_idx, rate=sink_rate, component=sink_component)
        
        if source_idx is not None:
            assert source_rate is not None, f"must provide a source rate along with source_idx {source_idx}"
            self.add_source(params=params, cell=source_idx, rate=source_rate, enthalpy=source_enthalpy, component=source_component)
        
        self.run_params.append(params)

    def add_sink(self, 
                 params: dict, 
                 rate: float, 
                 cell: Optional[int] = None, 
                 location: Optional[Tuple[int]] = None,
                 component: str = "water",
                 ) -> None:
        """
        Adds a sink to the simulation at a specified cell or location.

        Modifies the `params` dictionary in place by appending a new sink entry.

        Parameters:
        -----------
        params : dict
            The simulation parameters dictionary that will be modified to include the sink.
        rate : float
            The flow rate at the sink. The function ensures this value is stored as a negative. Units are kg/s if component is "water", 
            or J/s if component is "energy"
        cell : Optional[int], default=None
            The index of the cell where the sink is added. If not provided, `location` must be specified.
        location : Optional[Tuple[int]], default=None
            The spatial coordinates of the sink. If `cell` is not given, the function attempts to find the cell index
            corresponding to this location using `self.mesh.find(location).index`.
        component : str, default="water"
            Component being extracted, "water" or "energy".

        Raises:
        -------
        ValueError:
            If neither `cell` nor `location` is provided.
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
        Adds a moving sink to the simulation, defined by a list of cells or coordinates from a file.

        Parameters:
        -----------
        rate : float
            The flow rate of the sink. Units are kg/s if component is "water", or J/s if component is "energy".
        cells : Optional[List[int]], default=None
            A list of cell indices where a sink will be simulated. Either this or `coordinates_file` must be provided.
        coordinates_file : Optional[str], default=None
            A file containing spatial coordinates for determining the sink locations. Overrides `cells` argument.
        component : str, default="water"
            "water" or "energy"

        Raises:
        -------
        ValueError:
            If neither `cells` nor `coordinates_file` is provided.

        Notes:
        ------
        - If `coordinates_file` is given, the function uses `self._cells_from_file()` to extract the corresponding cells.
        - The function iterates through the `cells` list, adding metadata for each sink to `self.sink_terms`.
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
        """
        Adds a permanent sink to the simulation.

        This function calls `self.add_sink()` to modify `self.base_input` in place, ensuring that the sink remains
        for the entire duration of the simulation.

        Parameters:
        -----------
        rate : float
            The flow rate of the sink. Units are kg/s if component is "water", or J/s if component is "energy".
        cell : Optional[int], default=None
            The index of the cell where the sink is added. If not provided, `location` must be specified.
        location : Optional[Tuple[int]], default=None
            The spatial coordinates of the sink. If `cell` is not given, the function attempts to determine the cell
            from `location` using `self.mesh.find(location)`.
        component : str, default="water"
            "water" or "energy"

        Returns:
        --------
        bool
            Always returns `True` after successfully adding the sink.
        """
        self.add_sink(cell=cell, location=location, rate=rate, params=self.base_input, component=component)
        return True

    def permanent_source(self, 
                         rate: float, 
                         cell: Optional[int] = None, 
                         location: Optional[Tuple[int]] = None,
                         enthalpy: float = 84.9e3, 
                         component: str = "water"
                         ) -> bool:
        """
        Adds a permanent source to the simulation.

        This function calls `self.add_source()` to modify `self.base_input` in place, ensuring that the source remains
        for the entire duration of the simulation.

        Parameters:
        -----------
        rate : float
            The flow rate of the source. Units are kg/s if component is "water", or J/s if component is "energy".
        cell : Optional[int], default=None
            The index of the cell where the source is added. If not provided, `location` must be specified.
        location : Optional[Tuple[int]], default=None
            The spatial coordinates of the source. If `cell` is not given, the function attempts to determine the cell
            from `location` using `self.mesh.find(location)`.
        component : str, default="water"
            "water" or "energy"

        Returns:
        --------
        bool
            Always returns `True` after successfully adding the source.
        """
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
        """
        Adds a source to the simulation at a specified cell or location.

        Modifies the `params` dictionary in place by appending a new source entry.

        Parameters:
        -----------
        params : dict
            The simulation parameters dictionary that will be modified to include the source.
        rate : float
            The flow rate of the source. Units are kg/s if component is "water", or J/s if component is "energy".
        cell : Optional[int], default=None
            The index of the cell where the source is added. If not provided, `location` must be specified.
        location : Optional[Tuple[int]], default=None
            The spatial coordinates of the source. If `cell` is not given, the function attempts to find the cell index
            corresponding to this location using `self.mesh.find(location)`.
        component : str, default="water"
            "water" or "energy"

        Raises:
        -------
        ValueError:
            If neither `cell` nor `location` is provided.
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
        Adds a moving source to the simulation, defined by a list of cells or coordinates from a file.

        Parameters:
        -----------
        rate : float
            The flow rate of the source. Units are kg/s if component is "water", or J/s if component is "energy".
        cells : Optional[List[int]], default=None
            A list of cell indices where the source will be placed. Either this or `coordinates_file` must be provided.
        coordinates_file : Optional[str], default=None
            A file containing spatial coordinates for determining the cells of the sources.
            If provided, the function reads this file to determine the `cells` list.
        component : str, default="water"
            "water" or "energy"

        Raises:
        -------
        ValueError:
            If neither `cells` nor `coordinates_file` is provided.

        Notes:
        ------
        - If `coordinates_file` is given, the function uses `self._cells_from_file()` to extract the corresponding cells.
        - The function iterates through the `cells` list, adding metadata for each source to `self.source_terms`.
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
        """
        cells = []
        cell_lines = {}
        with open(file_path) as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if "#" in line:
                continue

            try:
                line_coords = tuple([float(i) for i in line.split(",")])
                assert len(line_coords) == 3
            except:
                warnings.warn(f"could not convert line {i+1} to 3D coordinates: {line}")
                continue
                
            cell = self.mesh.find(line_coords)
            if cell is None:
                warnings.warn(f"could not match coordinates from line {i+1} [{line}] with a cell in the mesh")
                continue

            if cell.index in cells:
                line_number = cell_lines[cell.index]
                warnings.warn(f"cell matching location in line {i+1} is already included from line {line_number}, ignoring")
                continue

            cells.append(cell.index)
            cell_lines[cell.index] = i + 1 # save the line number 

        return cells
        
    def _make_dirs(self, workdir: str):
        wd = Path(workdir)
        if wd.is_absolute():
            raise ValueError("workdir path cannot be absolute, must be a subdir")
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
            raise FileNotFoundError(f"could not find mesh .exo file at {mesh_exo_path}")
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
            raise FileNotFoundError(f"could not find initial state file at {initial_state_path}")
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
            return self.sink_terms + self.standalone_runs
        if len(self.sink_terms) == 0:
            return self.source_terms + self.standalone_runs
        
        # else we need to create a product of all terms
        combined_terms = []
        for source, sink in itertools.product(self.source_terms, self.sink_terms):
            combined = source.copy()
            combined.update(sink)
            combined_terms.append(combined)

        combined_terms += self.standalone_runs

        return combined_terms

    def _get_run_name(self, meta: dict) -> str:
        predefined_name = meta.get("run_name", None)
        if predefined_name:
            return predefined_name

        sink_name = meta.get("sink_name", "")
        source_name = meta.get("source_name", "")

        name = f"{sink_name}_{source_name}"
        if name.startswith("_"):
            name = name[1:]
        if name.endswith("_"):
            name = name[:-1]
        if name == "":
            raise ValueError("failed to assign name")
        return name

    def output_run_files(self):
        """
        Set the run parameters outputs and output the parameters to a json file to be 
        used by Waiwera. Also creates the self.meta DataFrame that holds information
        about run parameters
        """
        self.combined_terms = self._get_combined_terms() # populates self.combined_terms
        for terms in self.combined_terms:
            self._commit_run(**terms) # adds to self.run_params

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
        self._meta = metadata
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

    def list_run_output_files(self) -> List[Path]:
        """
        Returns a list of the Waiwera HDF5 output files (if any)
        """
        return list(self.output_dir.glob("*.h5"))

    def compute_loss(self, loss_function: Callable):
        """
        Applies `loss_function` to every run computed during the simulation phase. The function
        should have a signature like `loss_function(run: WaiweraRun) -> float`. The return values 
        are collected to a list and finally added to the `self.meta` DataFrame.

        Parameters:
        -----------
        loss_function : Callable
            a function with a signature like `f(run: WaiweraRun) -> float`
        """
        if self.meta is None:
            raise ValueError("No metadata found. Has .execute() been called?")
        calculated_loss = []
        for run in tqdm(self.runs):
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
        
        Parameters:
        -----------
        nproc : int, default=4
            Number of processes used by Waiwera
        """
        self.output_run_files()
        self._sequential_run(nproc=nproc)
        self.meta.to_csv(f"{self.workdir}/meta.csv")

    def adjust_timesteps(self, step_size: int, n_steps: int) -> None:
        """
        Adjust the timesteps for the production simulation. Disables adaptable timestepping.

        Parameters:
        -----------
        step_size : int, time step in seconds
        n_steps : int, number of time steps
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

    def _select(self, cond) -> DataFrame:
        keep = self.meta.mask(cond).isna()
        keep = keep.json_file
        return self.meta[keep]

    def run_by_sink_index(self, index: int) -> List[WaiweraRun]:
        select = self.meta[self.meta.sink_idx.eq(index)]
        if select.empty:
            return []

        runs = []
        for item in select.itertuples():
            runs.append(self.runs[item.run_index])
        return runs

    @classmethod
    def from_workdir(cls, workdir: str):
        raise NotImplemented

