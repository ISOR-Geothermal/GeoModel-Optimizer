# Borehole location optimizer

## Introduction

TODO

## Setup

The GeoModel Borehole optimizer package can be installed in a Python environment with pip:
pip install geomodel_optimizer

This will install all the required packages. It is recommended to do this in a virtual environment.
This package uses the pywaiwera package to run Waiwera using a Docker image. A Docker service must be installed and running. Docker desktop is recommended to those unfamiliar with Docker.

## Preparation

To use the optimizer, the following inputs are required:

**Mesh**: This module uses the layermesh package to work with the simulation mesh. For compatibility reasons two different mesh formats are required. The mesh will be saved as both a HDF5 (.h5) file as well as an EXODUS format file (.exo). The h5 file will be loaded by layermesh and used in python operations, while the .exo file will be read by Waiwera to run simulations. 

**Natural state**: Prior to the optimization step, the natural equilibrium state of the region needs to  be found. This is achieved using Waiwera, but the exact process is outside the scope of this documentation. For instructions on how to achieve this see the Waiwera documentation [https://waiwera.readthedocs.io/en/latest/]. The natural state serves as the starting point for optimization, representing the untouched geothermal state prior to any utilization. The output of the natural state calculation (h5 file) serves as the input for the optimization runs. 

**Base parameters**: This is a JSON file containing the parameters that resulted in the natural state. Copies of this file will be modified to serve as input files for individual optimization runs. This parameter file references the mesh file used to obtain the natural state by path, but the operator does not need to worry about fixing this path when moving to a new work directory for optimization, as the LocationOptimizer class does this automatically.

## Simple experiment: Moving sink

Having obtained a result for the natural state in our defined mesh with some natural background sources we can now move on to the optimization process. We will first go through an example to familiarize ourselves with the main methods of the LocationOptimizer class. 
The goal of the simulation is to compute the total output of a production well (sink) in multiple locations with a fixed reinjection well (source).

To start with, import the Optimizer class and instantiate it:
```
from geomodel_optimizer import LocationOptimizer

opt = LocationOptimizer(
    base_input='./natural_state.json',
    workdir='run_1',
    mesh_path_h5='./mesh1.h5',
    mesh_path_exo='./mesh1.exo',
    initial_state_path='./natual_state.h5'
)
```
We have supplied the constructor with several required arguments:
- `base_input`: The Waiwera parameters file. This can also be supplied as a dictionary. This is mostly the same file as the parameters used to find the natural state, except some changes might be made to the output parameters and timestepping
- `workdir`: the script will create this directory and set up the directory structure of the optimizer outputs
- `mesh_path_h5`: The mesh in HDF5 (.h5) format, this will be read by layermesh to be used by Python
- `mesh_path_exo`: The mesh in an EXODUS (.exo) format, this will be used by Waiwera
- `initial_state_path`: The output of the Waiwera run that computed the natural state in HDF5 format
When this constructor is called, all the input files will be copied to the working directory. This is because when Waiwera is called it will mount the work directory as a volume inside a container, so all Waiwera inputs must be contained in a subdirectory thereof.

The directory structure should look like this:
TODO
The goal is to simulate 50 years of extraction. The timesteps can be modified with the adjust_timesteps method:
```
opt.adjust_timesteps(step_size=opt.YEAR_IN_SECONDS, n_steps=50)
```
The step size is given in seconds, and for convenience the LocationOptimizer has a class attribute containing the number of seconds in a year.

The reinjection well will be stationary in all tests, and can be added with the permanent_source method:
```
opt.permanent_source(location=[1200, 1200, -1500], rate=100, enthalpy=160e3)
```
The location of the well can be specified with either a 3D coordinate or the index of a cell in the mesh. An analogous method permanent_sink also exists.

The previous two actions have modified the base parameters that will be common to all simulations, and so must be performed before adding any runs.

For this example we have predefined some sink locations in a file `sinks.txt`. We want to test each of the locations defined in the file, so each row will result in a Waiwera simulation being run with a sink at those coordinates. The file looks like this:
```
150, 150, -1200
650, 350, -500
# 1150, 150, -400
650, 850, -800
```
Each line represents X, Y, Z coordinates in the mesh in meters. We have commented out line 3 with a "#" sign, it will not be processed. Note the negative sign for the depth. 
To add these runs to the schedule, call the moving_sink method:
```
opt.moving_sink(rate=100, coordinates_file='sinks.txt')
```
The rate has been specified as 100 kg/s for every location. If you want to have different rates, you will need to use a different method shown later. 
We can see that the runs have been added to the `opt.sink_terms` attribute. 

```
>>> opt.sink_terms
[{'sink_idx': 2205,
  'sink_rate': 100,
  'sink_coords': array([  250.,   250., -1100.]),
  'sink_name': 'sink_2205'},
 {'sink_idx': 883,
  'sink_rate': 100,
  'sink_coords': array([ 700.,  250., -500.]),
  'sink_name': 'sink_883'},
 {'sink_idx': 1345,
  'sink_rate': 100,
  'sink_coords': array([ 700.,  700., -700.]),
  'sink_name': 'sink_1345'}]
```
The `sink_idx` variable refers to the cell index in the mesh. 
An analogous method moving_source also exists.

Since these are all the runs we wish to do at this time, we can now call the execute method:
```
opt.execute()
```
This method outputs all the optimization run parameters we have requested and starts running them with Waiwera through Docker. This step might take a while, especially if many runs are requested. 
If we want to inspect the individual parameter files for each run before starting computation, the output_run_files method will create all the files without starting the simulation. 

## Outputs

After Waiwera has processed all input files, we can take a look at what our work directory looks like.

For every run that was performed we have a waiwera control file in JSON format, the Waiwera output in HDF5 format, and the log file in YAML format. The runs have a name based on the sink location. 
The LocationOptimizer instance stores all our runs in its `opt.runs` attribute. This is a list containing each run as a `WaiweraRun` instance. This is a utility class that helps interact with the numerical output from Waiwera. 

### WaiweraRun

This class represents one run of Waiwera, keeping tabs on the relevant output and input files as well as some supplementary information. It is not necessary to work directly with runs as instances of this class, but knowing its structure could be helpful for inspecting single runs. The most important properties of this class are outlined below. 

Methods:
- `plots(value, depth)`: Returns a plot of the property specified by value (either pressure or temperature (default)). A depth slice and XY-slice are plotted. The layer to plot on the XY-slice is chosen by the depth parameter.
- `execute_function(func)`: Applies the function func to the data contained within the HDF5 file. The function must take a h5py dataset as input. 
- `get_values(value, filename, timestep)`: Read values of type value at timestep timestep from the file. This method takes care to reorder the values according to the cell index. This is important if multiple processes are used during computation as the cells will not be in the correct order. 
- `get_temperature(filename, timestep)`: Convenience method to call the get_values method for temperature
- `get_pressure`(filename, timestep): Convenience method to call the get_values method for pressure


Attributes:
- `temperature`: Upon construction the class will read the final temperature distribution from the file and set it as the temperature attribute. Only the final step is kept in memory to keep memory costs down
- `pressure`: The pressure distribution of the final step
- `meta`: If supplied, some metadata can be stored such as the sink location and properties that were special to this run
- `params`: a dict object containing the run parameters
- `json_file`: Path to the JSON file used as input
- `mesh`: The mesh object
- `h5_path`: Path to the Waiwera output HDF5

## More complex optimizations
If we want to optimize with both a moving sink and a moving source, this can be accomplished with sequential calls to those methods:

```
opt.moving_sink(rate=100, coordinates_file='sinks.txt')
opt.moving_source(rate=100, coordinates_file='sources.txt')
```
After calling both methods, both the `opt.sink_terms` and `opt.source_terms` will be populated. When `output_run_files` or `execute` is called the product of both lists will be used to create the run files. If the file `sinks.txt` contains *N* locations and the file `sources.txt` contains *M* locations there will be *N*x*M* runs. 
If you want a more complex optimization with variable rates or very specific combinations of sink and source parameters you can use the `add_run` method of the LocationOptimizer. 