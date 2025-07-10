# MECHBench: A Set of Optimization Benchmarks inspired by Structural Mechanics

## Overview

This repository provides an easy-to-use framework for evaluating optimization algorithms on structural mechanics problems using [OpenRadioss](https://openradioss.org). The optimization problem setup is defined in `main.py`, which controls simulation calls, design variables, objectives, and constraints.

We aim to provide a framework for testing algorithms and models on various structural optimization problems. It includes modules for defining optimization problems, generating finite element meshes, solving Explicit Finite Element Method (FEM) models, and post-processing results.


## Mechanical Test Cases

The repository currently contains three mechanical benchmark problems: a star-shaped crash box, a three-point bending problem, and a crash tube with trigger optimization. Each of these numerical models provides multiple outputs, enabling the definition of various optimization tasks (single/multi-objective, unconstrained/constrained scenarios). This diversity allows researchers to evaluate optimization algorithms on problems with different characteristics, objectives, and constraints within structural mechanics.

<table>
  <tr>
    <td><img src="starbox_diagram.png" alt="Star Box Problem" width="300"/></td>
    <td><img src="three_point_bending_diagram.png" alt="Three Point Bending Problem" width="300"/></td>
    <td><img src="crashtube_diagram.png" alt="Crash Tube Problem" width="300"/></td>
  </tr>
  <tr>
    <td>Starbox (Problem 1)</td>
    <td>Three-Point Bending (Problem 2)</td>
    <td>CrashTube (Problem 3)</td>
  </tr>
</table>




## Repository Structure

The project is structured as follows:

```
MECHBench/
├── main.py                      # Main optimization problem script
├── results/                     # Example output results
├── src/
│   └── sob/                        # Directory containing the core modules of the project
│       ├── __init__.py                # Initializes the Python package
│       ├── problems.py                # Defines optimization problem classes
│       ├── fem.py                     # Implements finite element analysis functionality
│       ├── mesh.py                    # Handles mesh generation and manipulation
│       ├── solver.py                  # A function to perform calls to OpenRadioss (starter and engine)
│       ├── post_processing.py         # Provides tools for post-processing FEM simulation results
│       ├── utils/                     # A directory containing utilities for setting up the simulations
│       │   ├── run_openradioss.py       # Contains RunOpenRadioss class to generates the commands to run OpenRadioss as well as external executables (TH_TO_CSV,ANIM_TO_VTK)
│       │   └── solver_setup.py          # Contains RunnerOptions class for simulation setup configuration. Extension of a dictionary to set up the runs.
│       └── lib/                       # A directory containing supplementary files and libraries
│           ├── py_mesh.py                   # Mesh functions for starbox model (Deprecated)
│           ├── py_mesh_v2.py                # Modified py_mesh for crash tube model, trigger generation variables added (Deprecated)
│           ├── gmsh_base_meshes.py          # Base class to generate GMSH meshes via __call__
│           ├── starbox_gmsh.py              # Mesh generation for starbox model with GMSH framework
│           ├── crashtube_gmsh.py            # Mesh generation for crash tube model with GMSH framework
│           └── three_point_bending_gmsh.py  # Mesh generation for three-point bending model with GMSH framework
├── tests.py                      # Contains multiple functions for checking the module functionality
├── README.md                     # Project documentation (this file)
└── requirements.txt              # Python dependencies
```

## Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/BayesOptApp/MECHBench.git
cd MECHBench
```
### 2. Install dependencies
Create a new environment (recommended) and install required packages:
```bash
pip install -r requirements.txt
```
:warning: Ensure that OpenRadioss is installed and accessible from your environment. For installation instructions, refer to [OpenRadioss Documentation](https://github.com/OpenRadioss/OpenRadioss/blob/main/INSTALL.md). There's also an official [YouTube video](https://www.youtube.com/watch?v=ddvH2CNfaYg&list=LL&index=4) for building it from source code. You can either download the stable version or build the binaries yourself. You should have the tools `TH_TO_CSV` and `ANIM_TO_VTK` in the directory `../OpenRadioss/exec` as well as the starter and engine executables.
Also, the `main.py` should point the path to the OpenRadioss directory in order to run the code (further instructions later).


## Usage

The optimization setup is controlled via `main.py`. Modify this file to:

- Set the path to your local OpenRadioss installation.
- Select the test case and the dimensionality of your problem.
- Configure objectives and constraints as needed for your problem.
- Specify Runner options and designated output folder.

### Running the benchmark

```bash
python main.py
```
This will:

1. Initialize the optimization problem.
2. Evaluate the objective function(s) by running OpenRadioss simulations.
3. Store the results in the designated output folder.

### Problem Setup in `main.py`

Below is a detailed walkthrough of the execution module `main.py`.

#### Operating System Detection

```python
if linux_system:
    orss_main_path = "/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/"
else:
    orss_main_path = "C:/Users/iolar/Downloads/OpenRadioss_win64_2/OpenRadioss/"
```
This code snippet sets the path to your OpenRadioss installation folder based on the operating system.


#### Runner Options

```python
runnerOptions = {
    "open_radioss_main_path": orss_main_path,
    "write_vtk": True,
    "np": 4,
    "nt": 1,
    "h_level": 1,
    "gmsh_verbosity": 0,
}
```

Configures how the simulation is run:

| Option               | Description                                    |
|-----------------------|------------------------------------------------|
| `open_radioss_main_path` | Path to OpenRadioss binaries.               |
| `write_vtk`             | Whether to output VTK files for visualization. |
| `np`                   | Number of processors for parallel runs.       |
| `nt`                   | Number of threads per process.                |
| `h_level`              | Mesh refinement level (if used).              |
| `gmsh_verbosity`       | Verbosity level for GMSH mesh generation.     |

#### Main function


Usage example of the main function:
```python
def main():
    sim_id = 253
    vector = [2.5, -1]
    f = sob.get_problem(1, 2, "intrusion", runnerOptions, sequential_id_numbering=False)
    obj_value = f(vector, sim_id)
    print(obj_value)
```

**Explanation:**

- `sim_id = 253`: Unique identifier for this simulation run. Results will be saved in folders named accordingly.
- `vector = [2.5, -1]`: Design variable vector to evaluate (must match problem dimensionality).
- `sob.get_problem(1, 2, "intrusion", runnerOptions, sequential_id_numbering=False)`: Retrieves objective function handle f:
     - `1`: Problem index or category ID.
     - `2`: Number of design variables.
     - `"intrusion"`: Target metric (e.g., intrusion depth).
     - `runnerOptions`: Simulation configuration.
     - `sequential_id_numbering`: Uses sequential numbering if True.
- `obj_value = f(vector, sim_id)`: Evaluates the objective function by:
     - Running OpenRadioss simulation with the design vector.
     - Storing results linked to sim_id.
     - Returning the computed objective value(s).
- `print(obj_value)`: Prints the output for inspection.


## Contact constributors
- Ivan Olarte Rodriguez (ivan.olarte.rodriguez@liacs.leidenuniv.nl) 
- Maria Laura Santoni (maria-laura.santoni@lip6.fr) 
- Elena Raponi (e.raponi@liacs.leidenuniv.nl)

## Other contributors
- Feifan Li (feifan.li@tum.de)

