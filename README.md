# MECHBench: A Set of Optimization Benchmarks inspired by Structural Mechanics

## Overview

This repository provides an easy-to-use framework for evaluating optimization algorithms on structural mechanics problems using [OpenRadioss](https://openradioss.org). The optimization problem setup is defined in `main.py`, which controls simulation calls, design variables, objectives, and constraints.

We aim to provide a framework for testing algorithms and models on various structural optimization problems. It includes modules for defining optimization problems, generating finite element meshes, solving Explicit Finite Element Method (FEM) models, and post-processing results.


## Mechanical Test Cases

The repository currently contains three mechanical benchmark problems: a star-shaped crash box, a three-point bending problem, and a crash tube with trigger optimization. Each of these numerical models provides multiple outputs, enabling the definition of various optimization tasks (single/multi-objective, unconstrained/constrained scenarios). This diversity allows researchers to evaluate optimization algorithms on problems with different characteristics, objectives, and constraints within structural mechanics.

<table>
  <tr>
    <td><img src="figures/starbox_diagram.png" alt="Star Box Problem" width="300"/></td>
    <td><img src="figures/three_point_bending_diagram.png" alt="Three Point Bending Problem" width="300"/></td>
    <td><img src="figures/crashtube_diagram.png" alt="Crash Tube Problem" width="300"/></td>
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
├── main.py
├── main_sequential.py
├── main_experiments.py
├── results/
├── src/
│   └── sob/
│       ├── __init__.py
│       ├── observer.py
│       ├── sampler.py
│
│       ├── physical_models/
│       │   ├── __init__.py
│       │
│       │   # -------------------------
│       │   # Core physical models
│       │   # -------------------------
│       │   ├── abstractPhysicalModel.py
│       │   ├── crashTube.py
│       │   ├── starBox.py
│       │   └── threePointBending.py
│
│       │   # -------------------------
│       │   # FEM SETTINGS
│       │   # -------------------------
│       │   ├── fem_settings/
│       │   │   ├── __init__.py
│       │   │   ├── abstractFEMSettings.py
│       │   │   ├── crashTubeModel.py
│       │   │   ├── starBoxModel.py
│       │   │   └── threePointBendingModel.py
│
│       │   # -------------------------
│       │   # MESHES
│       │   # -------------------------
│       │   ├── meshes/
│       │   │   ├── __init__.py
│       │   │   ├── abstractMeshSettings.py
│       │   │   ├── crashTubeMesh.py
│       │   │   ├── starBoxMesh.py
│       │   │   ├── threePointBendingMesh.py
│       │   │   └── routines/              # (new submodule)
│       │   │       └── gmsh/
│       |   |           ├── crashtube_gmsh.py
│       │   │           ├── gmsh_base_meshes.py
│       │   │           ├── starbox_gmsh.py
│       │   │           └── three_point_bending_gmsh.py
│       │   │       ├── legacy/
│       │   |           ├── py_mesh.py
│       │   |           └── py_mesh_v2.py
│       │   # -------------------------
│       │   # SOLVERS
│       │   # -------------------------
│       │   ├── solvers/
│       │   │   └── openRadioss_runner.py
│
│       │   # -------------------------
│       │   # UTILS (model-specific)
│       │   # -------------------------
│       │   ├── utils/
│       │   │   ├── platform_det.py
│       │   │   ├── run_openradioss.py
│       │   │   └── solver_setup.py
│
│       │   # -------------------------
│       │   # LIB (input decks, etc.)
│       │   # -------------------------
│       │   ├── lib/
│       │   │   ├── ThreePointBending_0001.rad
│       │   │   └── ThreePointBending_base.rad
│
│       ├── problems/
│       │   ├── IOH_wrappers/
│       │   │   ├── __init__.py
│       │   │   ├── constrained_single_objective.py
│       │   │   ├── real_constraint.py
│       │   │   └── single_objective.py
│       │   └── __init__.py
│
├── tests.py
├── README.md
└── requirements.txt
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


## Usage

The optimization setup is controlled via `main.py`. Modify this file to:

- Select the test case and the dimensionality of your problem.
- Configure objectives and constraints as needed for your problem.

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


#### ⚙️ Runner Options

The behavior of simulations in MECHBench is controlled through the `runnerOptions` dictionary, which is passed to `sob.get_problem(...)`.

---

##### Default Configuration

If no options are provided, MECHBench uses the following defaults:

```python
runnerOptions = {
    "np": 1,
    "nt": 1,
    "h_level": 1,
    "gmsh_verbosity": 0,
    "write_vtk": False,
}
```

These settings provide a stable, low-resource baseline that works on most systems.

##### Basic Usage
You can override the defaults by defining your own configuration:

```python
runnerOptions = {
    "np": 1,
    "nt": 8,
    "h_level": 1,
    "gmsh_verbosity": 0}
```

Then pass it to the problem
```python
f = sob.get_problem(problem_id, dim, runner_options=runnerOptions)
```

##### Available Options
| Option                                | Type | Default | Description                                                            |
| ------------------------------------- | ---- | ------- | ---------------------------------------------------------------------- |
| `np`                                  | int  | 1       | Number of processes used for the simulation (MPI parallelism).         |
| `nt`                                  | int  | 1       | Number of threads per process (OpenMP parallelism).                    |
| `h_level`                             | int  | 1       | Mesh refinement level. Higher values increase resolution (and cost).   |
| `gmsh_verbosity`                      | int  | 0       | Controls GMSH output verbosity (0 = silent).                           |
| `write_vtk`                           | bool | False   | If `True`, exports simulation results in VTK format for visualization. |
| `open_radioss_main_path` *(optional)* | str  | None    | Path to a local OpenRadioss installation.                              |

##### Automatic OpenRadioss Handling
You do not need to manually install OpenRadioss.

If `open_radioss_main_path` is not provided, MECHBench will:

1. Check if OpenRadioss binaries are available locally
2. If not found, automatically download the appropriate binaries

##### Advanced: Custom OpenRadioss Path
For full control over the solver (e.g., custom builds or HPC environments):
```python
runnerOptions = {
    "open_radioss_main_path": "/path/to/OpenRadioss/",
    "np": 4,
    "nt": 2,
}
```

##### Performance Tips
- Increase np for multi-process execution (MPI)
- Increase nt for multi-threading (OpenMP)
- Use higher h_level for finer meshes (at higher computational cost)
- Keep write_vtk=False unless visualization is needed

#### Main function


The main entry point for running a simulation or evaluating a design is defined in `main.py`.

Below is a minimal working example:
```python
def main():
    sim_id = 238
    dim = 5
    problem_id = 3  # 1: star box, 2: three point bending, 3: crash tube

    vector = np.random.uniform(0, 0, (dim,)).tolist()
    print(f"Evaluating vector: {vector}")

    f = sob.get_problem(problem_id, dim, runner_options=runnerOptions)
    obj_value = f(vector, sim_id)

    print(obj_value)
```

**Explanation:**

- `sim_id`
Unique identifier for the simulation.
Results are stored in folders associated with this ID, allowing reproducibility and tracking.
- `problem_id`
Selects the mechanical benchmark problem:
- 1 → Star Box
- 2 → Three-Point Bending
- 3 → Crash Tube
- `dim`
Number of design variables.
This must match the dimensionality expected by the selected problem.
- `vector`
Design variable vector to evaluate.
It can be randomly generated (as shown) or manually defined.

##### Creating the objective function
```python
f = sob.get_problem(problem_id, dim, runner_options=runnerOptions)
```

This function:

- Instantiates the selected _physical model_
- Prepares:
- - Mesh generation
- - FEM configuration
- - Solver pipeline
- Returns a callable function `f`

##### Evaluating a design
```python
obj_value = f(vector, sim_id)
```

This step:

1. Takes the input design vector
2. Generates the corresponding mesh and simulation setup
3. Runs the OpenRadioss simulation (if required)
4. Post-processes results
5. Stores outputs using sim_id
6. Returns the computed objective value(s)

###### Notes
- The model is initialized lazily: mesh and FEM data are only generated when the function `f` is called with a vector.
- The returned value can be:
  - A single scalar (single-objective)
  - A tuple/list (multi-objective or multiple metrics)

## Metrics to Extract

In the following table you may find all the possible metrics you can extract for the problems. You can either input a sole target metric or define a list (or tuple) with the names of the metrics to extract multiple. With this in mind, if you want to define different optimization problems (extra constraints, or multi-objective) you can extract the variables/metrics in order to tailor your own problem.

| Metric                    | Requires FEM Simulation | Notes                                                              |
|---------------------------|--------------------------|--------------------------------------------------------------------|
| `mass`                    | No                       |                                                                    |
| `absorbed_energy`         | No                       |                                                                    |
| `intrusion`               | Yes                      | Requires running FEM simulation                                    |
| `specific_energy_absorbed`| No                       | Computed from `mass` and `absorbed_energy`                         |
| `mean_impact_force`       | Yes                      | Requires running FEM simulation                                    |
| `peak_impact_force`       | Yes                      | Requires running FEM simulation                                    |
| `usage_ratio`             | Yes                      | Requires running FEM simulation, Computed from `mean_impact_force`  and `peak_impact_force`  |
| `load_uniformity`         | Yes                      | Requires running FEM simulation, Computed from `mean_impact_force`  and `peak_impact_force`  |
| `penalized_sea`         | Yes                      | Requires running FEM simulation, Computed from `mass` and `absorbed_energy`, only available for Problem 1  |
| `penalized_mass`         | Yes                      | Requires running FEM simulation, Computed from `mass` and `absorbed_energy`, only available for Problem 2  |

### ⚠️ Important: Force Computation Differences (Windows vs Linux)

When working with force-related metrics (e.g., `mean_impact_force`, `peak_impact_force`, `load_uniformity`, `usage_ratio`), please be aware of a platform-dependent discrepancy between OpenRadioss binaries:

#### Key Difference

- **Linux binaries**  
  → Directly output **force values**

- **Windows binaries**  
  → Output **impulse**, which is then internally converted to force

---

#### Implications

- Results involving force-related metrics **may differ slightly** between Windows and Linux systems  
- This discrepancy arises from:
  - Numerical differentiation of impulse (Windows)
  - Direct force extraction (Linux)
- The difference is typically small but can be noticeable in:
  - Peak force values
  - Derived metrics (e.g., load uniformity, usage ratio)

---

#### Recommendations

- For **consistent benchmarking**, use the **same operating system across all experiments**
- Prefer **Linux environments** when possible, as they provide direct force outputs
- Avoid mixing results from Windows and Linux runs in the same analysis

---

#### Affected Metrics

- `mean_impact_force`
- `peak_impact_force`
- `usage_ratio`
- `load_uniformity`

---

This behavior is inherent to the OpenRadioss binaries and not specific to MECHBench.

## Contact constributors
- Ivan Olarte Rodriguez (ivan.olarte.rodriguez@liacs.leidenuniv.nl) 
- Maria Laura Santoni (maria-laura.santoni@lip6.fr) 
- Elena Raponi (e.raponi@liacs.leidenuniv.nl)

## Other contributors
- Feifan Li (feifan.li@tum.de)

