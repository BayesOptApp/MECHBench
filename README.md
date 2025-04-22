# Structural Optimization Benchmark Project

This repository contains a structural optimization benchmark project using OpenRadioss as solver.

## Overview

The structural optimization benchmark project aims to provide a framework for implementing and evaluating various structural optimization problems. It includes modules for defining optimization problems, generating finite element meshes, solving FE models, and post-processing results.

## Project Structure

The project is structured as follows:

- `tests.py`: Containing multiple functions for checking the module functionality.
- `src/sob/`: The directory containing the core modules of the project.
  - `__init__.py`: Initializes the Python package.
  - `problems.py`: Defines optimization problem classes.
  - `fem.py`: Implements finite element analysis functionality.
  - `mesh.py`: Handles mesh generation and manipulation.
  - `solver.py`: OpenRadioss.
  - `post_processing.py`: Provides tools for post-processing FEM simulation results.
  - `lib/`: A directory containing supplementary files and libraries.
      - `py_mesh.py`: A file with functions related to mesh generation (for starbox model).
      - `py_mesh_v2.py`: Modified version of py_mesh_v2, trigger generation variables added (for crash tube model).

## Instructions for Running OpenRadioss Solver
1. You should have OpenRadioss installed. About how to install it, see: https://github.com/OpenRadioss/OpenRadioss/blob/main/INSTALL.md. There's also an official YouTube video for building it from sourse code: https://www.youtube.com/watch?v=ddvH2CNfaYg&list=LL&index=4
2. In order to run the OpenRadioss solver within the module, make sure to have the job submission tool folder downloaded in the root directory of your OpenRadioss installation. You can obtain this folder from the webpage: https://openradioss.atlassian.net/wiki/spaces/OPENRADIOSS/pages/58785793/Latest+gui+scripts

## Contact constributors
- Maria Laura Santoni
- Feifan Li (feifan.li@tum.de)
- Ivan Olarte Rodriguez (ivan.olarte.rodriguez@liacs.leidenuniv.nl)
