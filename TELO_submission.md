# TELO submission environment

This guide explains how to install MECHBench and the optional dependency groups
declared in `pyproject.toml`.

## Install MECHBench

From the repository root, create and activate a Python environment, then install
the project and its core dependencies:

```bash
python -m venv .venv
```

Activate it before continuing:

```bash
# Linux/macOS
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install MECHBench in editable mode:

```bash
python -m pip install -e .
```

## Install optional dependencies

Optional groups are declared under `[project.optional-dependencies]` in
`pyproject.toml`. Install a group by appending its name in square brackets to
the project install command:

```bash
# HDBO and related optimizer dependencies
python -m pip install -e ".[hdbo]"

# ELA (Exploratory Landscape Analysis) dependencies
python -m pip install -e ".[ela]"

# Both optional groups
python -m pip install -e ".[hdbo,ela]"
```

The `hdbo` group includes dependencies for the HDBO and related optimizers,
including BoTorch, GPyTorch, PyTorch, HEBO, and CMA-ES. The `ela` group installs
`pflacco`.

If MECHBench is already installed in the active environment, rerunning the
appropriate command installs the selected optional dependencies as well.

## PyTorch and cluster environments

For GPU runs, install a PyTorch build compatible with the cluster's CUDA
version. On clusters that provide PyTorch as a module, follow the cluster's
environment instructions and confirm the activated Python environment can
import `torch`; a module and a virtual environment do not always share packages.

The repository's `requirements.txt` lists core dependencies. To install optional
groups, use the `pyproject.toml` commands above rather than assuming
`pip install -r requirements.txt` includes them.
