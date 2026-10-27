#!/bin/bash
#SBATCH --job-name=experiments_baxus_botorch
#SBATCH --mail-user=olarterodriguezi@vuw.leidenuniv.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/olarterodriguezi/logs/experiments_baxus_botorch_%A_%a.out
#SBATCH --error=/home/olarterodriguezi/logs/experiments_baxus_botorch_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=17
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2048M
#SBATCH --partition=cpu-zen4
#SBATCH --time=1-23:59:59
#SBATCH --array=26


#--old SLURM --- --partition=cpu-medium
#--new SLURM --- --partition=cpu-zen4
#--new SLURM - heavy-duty --- --partition=cpu-skylake

# Purge modules
module purge
# load modules (assuming you start from the default environment)
# we explicitly call the modules to improve reproducibility
# in case the default settings change

#FOSS
module load ALICE/default
module load slurm

# ACTIVATE OPENMPI
module load OpenMPI/4.1.5-GCC-12.3.0

#Activate libGLU (This is a helper if ALICE/default does not load it)
module load libGLU/9.0.3-GCCcore-12.3.0

#PYTHON
module load Python/3.11.3-GCCcore-12.3.0

#PYTORCH
module load PyTorch/2.1.2-foss-2023a-CUDA-12.1.1

# Activate the Python virtual environment
source $HOME/BO_torch/bo-torch-run/bin/activate

# Change source directory to the home directory
cd $HOME/MECHBench

# Add executable permissions to the open_radioss binaries
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/starter_linux64_gf
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/engine_linux64_gf
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/engine_linux64_gf_ompi
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/th_to_csv_linux64_gf
chmod +x /home/olarterodriguezi/OpenRadioss_linux64/OpenRadioss/exec/anim_to_vtk_linux64_gf

# SPECIFY PROPERTIES
#ALGORITHM="HEBO"
#ALGORITHM="pyCMA"
#ALGORITHM="turbo1"
#ALGORITHM="BO_botorch"
#OBJECTIVE_FUNCTION="penalized_mass"
#OBJECTIVE_FUNCTION="penalized_sea"
#OBJECTIVE_FUNCTION="load_uniformity"
#ALGORITHM="BAxUS_botorch"

#SPECIFY PROPERTIES
ALGORITHM="BAxUS_botorch"
DIMENSION=20
PROBLEM_TYPE=2
RUN_ID=$SLURM_ARRAY_TASK_ID
MAX_ITER=600
NP=16
OBJECTIVE_FUNCTION="penalized_mass"
H_LEVEL=1



echo "Running topology optimization with Optimizer: $ALGORITHM, Dimension: $DIMENSION, Problem Type: $PROBLEM_TYPE, Run ID: $RUN_ID, Max Iter: $MAX_ITER, NP: $NP, H_LEVEL: $H_LEVEL, Objective Function: $OBJECTIVE_FUNCTION"
# Replace the following line with your actual command
python3 $HOME/MECHBench/main_experiments.py --problem_type $PROBLEM_TYPE --dimension $DIMENSION --optimizer $ALGORITHM --run_id $RUN_ID --max_iter $MAX_ITER --np $NP --objective_function $OBJECTIVE_FUNCTION --h_level $H_LEVEL --nt 1