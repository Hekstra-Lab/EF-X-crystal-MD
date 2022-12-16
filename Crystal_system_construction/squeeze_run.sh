#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p gpu
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o squeeze_gpu_%A_%a.out
#SBATCH -e squeeze_gpu_%A_%a.err

source $LOGIN/.bashrc
source $LOGIN/.bash_addons 
conda activate openmm
PENMM_DEFAULT_PLATFORM=CUDA
python squeeze_run.py

