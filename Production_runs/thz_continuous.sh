#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p gpu
#SBATCH -c 32
#SBATCH --mem=128G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o thz_cont_gpu_%A_%a.out
#SBATCH -e thz_cont_gpu_%A_%a.err

source /n/home11/ziz531/.bashrc
source  $LOGIN/.bash_addons
conda activate openmm
OPENMM_DEFAULT_PLATFORM=CUDA
python /n/holyscratch01/hekstra_lab/ziyuan/EF-X-crystal-MD/Production_runs/thz_continuous_simulation.py -n 4 -o EF_10MV_cm -t0 0 -t1 200 -t2 100 -r $1&
python /n/holyscratch01/hekstra_lab/ziyuan/EF-X-crystal-MD/Production_runs/thz_continuous_postproc.py -n 4 -o EF_10MV_cm -t0 0 -t1 200 -t2 100 -r $1

