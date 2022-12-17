#!/bin/bash
#SBATCH -t 0-1:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test
#SBATCH -c 8
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o avg_sf_%A_%a.out
#SBATCH -e avg_sf_%A_%a.err

source /n/home11/ziz531/.bashrc
source  $LOGIN/.bash_addons
conda activate openmm
mkdir processed_data/$1
python average_structure_factors.py -i ./$1/data/EF_10MV_cm -N 200
cp ./$1/data/EF_10MV_cm_chainwise_*_subtraj_*_avg.mtz processed_data/$1

