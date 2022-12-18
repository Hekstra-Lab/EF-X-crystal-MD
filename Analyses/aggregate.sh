#!/bin/bash
#SBATCH -t 0-1:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test
#SBATCH -c 16
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o aggregate_%A_%a.out
#SBATCH -e aggregate_%A_%a.err

source /n/home11/ziz531/.bashrc
source  $LOGIN/.bash_addons
conda activate openmm
mkdir processed_data/aggregate
for i in 30 60 90 120 150;
do
mkdir processed_data/aggregate/$i
done
python aggregate.py
