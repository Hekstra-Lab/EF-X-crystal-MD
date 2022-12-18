#!/bin/bash
#SBATCH -t 0-1:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p test
#SBATCH -c 16
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o avg_sf_%A_%a.out
#SBATCH -e avg_sf_%A_%a.err

source /n/home11/ziz531/.bashrc
source  $LOGIN/.bash_addons
conda activate openmm
#mkdir -p processed_data/$1/full
for EPOCHS in 30 60 90 120 150 
do 
mkdir -p processed_data/$1/$EPOCHS 
done
echo "Finished creating dirs"

python average_structure_factors.py -i ./$1/data/EF_10MV_cm -N 150 -dN 30

#cp ./$1/data/EF_10MV_cm_chainwise_*_subtraj_*_avg.mtz processed_data/$1/full
for EPOCHS in 30 60 90 120 150 
do 
cp ./$1/data/EF_10MV_cm_partial_${EPOCHS}_chainwise_*_subtraj_*_avg.mtz processed_data/$1/$EPOCHS
done

#cp ./$1/data/EF_10MV_cm_diff_*.mtz processed_data/$1/full
for EPOCHS in 30 60 90 120 150 
do 
cp ./$1/data/EF_10MV_cm_partial_${EPOCHS}_diff_*.mtz processed_data/$1/$EPOCHS
done