mkdir $1 
cd $1
ln -s ../../CRO_parametrization/cro.xml .
ln -s ../../Crystal_system_construction/neutralized.pdb .
ln -s ../../Crystal_system_construction/squeeze_run.py .
ln -s ../../Crystal_system_construction/squeeze_run.sh
sbatch squeeze_run.sh