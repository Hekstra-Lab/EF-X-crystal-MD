cd $1
mkfifo fifo_pipe
mkdir data
ln -s ../asu_ref.h5 .
ln -s ../atoms_for_alignment.npy .
sbatch ../thz_continuous.sh $1
