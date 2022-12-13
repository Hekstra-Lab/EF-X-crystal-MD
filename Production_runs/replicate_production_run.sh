cd $1
mkfifo fifo_pipe
sbatch ../thz_continuous.sh
