from simtk.unit import *
from tqdm import tqdm
import mdtraj
import mdtools
from mdtools.utils import *
import argparse
import subprocess
import os
import multiprocessing as mp
from itertools import product

fifo_name = 'fifo_pipe'
cpu_count = mp.cpu_count()
print(f'Running on {cpu_count} cpus')

# parameters
parser = argparse.ArgumentParser()
parser.add_argument("-E", "--E", type=str, help="Field strength of the THz pulse (MV/cm)", default=10e8)
parser.add_argument("-t0", "--t0", type=int, help="Initial duration of NVT equilibration (ns)", default=10)
parser.add_argument("-t1", "--t1", type=int, help="Number of pulse cycles for the production run")
parser.add_argument("-t2", "--n_pulses", type=int, help="Number of pulses per cycle", default=100)
parser.add_argument("-i", "--input", type=str, help="Input file for the crystal system", default="squeezed.pdb")
parser.add_argument("-o", "--output", type=str, help="Prefix for the output trajectory and state files")
parser.add_argument("-n", "--n_chains", type=int, help="Number of protein chains in the system", default=4)
parser.add_argument("-e", "--epoch_offset", type=int, help="The epoch to start from", default=0)
parser.add_argument("-r", "--replicate", type=int, help="Index number of the current replicate")
args = parser.parse_args()

def get_file_path(file):
	root = f'/n/holyscratch01/hekstra_lab/ziyuan/EF-X-crystal-MD/Production_runs/{args.replicate}'
	for (root, dirs, files) in os.walk(root):
		if file in files:
			return os.path.join(root, file)

def aux(epoch, phase, k):
    fname = args.output+f'_epoch_{epoch}_chainwise_{phase}_subtraj_{k}'
    save_snapshots_from_traj(mdtraj.load(f'{fname}.h5'), output_name=fname, frame_offset=0, d_frame=1) # should give 100 frames
    batch_annotate_spacegroup(fname, args.n_pulses, "P 21 21 21")
    batch_fmodel(fname, max_frame=args.n_pulses, resolution=1.5,
                 phenix_command='source /n/hekstra_lab/people/ziyuan/egfp/phenix_env.sh; phenix.fmodel')
    average_structure_factors(fname, max_frame=args.n_pulses)

asu_ref = mdtraj.load(get_file_path('asu_ref.h5'))
# unitcell_ref = mdtraj.load(get_file_path('unitcell_ref.h5'))
atom_selection = np.load(get_file_path('atoms_for_alignment.npy'))
print("ASU reference retrieved from", get_file_path('asu_ref.h5'))
print("Atom selection retrieved from", get_file_path('atoms_for_alignment.npy'))
print("Loaded auxiliary data files.")

while True:
    with open(fifo_name, 'r') as f:
        epoch = int(f.readline())
    # remove solvent, unwrap, align
    traj = mdtraj.load(f'{args.output}.h5')
    traj.remove_solvent()
    unwrap_time_axis(traj)
    def aux_(tupl):
        offset, phase = tupl
        fname=args.output+f'_epoch_{epoch}_chainwise_{phase}'
        align_and_split_by_chain(traj[offset::100], fname,
                                unitcell_ref=None, asu_ref=asu_ref,
                                sg=19, chainwise_alignment=True,
                                atom_selection=atom_selection)

    with mp.Pool(processes=cpu_count) as pool:
        pool.map(aux_, [(9, 'pos'), (19, 'neg'), (99, 'zero')])

    print("Aligned and split into subtrajs")
   
    with open(fifo_name, 'w') as f:
        f.write("pass\n")
    # convert to snapshots, calculate structural factors, and average
    def aux_(tupl):
        phase, k = tupl
        aux(epoch, phase, k)

    with mp.Pool(processes=cpu_count) as pool:
        pool.map(aux_, product(['pos', 'neg', 'zero'], np.arange(args.n_chains)))
    print("Computed average reflection for pos, neg, and zero parts")

    # cleanup
    subprocess.run('mv *avg* ./data', shell=True)
    subprocess.run('rm *subtraj*', shell=True)
    # subprocess.run(f'rm {args.output}.h5', shell=True) --> we write on it over and over again!
    print("Cleaned up intermediate files")

print("Finished terahertz pulse production run", flush=True)


