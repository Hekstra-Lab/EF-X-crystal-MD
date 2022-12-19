import reciprocalspaceship as rs
from mdtools.utils import *
import argparse
from copy import deepcopy

# parameters
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input file for the crystal system")
# parser.add_argument("-o", "--output", type=str, help="Prefix for the output trajectory and state files")
parser.add_argument("-n", "--n_chains", type=int, help="Number of protein chains in the system", default=4)
parser.add_argument("-N0", "--n_epochs_offset", type=int, default=0)
parser.add_argument("-N", "--n_epochs", type=int, help="Number of epochs to average over")
parser.add_argument("-dN", "--step", type=int, help="Number of steps between consecutive saves (useful for convergence test), set to 0 to save just the final result", default=0)
args = parser.parse_args()

dataset = rs.read_mtz(args.input+f'_epoch_0_chainwise_pos_subtraj_0_avg.mtz')
n_reflections = dataset.shape[0]

# need to do this so for phenix refinemenet, a bit of a hack...
def reformat_dataset(dataset, sigma=0.05):
    dataset['FOBS'] = dataset['FMODEL']
    dataset['SIGFOBS'] = (sigma*dataset['FOBS']).astype("Q")
    dataset = dataset[['FOBS', 'SIGFOBS']]
    
def to_dataset(dataset, complex_reflections):
    ds = deepcopy(dataset)
    ds[:] = cartesian_arr_to_polar(complex_reflections)
    ds.infer_mtz_dtypes(inplace = True)
    reformat_dataset(ds)
    return ds

for phase in ['pos', 'neg', 'zero']:
    for k in range(args.n_chains):
        complex_reflections = np.zeros(n_reflections, dtype='complex128')

        for epoch in range(args.n_epochs_offset, args.n_epochs):
            dataset = rs.read_mtz(f'{args.input}_epoch_{epoch}_chainwise_{phase}_subtraj_{k}_avg.mtz')
            complex_reflections = complex_reflections * (1 - 1/(epoch + 1)) + mtz_to_cartesian_arr(dataset) / (epoch + 1)
            if args.step > 0 and (epoch+1) % args.step == 0:
                intermediate = to_dataset(dataset, complex_reflections)
                dataset.write_mtz(f'{args.input}_partial_{epoch+1}_chainwise_{phase}_subtraj_{k}_avg.mtz')
        dataset = to_dataset(dataset, complex_reflections)
        dataset.write_mtz(f'{args.input}_chainwise_{phase}_subtraj_{k}_avg.mtz')
                

for epoch in range(args.n_epochs_offset, args.n_epochs):     
    if args.step > 0 and (epoch+1) % args.step == 0:
        compute_all_difference_maps(f'{args.input}_partial_{epoch+1}')
# compute_all_difference_maps(args.input)
