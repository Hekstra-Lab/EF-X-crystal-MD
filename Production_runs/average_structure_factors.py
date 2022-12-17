import reciprocalspaceship as rs
from mdtools.utils import *
import argparse

# parameters
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input file for the crystal system")
# parser.add_argument("-o", "--output", type=str, help="Prefix for the output trajectory and state files")
parser.add_argument("-n", "--n_chains", type=int, help="Number of protein chains in the system", default=4)
parser.add_argument("-N0", "--n_epochs_offset", type=int, default=0)
parser.add_argument("-N", "--n_epochs", type=int, help="Number of epochs to average over")
args = parser.parse_args()

dataset = rs.read_mtz(args.input+f'_epoch_0_chainwise_pos_subtraj_0_avg.mtz')
n_reflections = dataset.shape[0]

for phase in ['pos', 'neg', 'zero']:
    for k in range(args.n_chains):
        complex_reflections = np.zeros(n_reflections, dtype='complex128')

        for epoch in range(args.n_epochs_offset, args.n_epochs):
            dataset = rs.read_mtz(args.input+f'_epoch_{epoch}_chainwise_{phase}_subtraj_{k}_avg.mtz')
            complex_reflections = complex_reflections * (1 - 1/(epoch + 1)) + np.array([amp*np.exp(np.pi*phase/180 * 1j) for [amp, phase] in dataset.to_numpy()]) / (epoch + 1)


        dataset[:] = np.stack([np.abs(complex_reflections), np.angle(complex_reflections) / np.pi * 180]).T
        dataset.infer_mtz_dtypes(inplace = True)
        dataset.write_mtz(args.input+f'_chainwise_{phase}_subtraj_{k}_avg.mtz')

compute_all_difference_maps(args.input)
