from openmm.app import PDBFile, ForceField
from openmmforcefields.generators import GAFFTemplateGenerator
from simtk.unit import *
from tqdm import tqdm
import mdtraj
from mdtraj.reporters import HDF5Reporter
import mdtools
from mdtools.utils import *
import pickle
import argparse
import subprocess
import os

#
fifo_name = 'fifo_pipe'

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
args = parser.parse_args()

def get_file_path(file):
    root = '/n/hekstra_lab/people/ziyuan'
    for (root, dirs, files) in os.walk(root):
        if file in files:
            return os.path.join(root, file)

# initialize forcefield
# need to have openmmforcefields installed first
forcefield = ForceField('amber/ff14SB.xml',
                        'amber/tip3p_standard.xml',
                        'amber/tip3p_HFE_multivalent.xml',
                        get_file_path('cro.xml'))
efx = getFieldStrength(args.E * volts/meters)
print("Field strength is", efx)

crystal = PDBFile(get_file_path(args.input))
mdsystem = mdtools.LatticeMDSystem(crystal.topology,
                                   crystal.positions,
                                   forcefield, "P 21 21 21")
print("Started simulation.", flush=True)
mdsystem.buildSimulation(ensemble="NVT", posre=True,
                         saveTrajectory=False, saveStateData=False,
                         posre_sel="not water and not (element Na or element Ca) and not element H",
                         dt=0.002*picoseconds)
mdsystem.equilibrate(args.t0*nanoseconds, posre=True)
print("Finished NVT equilibration.", flush=True)
mdsystem.buildSimulation(ensemble="NVT",  filePrefix=args.output,
                         saveTrajectory=True, saveStateData=False,
                         trajInterval=50, stateDataInterval=50,
                         dt=0.002*picoseconds, efx=True) # record every 0.1ps

mdsystem.simulation.reporters[0].close()

for epoch in tqdm(range(args.epoch_offset, args.t1)):
    mdsystem.simulation.reporters[0] = HDF5Reporter(f"{args.output}.h5", 50)
    for j in tqdm(range(args.n_pulses)):
        # Up
        mdsystem.simulation.context.setParameter('Ex', efx)
        mdsystem.simulate(1*picoseconds) # 1ps

        # Down
        mdsystem.simulation.context.setParameter('Ex', -1*efx)
        mdsystem.simulate(1*picoseconds) # 1ps

        # Off
        mdsystem.simulation.context.setParameter('Ex', 0.0)
        mdsystem.simulate(8*picoseconds) # 8ps

    # post-process this segment
    mdsystem.simulation.reporters[0].close()

    #
    with open(fifo_name, 'w') as f:
        f.write(f'{epoch}\n')
    with open(fifo_name, 'r') as f:
        assert(f.readline().strip()=="pass")

print("Finished terahertz pulse production run", flush=True)


