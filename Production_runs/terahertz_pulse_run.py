# import stuff
from openmm.app import PDBFile, ForceField
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from simtk.unit import *
import mdtraj
import mdtools
import pickle
import argparse
import os

def getFieldStrength(e):
    """
    Convert from given unit of electric field strength to a value
    in OpenMM's standard unit
    """
    def convert(v):
        return (v*AVOGADRO_CONSTANT_NA).value_in_unit(kilojoule_per_mole / elementary_charge / nanometer)

    if isinstance(e, list):
        return [convert(e_c) for e_c in e]
    else:
        return convert(e)

    
parser = argparse.ArgumentParser()
parser.add_argument("-E", "--E", type=str, help="Field strength of the THz pulse (MV/cm)", default=10e8)
parser.add_argument("-t0", "--t0", type=int, help="Initial duration of NVT equilibration (ns)", default=10)
parser.add_argument("-t1", "--t1", type=int, help="Number of pulse cycles for the production run")
parser.add_argument("-o", "--output", type=str, help="Prefix for the output trajectory and state files")
parser.add_argument("-r", "--replicate", type=int, help="Index number of the current replicate")
args = parser.parse_args()

def get_file_path(file):
    root = f'/n/holyscratch01/hekstra_lab/ziyuan/EF-X-crystal-MD/Production_runs/{args.replicate}'
    for (root, dirs, files) in os.walk(root):
        if file in files:
            return os.path.join(root, file)
        
# initialize forcefield
# need to have openmmforcefields installed first
forcefield = ForceField('amber/ff14SB.xml',
                        'amber/tip3p_standard.xml',
                        'amber/tip3p_HFE_multivalent.xml',
                        get_file_path('cro.xml'))
acetate = Molecule.from_smiles('CC(=O)[O-]')
gaff = GAFFTemplateGenerator(molecules=acetate)
forcefield.registerTemplateGenerator(gaff.generator)

efx = getFieldStrength(args.E * volts/meters)
print("Field strength is", efx)

crystal = PDBFile(get_file_path("squeezed.pdb"))
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
                         saveTrajectory=True, saveStateData=True,
                         trajInterval=50, stateDataInterval=50,
                         dt=0.002*picoseconds, efx=True) # record every 0.1ps
#mdsystem.simulate(10*nanoseconds) # a pulse sequence lasts 10ns
for i in range(args.t1):
    # Up
    mdsystem.simulation.context.setParameter('Ex', efx)
    mdsystem.simulate(1*picoseconds) # 1ps

    # Down
    mdsystem.simulation.context.setParameter('Ex', -1*efx)
    mdsystem.simulate(1*picoseconds) # 1ps

    # Off
    mdsystem.simulation.context.setParameter('Ex', 0.0)
    mdsystem.simulate(8*picoseconds) # 8ps
print("Finished terahertz pulse production run", flush=True)

