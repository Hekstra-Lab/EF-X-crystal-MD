# import stuff
from openmm.app import PDBFile, ForceField
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from simtk.unit import *
import mdtraj
import mdtools
import pickle 


# initialize forcefield
# need to have openmmforcefields installed first
forcefield = ForceField('amber/ff14SB.xml', 
                        'amber/tip3p_standard.xml',
                        'amber/tip3p_HFE_multivalent.xml',
                        'cro.xml')
acetate = Molecule.from_smiles('CC(=O)[O-]')
gaff = GAFFTemplateGenerator(molecules=acetate)
forcefield.registerTemplateGenerator(gaff.generator)

# load the crystal system
crystal = PDBFile("neutralized.pdb")
mdsystem = mdtools.LatticeMDSystem(crystal.topology,
                                   crystal.positions, 
                                   forcefield, "P 21 21 21")

mdsystem.buildSimulation()
print("System initialized for squeezing. Currently running on platform ", mdsystem.simulation.context.getPlatform().getName())
mdsystem.squeeze(tolerance=0.0005, maxIterations=20, dt=0.002*picoseconds, initial_water_perturb=1000, dn=500)
mdsystem.save("squeezed.pdb")
