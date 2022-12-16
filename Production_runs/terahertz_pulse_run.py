# import stuff
from openmm.app import PDBFile, ForceField
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from simtk.unit import *
import mdtraj
import mdtools
import pickle

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


# initialize forcefield
# need to have openmmforcefields installed first
forcefield = ForceField('amber/ff14SB.xml',
                        'amber/tip3p_standard.xml',
                        'amber/tip3p_HFE_multivalent.xml',
                        'cro.xml')
efx = getFieldStrength(10e8 * volts/meters)
print("Field strength is", efx)

crystal = PDBFile("squeezed.pdb")
mdsystem = mdtools.LatticeMDSystem(crystal.topology,
                                   crystal.positions,
                                   forcefield, "P 21 21 21")
print("Started simulation.", flush=True)
mdsystem.buildSimulation(ensemble="NVT", posre=True,
                         saveTrajectory=False, saveStateData=False,
                         posre_sel="not water and not (element Na or element Ca) and not element H",
                         dt=0.002*picoseconds)
mdsystem.equilibrate(10*nanoseconds, posre=True)
print("Finished NVT equilibration.", flush=True)
mdsystem.buildSimulation(ensemble="NVT",  filePrefix=f"EF_10MV_cm_2",
                         saveTrajectory=True, saveStateData=True,
                         trajInterval=50, stateDataInterval=50,
                         dt=0.002*picoseconds, efx=True) # record every 0.1ps
#mdsystem.simulate(10*nanoseconds) # a pulse sequence lasts 10ns
for i in range(5000):
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

