"""

Script based on Jack Greisman's code, written by Ziyuan on Feb 7, 2022.

"""
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import mdtraj
import os
import sys, getopt


try:
    opts, args = getopt.getopt(sys.argv[1:], "Dvd:i:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

dirname = ""
verbose = False
dryrun = False
inputname = "iter"
outputname = "squeeze_analysis_plot.png"

for o, a in opts:
    if o == "-d":
        dirname = a
    elif o == "-v":
        verbose = True
    elif o == "-D":
        dryrun = True
    elif o == "-i":
        inputname = a
    elif o == "-o":
        outputname = a

sns.set_context("notebook", font_scale=1.3)

if dirname != "":
    try: 
        os.mkdir(dirname)
    except:
        pass # ignore err due to existing dir 
    os.chdir(dirname)

h5trajs = sorted(glob.glob(inputname+"*.h5"))

sns.set_palette(sns.cubehelix_palette(len(h5trajs), start=.5, rot=-.75))
plt.figure(figsize=(9, 6))

if dryrun:
    sys.exit(1)

for h5traj in h5trajs:
   traj = mdtraj.load(h5traj)
   n = int(len(traj.topology.select("water")) / 3)
   plt.plot(traj.time, traj.unitcell_volumes, label=f"{n} waters")

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Time (ps)")
plt.ylabel(r"Box Volume (nm$^3$)")
plt.xlim(0, traj.time[-1])
plt.tight_layout()
plt.savefig(outputname)

if verbose:
    print("Figured successfully saved!")