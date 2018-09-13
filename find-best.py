#!/usr/bin/env python
from __future__ import print_function

import glob
import os
import math
import sys

from collections import namedtuple

# note: # modify below two lines 
# according to your directory structures
globdir = sys.argv[1] 
docked_pdbqt = sorted(glob.glob("%s/COM-*/*docked*.pdbqt" % globdir))

atom = namedtuple("atom", ["serial", "name", "resname", "x", "y", "z"])
pose = namedtuple("pose", ["rec", "model", "score", "atoms", "pdbqt","count"])

def get_rmsd(a, b):
    rmsd = 0.0
    n = 0
    for (atom1, atom2) in zip(a, b):
        assert atom1.name == atom2.name
        rmsd += (atom1.x - atom2.x) ** 2
        rmsd += (atom1.y - atom2.y) ** 2
        rmsd += (atom1.z - atom2.z) ** 2
        n += 1
    return math.sqrt(rmsd / n)

def aggregate_poses (poses, pdbqt):
    with open(pdbqt, "r") as f:
        # get identifiers
	# XXXX-COM-XXXX-docked-#
        prefix = os.path.basename(pdbqt).split(".")[0] 
        rec,comp,ser,docked,repeat = prefix.split("-")
        lig = comp+"-"+ser

        for line in f:
            if line.startswith("MODEL"):
                model = int(line.rstrip().split()[1])
                atoms = []

            if line.startswith("REMARK VINA RESULT:"):
                score = float(line.split()[3])

            if line.startswith("HETATM"):
                (hetatm_, serial, name, resname, chain,
                 comp, y, z, d1_, d2_, q, atomtype) = \
                    line.rstrip().split()
                atoms.append(atom(serial=int(serial),
                                  name=name,
                                  resname=resname,
                                  x=float(x),
                                  y=float(y),
                                  z=float(z),
                                  )
                             )
            if line.startswith("ENDMDL"):
                if not lig in poses:
                    poses[lig] = [[rec,model,score,atoms,pdbqt,1]]
                else:
                    pose_added = False
                    for prev in poses[lig]:
                        if prev[0] != rec :
                            continue
                        rmsd = get_rmsd(prev[3], atoms)
                        if rmsd < 2.0 and abs(prev[2] - score) < 0.5:
                            prev[5] += 1
                            pose_added = True
                            break
                    if not pose_added:
                        poses[lig].append([rec,model,score,atoms,pdbqt,1])
def main() :
    poses = {}
    for pdbqt in docked_pdbqt:
        aggregate_poses(poses, pdbqt)

    for lig in sorted(poses):
        print()
        ranked = sorted(poses[lig], key=lambda x: x[2])
        for i, pose in enumerate(ranked):
            #if i >= 10: continue
            print("{0:10} {1:6.1f} count {2:2d} {3:40} model {4}".format(lig, pose[2], pose[5], pose[4], pose[1]))

if __name__ == "__main__" :
    main()
