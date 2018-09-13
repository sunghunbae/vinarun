#!/usr/bin/env /home/shbae/anaconda3/envs/work/bin/python

import glob
import sys
import os
import subprocess

# modify below three lines 
# according to your directory structures
vina_path = "/home/shbae/bin/vina-1.1.2/vina"
receptors = glob.glob("./rec/xyz.pdbqt")
ligands = glob.glob("./lig/com-1234.pdbqt")

def dock(repeat=1, use_Rg_boxsize=True):
    for lig_pdbqt in ligands :
        lig = os.path.basename(lig_pdbqt).split(".")[0] # COM-XXX.pdbqt
        ligdir = os.path.dirname(lig_pdbqt)
        output_dir = os.path.join(ligdir,lig)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # read ligand pdbqt and get optimal box size
        boxsize = None
        with open(lig_pdbqt,"r") as f:
            for line in f:
                if line.startswith("REMARK  optimal box size"):
                    boxsize = line.strip().split()[4] # string

        for rec_pdbqt in receptors:
            rec= os.path.basename(rec_pdbqt)[:4]
            recdir= os.path.dirname(rec_pdbqt)
            config= os.path.join(recdir,rec+".txt")
            with open(config,"r") as f:
                box_setting = {}
                for line in f:
                    c = line.strip().split()
                    box_setting[c[0]] = c[2] # string

            docked = "%s/%s-%s-docked" % (output_dir, rec, lig)

            with open(docked+".txt","w") as f:
                for k in ["size_x","size_y","size_z"]:
                    if use_Rg_boxsize and boxsize:
                        f.write(k + " = " + boxsize +"\n")
                    else:
                        f.write(k + " = " + box_setting[k] + "\n")
                for k in ["center_x","center_y","center_z"]:
                    f.write(k + " = " + box_setting[k] + "\n")

            for n in range(0,repeat):

                cmd = [vina_path,
                       "--cpu","4",
                       "--receptor", rec_pdbqt,
                       "--ligand", lig_pdbqt,
                       "--config", docked+".txt",
                       "--out", docked+"-"+str(n+1)+".pdbqt",
                       "--log", docked+"-"+str(n+1)+".log",
                       ]

                x = subprocess.Popen(cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
                out = [x.decode("utf-8") for x in x.stdout.readlines()]
                err = [x.decode("utf-8") for x in x.stderr.readlines()]
                if len(out) > 1 or len(err) > 0:
                    print("".join(out))
                    print("".join(err))

if __name__ == '__main__':
    dock(repeat=10, use_Rg_boxsize=True)
