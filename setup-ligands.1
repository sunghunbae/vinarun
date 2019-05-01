#!/usr/bin/env /home/shbae/anaconda3/envs/work/bin/python

import sys
import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# modify below two lines according to your directory structures
pythonsh_path = "/home/shbae/bin/mgltools-1.5.7rc1/bin/pythonsh"
prepare_ligand4_path = "/home/shbae/bin/autodock4/prepare_ligand4.py"

def generate_3d(smi, name, output_dir, ff='MMFF', save_pdbqt=True):
    m = Chem.MolFromSmiles(smi)
    m.SetProp("_Name", name)

    # calculate energy minimized 3D coordinates
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m,randomSeed=0xf00d,
	useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    if ff == 'MMFF':
        AllChem.MMFFOptimizeMolecule(m,maxIters=2000)
    elif ff == 'UFF':
        AllChem.UFFOptimizeMolecule(m,maxIters=2000)

    # calculate radius of gyration
    # optimal box size for Autodock Vina is 2.9 times of Rg
    # Feinstein and Brylinski (2015), J. Cheminfo. 7:18
    # DOI: 10.1186/s13321-015-0067-5

    Rg = rdMolDescriptors.CalcRadiusOfGyration(m)
    boxsize = Rg * 2.857

    if save_pdbqt :
        pdb_file = os.path.join(output_dir,name+".pdb")
        pdbqt_file = os.path.join(output_dir,name+".pdbqt")
        Chem.MolToPDBFile(m, pdb_file, flavor=4)
        cmd  = [ pythonsh_path,
                 prepare_ligand4_path,
                 "-l", pdb_file ,
                 "-o", pdbqt_file, ]
        x=subprocess.Popen(cmd,stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        out = [x.decode("utf-8") for x in x.stdout.readlines()]
        err = [x.decode("utf-8") for x in x.stderr.readlines()]
        if len(out) > 1 or len(err) > 0 :
            print("".join(out))
            print("".join(err))
        else:
            with open(pdbqt_file,"a") as f:
                f.write("REMARK  optimal box size %.2f (Rg= %.3f)" % (boxsize,Rg))
                f.write(" Feinstein and Brylinski (2015)")
            print("%-10s Rg= %.3f size= %.2f  --->  %s" %
                  (name, Rg, boxsize, pdbqt_file))
    else:
        print("%-10s Rg= %.3f size= %.2f" % (name, Rg, boxsize))

def main() :
    for smifilename in sys.argv[1:]:
        try:
            f = open(smifilename, "r")
        except:
            print("cannot open SMILES file", smifilename)
            continue
        output_dir = os.path.dirname(smifilename)
        for line in f:
            try:
                smi, name = line.strip().split()
            except:
                smi = line.strip()
                name = os.path.basename(smifilename).split(".")[0]
            generate_3d(smi, name, output_dir)

if __name__ == '__main__':
    main()
