#!/usr/bin/python
import __main__
__main__.pymol_argv = ['pymol','-qc']

import pymol
import sys
import os

pymol.finish_launching()
pdbqt= sys.argv[1] # ex. rec-COMP-1234-docked-1.pdbqt
tmp= os.path.basename(pdbqt).split(".")[0].split("-")
rec= tmp[0]
lig= tmp[1]+"-"+tmp[2]
recpdb = "./rec/%s.pdb" % rec
outdir= os.path.dirname(pdbqt)
outpdb= os.path.join(outdir,"%s-%s.pdb" % (rec,lig))
pymol.cmd.load(recpdb,rec)
pymol.cmd.load(pdbqt,lig)
pymol.cmd.save(outpdb)
pymol.cmd.reinitialize()
pymol.cmd.quit()
