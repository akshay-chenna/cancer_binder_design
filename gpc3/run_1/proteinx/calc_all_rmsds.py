import numpy as np
import MDAnalysis as mda
import sys
from MDAnalysis.analysis.rms import RMSD

## TARGET
def target_bb_rmsd(u,v):
        u = mda.Universe(u)
        u = u.select_atoms('(resid 1A:291A or resid 325A:376A) and (name CA or name C or name N)')
        v = mda.Universe(v)
        v = v.select_atoms('(resid 1A:291A or resid 325A:376A) and (name CA or name C or name N)')

        a = RMSD(v,u).run().rmsd[:,2]
        return round(float(a),2)


## BINDER
def binder_bb_rmsd(u,v):
        u = mda.Universe(u)
        u = u.select_atoms('chainID B and (name CA or name C or name N)')
        v = mda.Universe(v)
        v = v.select_atoms('chainID B and (name CA or name C or name N)')

        a = RMSD(v,u).run().rmsd[:,2]
        return round(float(a),2)

## TARGET ALIGNED BINDER RMSD
def targetaligned_binder_bb_rmsd(u,v):
        u = mda.Universe(u)
        v = mda.Universe(v)
        a = RMSD(v,u,select='(resid 1A:291A or resid 325A:376A) and (name CA or name C or name N)', groupselections=['chainID B and (name CA or name C or name N)']).run().rmsd[:,3]
        return round(float(a),2)
        
t = target_bb_rmsd(sys.argv[1], sys.argv[2])
b = binder_bb_rmsd(sys.argv[1], sys.argv[2])
tb = targetaligned_binder_bb_rmsd(sys.argv[1], sys.argv[2])

a = sys.argv[1]
name = (a.split('/', 1)[1] if '/' in a else a).split('.', 1)[0] if '.' in (a.split('/', 1)[1] if '/' in a else a) else ''

print(name,t,b,tb,sep='\t')
