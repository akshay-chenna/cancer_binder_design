import MDAnalysis as mda
import numpy as np
u = mda.Universe('test.pdb','test.pdb')
x = u.select_atoms('(chainID B or chainID C) and (around 6 chainID A)')
print(len(np.unique(x.resids)))
print(np.unique(x.resids))
