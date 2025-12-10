import MDAnalysis as mda
import numpy as np
u = mda.Universe('model_0_neo.pdb','model_0_neo.pdb')
x = u.select_atoms('chainID B and (around 6 (chainID A and resid 4:8))')
print(len(np.unique(x.resids)))
print(np.unique(x.resids))
