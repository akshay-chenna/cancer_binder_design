import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
x = u.select_atoms('(chainID A) and (around 4.5 chainID B)')
x = np.unique(x.resids)
print(x)
