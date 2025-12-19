import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
y = np.array([46,85,87,92,93])
x = u.select_atoms('(chainID B) and (around 4.5 chainID A)')
x = np.unique(x.resids)
fraction = np.mean(np.isin(y,x))
print(f"{f} \t {fraction}")
