import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
y = np.array([31,33,70,71,75,76,81,86,87])
x = u.select_atoms('(chainID C) and (around 4.5 chainID A)')
x = np.unique(x.resids)
fraction = np.mean(np.isin(y,x))
print(f"{f} \t {fraction}")
