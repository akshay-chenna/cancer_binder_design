import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
y = np.array([232,242,243,404,405])
x = u.select_atoms('(chainID C) and (around 4.5 chainID A)')
x = np.unique(x.resids)
fraction = np.mean(np.isin(y,x))
print(f"{f} \t {fraction}")
