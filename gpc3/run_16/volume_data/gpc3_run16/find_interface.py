import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
y = np.array([67,71,147,153,156,275])
x = u.select_atoms('(chainID C) and (around 4.5 chainID A)')
x = np.unique(x.resids)
fraction = np.mean(np.isin(y,x))
print(f"{f} \t {fraction}")
