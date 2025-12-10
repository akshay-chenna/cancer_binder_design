import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
y = np.array([287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303])
x = u.select_atoms('(chainID B) and (around 4.5 chainID A)')
x = np.unique(x.resids)
fraction = np.mean(np.isin(y,x))
print(f"{f} \t {fraction}")
