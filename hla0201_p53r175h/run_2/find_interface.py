import MDAnalysis as mda
import numpy as np
import sys

f = sys.argv[1]

u = mda.Universe(f,f)
y = np.array([73,76,77,80,81,84,143,146,147,150,152,237,238,239])
z = np.array([237,238,239])
x = u.select_atoms('(chainID A) and (around 4.5 chainID B)')
x = np.unique(x.resids)
fraction = np.mean(np.isin(y,x))
fraction_peptide = np.mean(np.isin(z,x))
print(f"{f} \t {fraction} \t {fraction_peptide}")
