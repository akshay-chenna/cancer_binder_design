import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

# Define input and output file paths
PDB_FILE = "6W51_CA_cleaned.pdb"
CIF_FILE = "6W51_CA_cleaned.cif"  # or .mmcif

# --- 1. Read the PDB File ---
# Instantiate the PDBFile class for reading
pdb_file = pdb.PDBFile.read(PDB_FILE)

# Extract the structural data (AtomArray or AtomArrayStack)
# We assume a single model here. Omit 'model=1' to read all models.
atom_array = pdb_file.get_structure(model=1)

# --- 2. Write the Structure to a CIF File ---
# Instantiate the CIFFile class for writing
cif_file = pdbx.CIFFile()

# Write the AtomArray/AtomArrayStack into the CIFFile object
# This handles the conversion of fields from the AtomArray to CIF categories (like _atom_site)
pdbx.set_structure(cif_file, atom_array)

# Write the CIFFile object to the output file on disk
cif_file.write(CIF_FILE)

print(f"Successfully converted '{PDB_FILE}' to '{CIF_FILE}'")
