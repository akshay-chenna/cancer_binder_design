import biotite.structure.io as strucio
from atomworks.io.utils.io_utils import load_any
import sys
from pathlib import Path

input = sys.argv[1]
original_path = Path(input)
new_path = original_path.parent / (original_path.name.split('.',1)[0] + ".pdb")

atom_array = load_any(original_path)
strucio.save_structure(new_path, atom_array)
