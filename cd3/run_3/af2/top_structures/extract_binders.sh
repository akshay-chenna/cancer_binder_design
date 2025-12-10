cat top_visually_selected.txt | while read -r l ; do pdb_selchain -A ${l}.pdb | pdb_keepcoord >> chA_${l}.pdb ; done
