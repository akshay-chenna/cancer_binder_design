pdb_fetch 6w51 | pdb_selchain -A,C | pdb_delhetatm | pdb_keepcoord | pdb_selres -1:180 | pdb_rplchain -C:B | pdb_reatom | pdb_tidy > 6W51_AC_cleaned.pdb
