pdb_fetch 6w51 > 6W51.pdb
pdb_fetch 6w51 | pdb_selchain -A,C | pdb_delhetatm | pdb_keepcoord | pdb_selres -1:180 | pdb_reatom | pdb_rplchain -C:B | pdb_reres | pdb_tidy > 6W51_AC_cleaned.pdb
pdb_selchain -B 6W51_AC_cleaned.pdb | pdb_rplchain -B:A > 6W51_C_cleaned.pdb
pdb_selchain -A 6W51_AC_cleaned.pdb | pdb_rplchain -A:B > 6W51_A_cleaned.pdb
pdb_merge 6W51_C_cleaned.pdb 6W51_A_cleaned.pdb | pdb_reres | pdb_reatom | pdb_tidy > 6W51_CA_cleaned.pdb
