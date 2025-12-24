#modified_gpc3_out_1_l65_s619314_mpnn1_model2_relaxed.pdb is gpc3_out_1_l65_s599292_mpnn3_model2_relaxed.pdb upto A1-291 and B1-65
pdb_selchain -A modified_gpc3_out_1_l65_s619314_mpnn1_model2_relaxed.pdb | pdb_rplchain -A:B > chA2B.pdb
pdb_selchain -B modified_gpc3_out_1_l65_s619314_mpnn1_model2_relaxed.pdb | pdb_rplchain -B:A > chB2A.pdb
pdb_merge chB2A.pdb chA2B.pdb > reordered_modified_gpc3_out_1_l65_s619314_mpnn1_model2_relaxed.pdb 
