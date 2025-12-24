modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb is gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb upto A1-291 and B1-108
pdb_selchain -A modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb | pdb_rplchain -A:B > chA2B.pdb
pdb_selchain -B modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb | pdb_rplchain -B:A > chB2A.pdb
vi chB2A.pdb 
pdb_merge chB2A.pdb chA2B.pdb > reodered_modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb 
