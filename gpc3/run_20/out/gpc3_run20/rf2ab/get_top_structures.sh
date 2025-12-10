grep target_aligned_antibody_rmsd *pdb | sort -nk 3 | head -100 | cut -d : -f 1 | while read -r l
do
pdb_tofasta -multi $l | grep -Ev "PDB|GSDLQ|VRHAKN|VINTT|YILSL|SAYYP|DTLCW" | awk '
{ lines[NR] = $0 }  # Store line content in array 'lines'
END {
    if (NR >= 4) {
        # Concatenate line 2, line 3, the string, line 5, and line 6
        print lines[1] lines[2] "GGGGSGGGGSGGGGS" lines[3] lines[4]
    } else {
        print "Error: File too short."
    }
}' >> top_sequences.txt
done

grep interaction_pae *pdb | sort -nk 3 | head -100 | cut -d : -f 1 | while read -r l
do
pdb_tofasta -multi $l | grep -Ev "PDB|GSDLQ|VRHAKN|VINTT|YILSL|SAYYP|DTLCW" | awk '
{ lines[NR] = $0 }  # Store line content in array 'lines'
END {
    if (NR >= 4) {
        # Concatenate line 2, line 3, the string, line 5, and line 6
        print lines[1] lines[2] "GGGGSGGGGSGGGGS" lines[3] lines[4]
    } else {
        print "Error: File too short."
    }
}' >> top_sequences.txt
done

sort -u top_sequences.txt | nl -w 4 -n rz | sed 's/^/gpc3_rfa_scfv_1_ac_/' >> ../../../orders_instance_gpc3_run_20_dec2025.csv
