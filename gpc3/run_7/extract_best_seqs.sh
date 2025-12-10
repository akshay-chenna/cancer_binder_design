awk -F , '{print $3}' run_?/accepted.csv | sort -u | nl -w 4 -n rz | sed 's/^/gpc3_mber_nb_h1_ac_/'  > gpc3_mber_nb_h1.csv
