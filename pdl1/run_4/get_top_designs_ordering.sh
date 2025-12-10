# Get all designs from protein-hunter, i.e. no filtering
VOLUME=pdl1_run4
DESIGNS=1000
NAME=pdl1_ph_mb_2_ac_

awk -F , '{print $NF}' volume_data/${VOLUME}/job_*/summary_all_runs.csv | sort -u | nl -w 4 -n rz | head -${DESIGNS} | sed "s/^/${NAME}/" >> orders_instance_${VOLUME}_dec2025.csv

