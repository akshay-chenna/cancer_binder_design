# Get all designs from protein-hunter, i.e. no filtering
VOLUME=trca_run7
DESIGNS=2000
NAME=trca_ph_mb_1_ac_

awk -F , '{print $NF}' volume_data/${VOLUME}/job_*/summary_all_runs.csv | sort -u | nl -w 4 -n rz | head -${DESIGNS} | sed "s/^/${NAME}/" >> orders_instance_${VOLUME}_dec2025.csv

