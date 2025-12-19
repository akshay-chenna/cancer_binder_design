for i in job_*/high_iptm_pdb/*pdb ; do python find_interface.py $i >> hotspot_fraction.txt ; done
awk '$2 > 0' hotspot_fraction.txt | sort -rnk 2
