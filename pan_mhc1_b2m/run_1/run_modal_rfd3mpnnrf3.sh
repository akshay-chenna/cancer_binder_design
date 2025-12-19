VOLUME_NAME=panmhcb2m_run1
VOLUME_NAME=$VOLUME_NAME TIMEOUT_HOURS=24 GPU=L40S uvx modal run modal_rfd3_volume.py \
	--pdb-in 7uc5_cleaned_renum.pdb \
    	--contig "60-200,/0,A1-277,/0,B1-99" \
   	--hotspots "B46,B85,B87,B92,B93" \
    	--n-mpnn 6 \
    	--fixed-chains "B,C" \
    	--designable-chains "A" \
    	--num-designs 1650 \
    	--num-gpus 33

mkdir volume_data && uvx modal volume get $VOLUME_NAME / volume_data --force

cd volume_data && cat job_*/*csv | sort -ru  >> ../all_metrics.csv && cd -
awk -F , '{print $((NF-7)), $((NF-3)), $((NF-2)), $((NF))}' all_metrics.csv > iptm_vs_targetalignedrmsd.txt
#awk -F , '$((NF-6)) > 0.5' all_metrics.csv | awk -F , '{print $NF, $((NF-7)), $((NF-6)), $((NF-2))}' | sort -rnk 3 > top_designs.txt
#awk '{print $4}' top_designs.txt | sort -u | nl -w 4 -n rz | sed 's/^/gpc3_rfd3_smb_h1_ac_/' >> orders_instance_gpc3_run_21_dec2025.csv

