VOLUME_NAME=gpc3_run21
mkdir $VOLUME_NAME
VOLUME_NAME=$VOLUME_NAME TIMEOUT_HOURS=24 GPU=H100 uvx modal run modal_rfd3_volume.py \
	--pdb-in correct_uniprot_gpc3_59-477_renum.pdb \
    	--contig "80-200,/0,A1-419" \
   	--hotspots "A67,A71,A147,A153,A156,A275" \
    	--n-mpnn 6 \
    	--fixed-chains "B" \
    	--designable-chains "A" \
    	--num-designs 10000 \
    	--num-gpus 50

uvx modal volume get $VOLUME_NAME / volume_data --force

cd volume_data && cat job_*/*csv | sort -ru  >> ../all_metrics.csv && cd -
awk -F , '$((NF-6)) > 0.5' all_metrics.csv | awk -F , '{print $NF, $((NF-7)), $((NF-6)), $((NF-2))}' | sort -rnk 3 > top_designs.txt
awk '{print $4}' top_designs.txt | sort -u | nl -w 4 -n rz | sed 's/^/gpc3_rfd3_smb_h1_ac_/' >> orders_instance_gpc3_run_21_dec2025.csv
