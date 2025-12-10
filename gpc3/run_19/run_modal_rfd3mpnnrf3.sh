VOLUME_NAME=gpc3_run19
VOLUME_NAME=$VOLUME_NAME TIMEOUT_HOURS=1400 GPU=H100 uvx modal run modal_rfd3_volume.py \
	--pdb-in correct_uniprot_gpc3_59-477_renum.pdb \
    	--contig "50-150,/0,A1-419" \
   	--hotspots "A232,A242,A243,A404,A405" \
    	--n-mpnn 8 \
    	--fixed-chains "B" \
    	--designable-chains "A" \
    	--num-designs 2000 \
    	--num-gpus 40

uvx modal volume get $VOLUME_NAME / volume_data --force
cd volume_data && cat job_*/*csv | sort -ru  >> ../all_metrics.csv && cd -
awk -F , '$((NF-6)) > 0.5' all_metrics.csv | awk -F , '{print $NF, $((NF-7)), $((NF-6)), $((NF-2))}' | sort -rnk 3 > top_designs.txt

awk '{print $4}' top_designs.txt | sort -u | nl -w 4 -n rz | sed 's/^/gpc3_rfd3_mb_h2_ac_/' >> orders_instance_gpc3_run_19_dec2025.csv
