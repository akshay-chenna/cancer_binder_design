VOLUME_NAME=gpc3_test
VOLUME_NAME=$VOLUME_NAME TIMEOUT_MINUTES=1400 GPU=A100 uvx modal run modal_rfd3_volume.py \
	--pdb-in correct_uniprot_gpc3_59-477_renum.pdb \
    	--contig "50-150,/0,A1-419" \
   	--hotspots "A67,A71,A147,A153,A156,A275" \
    	--n-mpnn 8 \
    	--fixed-chains "B" \
	--designable-chains "A" \
    	--num-designs 4 \
    	--num-gpus 2
#uvx modal volume get $VOLUME_NAME / volume_data --force
