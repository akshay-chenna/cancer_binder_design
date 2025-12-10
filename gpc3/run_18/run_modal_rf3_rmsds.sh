VOLUME_NAME=gpc3_run18 TIMEOUT_HOURS=20 uvx modal run modal_rf3_rmsds.py \
    	--n-mpnn 8 \
    	--fixed-chains "B" \
    	--designable-chains "A" \
    	--num-designs 5000 \
    	--num-gpus 100
