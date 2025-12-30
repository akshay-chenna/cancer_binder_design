VOLUME_NAME=gpc3_test

uvx modal volume delete $VOLUME_NAME -y
uvx modal volume create --version=2 $VOLUME_NAME

VOLUME_NAME=$VOLUME_NAME TIMEOUT_MINUTES=1400 GPU=L40S uvx --with numpy modal run modal_rfd1-solmpnn-af2.py \
	--pdb-in reordered_modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb \
   	--hotspot-res "B36,B48,B59,B81,B86,B89,B99" \
	--partialt 5 \
	--noise-scale-ca 1 \
	--noise-scale-frame 1 \
	--guiding-potentials "type:binder_ROG,weight:1,min_dist:5" \
	--guide-scale 2 \
	--guide-decay "quadratic" \
	--n-mpnn 2 \
	--mpnn-t 0.1 \
	--fixed-chains "B" \
	--designable-chains "A" \
	--no-use-binder-template \
    	--num-designs 4 \
    	--num-gpus 2
rm -rf volume_data
mkdir volume_data
uvx modal volume get $VOLUME_NAME / volume_data --force
