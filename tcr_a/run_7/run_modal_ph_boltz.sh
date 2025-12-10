VOLUME_NAME=trca_run7
VOLUME_NAME=$VOLUME_NAME TIMEOUT_MINUTES=1400 GPU=H100 uvx modal run modal_proteinhunter.py \
	--protein-ids C \
	--protein-msas "" \
	--name trca_run7 \
	--percent-x 30 \
	--min-design-protein-length 50 \
	--max-design-protein-length 150 \
	--high-iptm-threshold 0.7 \
	--contact-residues 31,33,70,71,75,76,81,86,87 \
	--use-potentials \
	--num-designs 2000 \
	--num-cycles 20 \
	--num-gpus 50 \
	--protein-seqs IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCD
#	--template-path ./gpc3_af3.cif \
#	--template-chain-id A \
#	--protein-ids A \

#uvx modal volume get $VOLUME_NAME / volume_data --force
