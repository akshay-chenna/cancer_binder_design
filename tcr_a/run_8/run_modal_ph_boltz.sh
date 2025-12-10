VOLUME_NAME=trca_run8
VOLUME_NAME=$VOLUME_NAME TIMEOUT_MINUTES=1400 GPU=H100 uvx modal run modal_proteinhunter.py \
	--protein-ids C \
	--protein-msas "" \
	--name $VOLUME_NAME \
	--percent-x 30 \
	--min-design-protein-length 50 \
	--max-design-protein-length 150 \
	--high-iptm-threshold 0.7 \
	--contact-residues 139,141,178,179,183,184,189,194,195 \
	--use-potentials \
	--num-designs 2000 \
	--num-cycles 20 \
	--num-gpus 50 \
	--protein-seqs VEQDPGPLSVPEGAIVSLNCTYSNSAFQYFMWYRQYSRKGPELLMYTYSSGNKEDGRFTAQVDKSSKYISLFIRDSQPSDSATYLCAMSKGYSTLTFGKGTMLLVSPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCD
#	--template-path ./gpc3_af3.cif \
#	--template-chain-id A \
#	--protein-ids A \

uvx modal volume get $VOLUME_NAME / volume_data --force
