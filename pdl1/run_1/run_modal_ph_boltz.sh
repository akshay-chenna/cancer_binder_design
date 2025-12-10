VOLUME_NAME=pdl1_run1

VOLUME_NAME=${VOLUME_NAME} TIMEOUT_MINUTES=1400 GPU=H100 uvx modal run modal_proteinhunter.py \
	--protein-seqs FTVEVDSLSHVAEFYGDVTMGCRFQPGSWDPNLSVIWQRVQPLPDVEVYRLDNGQENLTSQNFQYRGRARLVSEELTNGWAKLHVSRLRINDSGVYRCLVEMGGADYKQTTLTVKATYKTIIKSMQRRGGGEVELACESEGYPLATINWRDKSLRNIKSNDTVVKTPNQLFHVTSKITVKYSEKNNYTCAFVEKGEAPKGPSARFDIP \
	--protein-ids B \
	--protein-msas "" \
	--template-path fold_pdl1_uniprot_model_0.cif \
	--template-cif-chain-id A \
	--binder-chain A \
	--min-design-protein-length 50 \
	--max-design-protein-length 150 \
	--high-iptm-threshold 0.7 \
	--contact-residues "8,42,43,44,62,63,64,99,107,142,143,144,198,163,164,167,171" \
	--use-potentials \
	--num-designs 1000 \
	--num-cycles 20 \
	--num-gpus 10 \
	--name ${VOLUME_NAME}

uvx modal volume get $VOLUME_NAME / volume_data --force
